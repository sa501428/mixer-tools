/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2011-2021 Rice University, Baylor College of Medicine, Aiden Lab
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 *  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 *  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 *  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 *  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 *  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 *  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 *  THE SOFTWARE.
 */

package mixer.utils.slice.cleaning;

import javastraw.reader.Dataset;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.type.NormalizationType;
import javastraw.tools.ExtractingOEDataUtils;
import javastraw.tools.HiCFileTools;
import mixer.MixerGlobals;
import mixer.utils.similaritymeasures.RobustCorrelationSimilarity;
import mixer.utils.similaritymeasures.SimilarityMetric;

import java.util.*;

public class IndexOrderer {

    private final Map<Chromosome, int[]> chromToReorderedIndices = new HashMap<>();
    private final int DISTANCE = 5000000;
    private final int resolution;
    private final int minDistanceThreshold;
    private final int numColsToJoin;
    private final int IGNORE = -1;
    private final int DEFAULT = -5;
    private final int CHECK_VAL = -2;
    private final float CORR_MIN = 0.2f;
    private final float INCREMENT = .1f;
    private final Random generator = new Random(0);
    private final Map<Integer, Integer> indexToRearrangedLength = new HashMap<>();
    private final Map<Integer, int[]> indexToWeights = new HashMap<>();
    private final int[] weights;

    public IndexOrderer(Dataset ds, Chromosome[] chromosomes, int resolution, NormalizationType normalizationType,
                        int numColumnsToPutTogether, GWBadIndexFinder badIndexLocations, long seed) {
        this.resolution = resolution;
        minDistanceThreshold = DISTANCE / resolution;
        numColsToJoin = numColumnsToPutTogether;
        generator.setSeed(seed);
        for (Chromosome chrom : chromosomes) {
            final MatrixZoomData zd = HiCFileTools.getMatrixZoomData(ds, chrom, chrom, resolution);
            try {
                float[][] matrix = HiCFileTools.getOEMatrixForChromosome(ds, zd, chrom, resolution,
                        normalizationType, 3f, ExtractingOEDataUtils.ThresholdType.TRUE_OE,
                        true, 1, 0);
                Set<Integer> badIndices = badIndexLocations.getBadIndices(chrom);

                matrix = IntraMatrixCleaner.clean(chrom, matrix, resolution, numColsToJoin, badIndices);
                int[] newOrderIndexes = getNewOrderOfIndices(chrom, matrix, badIndices);
                chromToReorderedIndices.put(chrom, newOrderIndexes);
            } catch (Exception e) {
                e.printStackTrace();
            }
            System.out.print(".");
        }
        weights = appendWeights(chromosomes);
    }

    private int[] appendWeights(Chromosome[] chromosomes) {
        int totalLength = 0;
        for (Chromosome chromosome : chromosomes) {
            totalLength += indexToWeights.get(chromosome.getIndex()).length;
        }
        int[] finalWeights = new int[totalLength];
        int offset = 0;
        for (Chromosome chromosome : chromosomes) {
            int[] region = indexToWeights.get(chromosome.getIndex());
            System.arraycopy(region, 0, finalWeights, offset, region.length);
            offset += region.length;
        }
        return finalWeights;
    }

    public Map<Integer, Integer> getIndexToRearrangedLength() {
        return indexToRearrangedLength;
    }

    public int[] get(Chromosome chrom) {
        return chromToReorderedIndices.get(chrom);
    }

    private int[] getNewOrderOfIndices(Chromosome chromosome, float[][] matrix, Set<Integer> badIndices) {
        int[] newIndexOrderAssignments = generateNewAssignments(matrix.length, badIndices);

        int gCounter = doFirstRoundOfAssignmentsByCentroids(matrix, newIndexOrderAssignments, chromosome.getName());
        gCounter = doSecondRoundOfAssignments(matrix, newIndexOrderAssignments, gCounter);
        indexToRearrangedLength.put(chromosome.getIndex(), gCounter);
        indexToWeights.put(chromosome.getIndex(), generateWeights(gCounter, newIndexOrderAssignments));
        return newIndexOrderAssignments;
    }

    private int[] generateNewAssignments(int length, Set<Integer> badIndices) {
        int[] newIndexOrderAssignments = new int[length];
        Arrays.fill(newIndexOrderAssignments, DEFAULT);
        for (int k : badIndices) {
            newIndexOrderAssignments[k] = IGNORE;
        }
        return newIndexOrderAssignments;
    }

    private int[] generateWeights(int maxlength, int[] assignments) {
        int length = (int) Math.ceil((double) maxlength / numColsToJoin);
        int[] weightsForChrom = new int[length];
        for (int i : assignments) {
            if (i > -1) {
                weightsForChrom[i / numColsToJoin]++;
            }
        }
        return weightsForChrom;
    }

    private int doFirstRoundOfAssignmentsByCentroids(float[][] matrix, int[] newIndexOrderAssignments, String chromName) {
        int numInitialCentroids = 10;
        float[][] centroids = new QuickCentroids(quickCleanMatrix(matrix, newIndexOrderAssignments), numInitialCentroids, generator.nextLong()).generateCentroids();

        if (MixerGlobals.printVerboseComments) {
            System.out.println("IndexOrderer: Planned centroids for " + chromName + ": " + numInitialCentroids + " Actual centroids: " + centroids.length);
        }
        SimilarityMetric corrMetric = RobustCorrelationSimilarity.SINGLETON;

        int vectorLength = newIndexOrderAssignments.length;
        int[] numDecentRelations = new int[centroids.length];
        float[][] correlationCentroidsWithData = new float[centroids.length][vectorLength];
        for (int k = 0; k < centroids.length; k++) {
            for (int z = 0; z < vectorLength; z++) {
                if (newIndexOrderAssignments[z] < CHECK_VAL) {
                    float corr = corrMetric.distance(centroids[k], matrix[z]);
                    correlationCentroidsWithData[k][z] = corr;
                    if (corr > CORR_MIN || corr < -CORR_MIN) {
                        numDecentRelations[k]++;
                    }
                } else {
                    correlationCentroidsWithData[k][z] = Float.NaN;
                }
            }
        }

        int maxIndex = 0;
        for (int k = 1; k < numDecentRelations.length; k++) {
            if (numDecentRelations[maxIndex] < numDecentRelations[k]) {
                maxIndex = k;
            }
        }

        int gCounter = doSequentialOrdering(correlationCentroidsWithData[maxIndex], newIndexOrderAssignments, 0);
        for (int c = 0; c < centroids.length; c++) {
            if (c == maxIndex) continue;
            gCounter = doSequentialOrdering(correlationCentroidsWithData[c],
                    newIndexOrderAssignments, gCounter);
        }

        return gCounter;
    }

    private float[][] quickCleanMatrix(float[][] matrix, int[] newIndexOrderAssignments) {
        List<Integer> actualIndices = new ArrayList<>();
        for (int z = 0; z < newIndexOrderAssignments.length; z++) {
            if (newIndexOrderAssignments[z] < CHECK_VAL) {
                actualIndices.add(z);
            }
        }

        float[][] tempCleanMatrix = new float[actualIndices.size()][matrix[0].length];
        for (int i = 0; i < actualIndices.size(); i++) {
            System.arraycopy(matrix[actualIndices.get(i)], 0, tempCleanMatrix[i], 0, tempCleanMatrix[i].length);
        }
        if (MixerGlobals.printVerboseComments) {
            System.out.println("New clean matrix: " + tempCleanMatrix.length + " rows kept from " + matrix.length);
        }
        return tempCleanMatrix;
    }

    private int doSequentialOrdering(float[] correlationWithCentroid, int[] newIndexOrderAssignments, int startCounter) {
        int counter = startCounter;
        for (float cutoff = 1 - INCREMENT; cutoff >= CORR_MIN; cutoff -= INCREMENT) {
            for (int z = 0; z < correlationWithCentroid.length; z++) {
                if (newIndexOrderAssignments[z] < CHECK_VAL && correlationWithCentroid[z] > cutoff) {
                    newIndexOrderAssignments[z] = counter++;
                }
            }
        }

        counter = getUpdatedNoMixIndex(counter);

        for (float cutoff = CORR_MIN; cutoff < 1; cutoff += INCREMENT) {
            for (int z = 0; z < correlationWithCentroid.length; z++) {
                float corr = correlationWithCentroid[z];
                float cutoff1 = -cutoff;
                float cutoff2 = cutoff1 - INCREMENT;
                if (newIndexOrderAssignments[z] < CHECK_VAL && corr < cutoff1 && corr >= cutoff2) {
                    newIndexOrderAssignments[z] = counter++;
                }
            }
        }

        return getUpdatedNoMixIndex(counter);
    }

    private int getUpdatedNoMixIndex(int counter) {
        int temp = (counter / numColsToJoin);
        if (counter % numColsToJoin > 0) {
            temp++;
        }
        return temp * numColsToJoin;
    }

    private int doSecondRoundOfAssignments(float[][] matrix, int[] newIndexOrderAssignments, int startCounter) {
        int vectorLength = newIndexOrderAssignments.length;
        int numRoundsThatHappen = 0;
        int counter = startCounter;
        SimilarityMetric corrMetric = RobustCorrelationSimilarity.SINGLETON;
        for (int cI = 0; cI < vectorLength; cI++) {
            // handle stuff
            if (newIndexOrderAssignments[cI] < CHECK_VAL) {
                numRoundsThatHappen++;
                newIndexOrderAssignments[cI] = counter++;

                for (int z = cI + 1; z < vectorLength; z++) {
                    if (newIndexOrderAssignments[z] < CHECK_VAL) {
                        float val = corrMetric.distance(matrix[cI], matrix[z]);
                        if (val >= CORR_MIN) {
                            newIndexOrderAssignments[z] = counter++;
                        }
                    }
                }
                counter = getUpdatedNoMixIndex(counter);
            }
            // else it has already been handled
            // or is a bad index, so skip
        }
        if (MixerGlobals.printVerboseComments) {
            System.out.println("Num rounds " + numRoundsThatHappen);
        }
        return counter;
    }

    public int[] getWeights() {
        return weights;
    }
}

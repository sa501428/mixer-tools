/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2011-2022 Rice University, Baylor College of Medicine, Aiden Lab
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

package mixer.utils.magic;

import javastraw.reader.Dataset;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.block.ContactRecord;
import javastraw.reader.mzd.Matrix;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.type.HiCZoom;
import javastraw.reader.type.NormalizationType;
import javastraw.tools.MatrixTools;
import mixer.utils.InterChromosomeRegion;
import mixer.utils.bed.BedFileMappings;
import mixer.utils.common.ArrayTools;
import mixer.utils.common.LogTools;
import mixer.utils.common.ZScoreTools;
import mixer.utils.drive.DriveMatrix;
import mixer.utils.magic.clean.EmptyRowCleaner;
import mixer.utils.slice.structures.SubcompartmentInterval;

import java.io.File;
import java.util.*;

public class MagicMatrix extends DriveMatrix {

    private final int numRows, numCols;
    private final float[][] matrix;
    private final int[] weights;
    private final Map<Integer, SubcompartmentInterval> rowIndexToIntervalMap;
    private final List<InterChromosomeRegion> regionsToIgnore;


    public MagicMatrix(Dataset ds, ChromosomeHandler chromosomeHandler, int resolution,
                       NormalizationType norm, File outputDirectory, long seed, BedFileMappings mappings,
                       List<InterChromosomeRegion> regionsToIgnore, boolean doScale, boolean useZscore) {
        numRows = mappings.getNumRows();
        numCols = mappings.getNumCols();

        weights = new int[numCols];
        this.regionsToIgnore = regionsToIgnore;
        this.rowIndexToIntervalMap = createIndexToIntervalMap(chromosomeHandler, resolution);
        matrix = populateMatrix(ds, chromosomeHandler, resolution,
                norm, mappings, doScale);
        System.out.println("MAGIC matrix loaded");

        if (useZscore) {
            ZScoreTools.inPlaceZscoreDownCol(matrix);
            System.out.println("MAGIC matrix zScored");
        }

        inPlaceScaleSqrtWeightCol();
        System.out.println("final magic matrix num rows: " + matrix.length);
    }

    private Map<Integer, SubcompartmentInterval> createIndexToIntervalMap(ChromosomeHandler handler,
                                                                          int resolution) {
        Map<Integer, SubcompartmentInterval> map = new HashMap<>();
        int counter = 0;
        Chromosome[] chromosomes = handler.getAutosomalChromosomesArray();
        for (Chromosome chromosome : chromosomes) {
            int maxGenomeLen = (int) chromosome.getLength();
            int chrBinLength = (int) (chromosome.getLength() / resolution + 1);
            for (int i = 0; i < chrBinLength; i++) {
                int x1 = i * resolution;
                int x2 = Math.min(x1 + resolution, maxGenomeLen);
                SubcompartmentInterval newRInterval = new SubcompartmentInterval(chromosome, x1, x2, counter);
                map.put(counter, newRInterval);
                counter++;
            }
        }
        return map;
    }

    private static void scaleMatrixColumns(float[][] matrix, int[] scalars) {
        for (int i = 0; i < matrix.length; i++) {
            for (int z = 0; z < matrix[i].length; z++) {
                matrix[i][z] *= scalars[z];
            }
        }
    }

    private static void normalizeMatrix(float[][] matrix, Map<Integer, int[]> chromIndexToNumTotalLoci, int[] binIndexToChromIndex) {
        for (int i = 0; i < matrix.length; i++) {
            if (binIndexToChromIndex[i] > -1) {
                divide(matrix[i], chromIndexToNumTotalLoci.get(binIndexToChromIndex[i]));
            }
        }
    }

    private int[] getSumOfAllLoci(Map<Integer, int[]> chromIndexToNumLoci) {
        int[] totalLoci = new int[numCols];
        for (int[] row : chromIndexToNumLoci.values()) {
            for (int z = 0; z < row.length; z++) {
                totalLoci[z] += row[z];
            }
        }
        return totalLoci;
    }

    private static void divide(float[] row, int[] totalLoci) {
        for (int k = 0; k < row.length; k++) {
            if (totalLoci[k] > 0) {
                row[k] = row[k] / totalLoci[k];
            } else if (row[k] > 0) {
                System.err.println("Impossible situation reached: row val: " + row[k] + " but expect no entries: " + totalLoci[k]);
            }
        }
    }

    private static void updateNumberOfLoci(int[] totalLoci, int[] lociForRegion) {
        for (int i = 0; i < totalLoci.length; i++) {
            totalLoci[i] += lociForRegion[i];
        }
    }

    public void export(String path) {
        MatrixTools.saveMatrixTextNumpy(path, matrix);
    }

    @Override
    public float[][] getData(boolean getCorrelationMatrix) {
        return matrix;
    }

    @Override
    public Map<Integer, SubcompartmentInterval> getRowIndexToIntervalMap() {
        return rowIndexToIntervalMap;
    }

    @Override
    public void inPlaceScaleSqrtWeightCol() {
        ZScoreTools.inPlaceScaleSqrtWeightCol(matrix, weights);
    }

    private static void populateMatrixFromIterator(float[][] matrix, Iterator<ContactRecord> iterator, int rowOffset, int colOffset,
                                                   int[] binToClusterID) {
        while (iterator.hasNext()) {
            ContactRecord cr = iterator.next();
            if (cr.getCounts() > 0) {
                int r = cr.getBinX() + rowOffset;
                int c = cr.getBinY() + colOffset;
                if (binToClusterID[c] > -1 && binToClusterID[r] > -1) {
                    matrix[r][binToClusterID[c]] += cr.getCounts();
                    matrix[c][binToClusterID[r]] += cr.getCounts();
                }
            }
        }
    }

    private boolean shouldSkipRegion(Chromosome c1, Chromosome c2) {
        for (InterChromosomeRegion region : regionsToIgnore) {
            if (region.is(c1, c2)) {
                return true;
            }
        }
        return false;
    }

    private float[][] populateMatrix(Dataset ds, ChromosomeHandler handler, int resolution,
                                     NormalizationType norm, BedFileMappings mappings,
                                     boolean doScale) {
        float[][] matrix = new float[numRows][numCols];
        int[] offsets = mappings.getOffsets();
        int[] indexToClusterID = mappings.getIndexToClusterID();
        Map<Integer, int[]> distributionForChrom = mappings.getChromIndexToDistributionForChromosome();

        // incorporate norm
        Chromosome[] chromosomes = handler.getAutosomalChromosomesArray();
        Map<Integer, int[]> genomewideDistributionForChrom = new HashMap<>();
        for (Chromosome chromosome : chromosomes) {
            genomewideDistributionForChrom.put(chromosome.getIndex(), new int[numCols]);
        }

        for (int i = 0; i < chromosomes.length; i++) {
            for (int j = i + 1; j < chromosomes.length; j++) {
                if (shouldSkipRegion(chromosomes[i], chromosomes[j])) continue;

                Matrix m1 = ds.getMatrix(chromosomes[i], chromosomes[j]);
                if (m1 == null) continue;
                MatrixZoomData zd = m1.getZoomData(new HiCZoom(resolution));
                if (zd == null) continue;

                if (norm.getLabel().equalsIgnoreCase("none")) {
                    populateMatrixFromIterator(matrix, zd.getDirectIterator(), offsets[i], offsets[j], indexToClusterID);
                } else {
                    populateMatrixFromIterator(matrix, zd.getNormalizedIterator(norm), offsets[i], offsets[j], indexToClusterID);
                }

                updateNumberOfLoci(genomewideDistributionForChrom.get(chromosomes[i].getIndex()),
                        distributionForChrom.get(chromosomes[j].getIndex()));
                updateNumberOfLoci(genomewideDistributionForChrom.get(chromosomes[j].getIndex()),
                        distributionForChrom.get(chromosomes[i].getIndex()));

                System.out.print(".");
            }
            System.out.println(".");
        }

        LogTools.simpleLogWithCleanup(matrix, Float.NaN);


        normalizeMatrix(matrix, genomewideDistributionForChrom, mappings.getBinIndexToChromIndex());

        float[] coverage = new float[numRows];
        updateCoverage(matrix, coverage);
        scaleCoverage(coverage);

        int[] totalDistribution = getSumOfAllLoci(distributionForChrom);
        System.arraycopy(totalDistribution, 0, weights, 0, Math.max(weights.length, totalDistribution.length));
        scaleMatrixColumns(matrix, totalDistribution);

        matrix = EmptyRowCleaner.cleanUpMatrix(matrix, rowIndexToIntervalMap, coverage);

        LogTools.simpleExpm1(matrix);

        if (doScale) {
            return FinalScale.scaleMatrix(new SymmLLInterMatrix(matrix), createTargetVector(totalDistribution));
        } else {
            return matrix;
        }
    }

    private int[] createTargetVector(int[] colSums) {
        int[] target = new int[numRows + numCols];
        Arrays.fill(target, 1);
        for (int i = 0; i < numCols; i++) {
            target[i] = colSums[i];
        }
        return target;
    }

    private void scaleCoverage(float[] coverage) {
        float mean = ArrayTools.getNonZeroMean(coverage);
        for (int k = 0; k < coverage.length; k++) {
            coverage[k] /= mean;
        }
    }

    private void updateCoverage(float[][] matrix, float[] coverage) {
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                coverage[i] += matrix[i][j];
            }
        }
    }
}

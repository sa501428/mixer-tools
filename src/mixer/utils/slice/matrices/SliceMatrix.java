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

package mixer.utils.slice.matrices;

import javastraw.reader.Dataset;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.block.ContactRecord;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.type.NormalizationType;
import javastraw.tools.HiCFileTools;
import mixer.MixerGlobals;
import mixer.algos.Slice;
import mixer.utils.slice.cleaning.BadIndexFinder;
import mixer.utils.slice.cleaning.IndexOrderer;
import mixer.utils.slice.structures.SubcompartmentInterval;

import java.io.File;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;

public class SliceMatrix extends CompositeGenomeWideMatrix {

    public SliceMatrix(ChromosomeHandler chromosomeHandler, Dataset ds, NormalizationType[] norms,
                       int resolution, File outputDirectory, long seed, BadIndexFinder badIndexLocations,
                       int maxClusterSizeExpected) {
        super(chromosomeHandler, ds, norms, resolution, outputDirectory, seed,
                badIndexLocations, maxClusterSizeExpected);
    }

    MatrixAndWeight makeCleanScaledInterMatrix(Dataset ds, NormalizationType interNorm) {
        // height/width chromosomes
        Map<Integer, Integer> indexToLength = calculateActualLengthForChromosomes(chromosomes);
        IndexOrderer orderer = new IndexOrderer(ds, chromosomes, resolution, norms[Slice.INTRA_SCALE_INDEX], badIndexLocations,
                generator.nextLong(), outputDirectory, maxClusterSizeExpected);
        Map<Integer, Integer> indexToCompressedLength = calculateCompressedLengthForChromosomes(orderer.getIndexToRearrangedLength());

        Dimension dimensions = new Dimension(chromosomes, indexToLength);
        Dimension compressedDimensions = new Dimension(chromosomes, indexToCompressedLength);

        int[] weights = getWeights(compressedDimensions, orderer);

        if (MixerGlobals.printVerboseComments) {
            System.out.println(dimensions.length + " by " + compressedDimensions.length);
        }

        if (dimensions.length == 0 || compressedDimensions.length == 0) {
            System.err.println("No matrix created; map is likely too sparse\n" +
                    "Try a lower resolution or higher compression");
            System.exit(9);
        }

        System.out.println("Indexing complete");
        float[][] interMatrix = new float[dimensions.length][compressedDimensions.length];

        for (int i = 0; i < chromosomes.length; i++) {
            Chromosome chr1 = chromosomes[i];

            for (int j = i; j < chromosomes.length; j++) {
                Chromosome chr2 = chromosomes[j];

                fillInChromosomeRegion(interMatrix, badIndexLocations, ds, resolution, i == j, orderer,
                        chr1, dimensions.offset[i], compressedDimensions.offset[i],
                        chr2, dimensions.offset[j], compressedDimensions.offset[j],
                        interNorm);
                System.out.print(".");
            }
        }
        System.out.println(".");
        System.out.println("Inter Matrix Built");

        return new MatrixAndWeight(interMatrix, weights);
    }

    private int[] getWeights(Dimension compressedDimensions, IndexOrderer orderer) {
        int[] weights = new int[compressedDimensions.length];
        for (int i = 0; i < chromosomes.length; i++) {
            Chromosome chr1 = chromosomes[i];
            Map<Integer, Integer> colPosChrom1 = makeLocalReorderedIndexMap(chr1,
                    badIndexLocations.getBadIndices(chr1), compressedDimensions.offset[i], orderer.get(chr1));
            updateWeights(weights, colPosChrom1);
        }
        return weights;
    }

    private void updateWeights(int[] weights, Map<Integer, Integer> colPosChrom) {
        for (int i : colPosChrom.values()) {
            weights[i]++;
        }
    }

    protected Map<Integer, Integer> calculateActualLengthForChromosomes(Chromosome[] chromosomes) {
        Map<Integer, Integer> indexToFilteredLength = new HashMap<>();
        for (Chromosome chrom : chromosomes) {
            indexToFilteredLength.put(chrom.getIndex(), (int) Math.ceil((float) chrom.getLength() / resolution) - badIndexLocations.getBadIndices(chrom).size());
        }
        return indexToFilteredLength;
    }

    private Map<Integer, Integer> calculateCompressedLengthForChromosomes(Map<Integer, Integer> initialMap) {
        Map<Integer, Integer> indexToCompressedLength = new HashMap<>();
        for (Integer key : initialMap.keySet()) {
            int val = (int) Math.ceil((double) initialMap.get(key));
            //System.out.println("size of " + key + " " + val + " was (" + initialMap.get(key) + ") num cols " + numColumnsToPutTogether);
            indexToCompressedLength.put(key, val);
        }

        return indexToCompressedLength;
    }

    private void fillInChromosomeRegion(float[][] matrix, BadIndexFinder badIndices,
                                        Dataset ds, int resolution, boolean isIntra, IndexOrderer orderer,
                                        Chromosome chr1, int offsetIndex1, int compressedOffsetIndex1,
                                        Chromosome chr2, int offsetIndex2, int compressedOffsetIndex2,
                                        NormalizationType interNorm) {
        Iterator<ContactRecord> recordIterator = null;
        try {
            if (!isIntra) {
                final MatrixZoomData zd = HiCFileTools.getMatrixZoomData(ds, chr1, chr2, resolution);
                recordIterator = zd.getNormalizedIterator(interNorm);
            }
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(99);
        }

        Map<Integer, Integer> rowPosChrom1 = makeLocalIndexMap(chr1, badIndices.getBadIndices(chr1), offsetIndex1);
        Map<Integer, Integer> rowPosChrom2 = makeLocalIndexMap(chr2, badIndices.getBadIndices(chr2), offsetIndex2);

        Map<Integer, Integer> colPosChrom1 = makeLocalReorderedIndexMap(chr1,
                badIndices.getBadIndices(chr1), compressedOffsetIndex1, orderer.get(chr1));
        Map<Integer, Integer> colPosChrom2 = makeLocalReorderedIndexMap(chr2,
                badIndices.getBadIndices(chr2), compressedOffsetIndex2, orderer.get(chr2));

        if (isIntra) {
            updateSubcompartmentMap(chr1, badIndices.getBadIndices(chr1), offsetIndex1, rowIndexToIntervalMap);
        }

        copyValuesToArea(matrix, recordIterator,
                rowPosChrom1, colPosChrom1, rowPosChrom2, colPosChrom2, isIntra);
    }


    private void copyValuesToArea(float[][] matrix, Iterator<ContactRecord> iterator,
                                  Map<Integer, Integer> rowPosChrom1, Map<Integer, Integer> colPosChrom1,
                                  Map<Integer, Integer> rowPosChrom2, Map<Integer, Integer> colPosChrom2,
                                  boolean isIntra) {
        if (isIntra) {
            for (int binX : rowPosChrom1.keySet()) {
                for (int binY : colPosChrom2.keySet()) {
                    matrix[rowPosChrom1.get(binX)][colPosChrom2.get(binY)] = Float.NaN;
                }
            }
        } else {
            while (iterator.hasNext()) {
                ContactRecord cr = iterator.next();
                float val = cr.getCounts();
                if (!Float.isNaN(val)) {
                    int binX = cr.getBinX();
                    int binY = cr.getBinY();

                    if (rowPosChrom1.containsKey(binX) && colPosChrom2.containsKey(binY)) {
                        matrix[rowPosChrom1.get(binX)][colPosChrom2.get(binY)] += val;
                    }
                    if (rowPosChrom2.containsKey(binY) && colPosChrom1.containsKey(binX)) {
                        matrix[rowPosChrom2.get(binY)][colPosChrom1.get(binX)] += val;
                    }
                }
            }
        }
    }

    private void updateSubcompartmentMap(Chromosome chromosome, Set<Integer> badIndices, int offsetIndex1, Map<Integer, SubcompartmentInterval> rowIndexToIntervalMap) {
        int counter = 0;
        int maxGenomeLen = (int) chromosome.getLength();
        int chrBinLength = (int) (chromosome.getLength() / resolution + 1);
        for (int i = 0; i < chrBinLength; i++) {
            if (badIndices.contains(i)) {
                continue;
            }
            int x1 = i * resolution;
            int x2 = Math.min(x1 + resolution, maxGenomeLen);
            SubcompartmentInterval newRInterval = new SubcompartmentInterval(chromosome, x1, x2, counter);
            rowIndexToIntervalMap.put(offsetIndex1 + counter, newRInterval);
            counter++;
        }
    }

    private Map<Integer, Integer> makeLocalIndexMap(Chromosome chrom, Set<Integer> badIndices, int offsetIndex) {
        Map<Integer, Integer> binToLocalMap = new HashMap<>();
        int counter = 0;

        int chrLength = (int) (chrom.getLength() / resolution + 1);
        for (int i = 0; i < chrLength; i++) {
            if (badIndices.contains(i)) {
                continue;
            }

            binToLocalMap.put(i, offsetIndex + (counter));
            counter++;
        }

        return binToLocalMap;
    }

    private Map<Integer, Integer> makeLocalReorderedIndexMap(Chromosome chrom, Set<Integer> badIndices,
                                                             int offsetIndex, int[] newOrder) {
        Map<Integer, Integer> binToLocalMap = new HashMap<>();

        int chrLength = (int) (chrom.getLength() / resolution + 1);
        for (int i = 0; i < chrLength; i++) {
            if (badIndices.contains(i) || newOrder[i] < 0) {
                continue;
            }

            binToLocalMap.put(i, offsetIndex + (newOrder[i]));
        }

        return binToLocalMap;
    }
}

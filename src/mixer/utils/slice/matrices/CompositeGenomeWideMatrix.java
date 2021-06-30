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

package mixer.utils.slice.matrices;

import javastraw.feature1D.GenomeWideList;
import javastraw.reader.Dataset;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.type.NormalizationType;
import javastraw.tools.MatrixTools;
import mixer.MixerGlobals;
import mixer.algos.Slice;
import mixer.utils.common.FloatMatrixTools;
import mixer.utils.similaritymeasures.RobustEuclideanDistance;
import mixer.utils.slice.cleaning.GWBadIndexFinder;
import mixer.utils.slice.cleaning.SimilarityMatrixTools;
import mixer.utils.slice.gmm.SimpleScatterPlot;
import mixer.utils.slice.kmeans.kmeansfloat.Cluster;
import mixer.utils.slice.structures.SliceUtils;
import mixer.utils.slice.structures.SubcompartmentInterval;

import java.io.File;
import java.util.*;

public abstract class CompositeGenomeWideMatrix {
    protected final NormalizationType[] norms;
    protected final int resolution;
    protected final Map<Integer, SubcompartmentInterval> rowIndexToIntervalMap = new HashMap<>();
    protected final Chromosome[] chromosomes;
    protected final Random generator = new Random(0);
    protected final File outputDirectory;
    private float[][] gwCleanMatrix, projectedData;
    protected final GWBadIndexFinder badIndexLocations;
    protected final int maxClusterSizeExpected;

    public CompositeGenomeWideMatrix(ChromosomeHandler chromosomeHandler, Dataset ds,
                                     NormalizationType[] norms,
                                     int resolution,
                                     File outputDirectory, long seed,
                                     GWBadIndexFinder badIndexLocations,
                                     int maxClusterSizeExpected) {
        this.maxClusterSizeExpected = maxClusterSizeExpected;
        this.norms = norms;
        System.out.println("Norms: " + Arrays.toString(norms));
        this.resolution = resolution;
        this.outputDirectory = outputDirectory;
        this.generator.setSeed(seed);
        this.badIndexLocations = badIndexLocations;

        chromosomes = chromosomeHandler.getAutosomalChromosomesArray();
        gwCleanMatrix = makeCleanScaledInterMatrix(ds, norms[Slice.INTER_SCALE_INDEX]);
        //gwCleanMatrix = makeComboNormMatrix(ds, norms);
    }

    private float[][] makeComboNormMatrix(Dataset ds, NormalizationType[] norms) {
        return FloatMatrixTools.concatenate(
                makeCleanScaledInterMatrix(ds, norms[Slice.INTER_SCALE_INDEX]),
                FloatMatrixTools.multiply(makeCleanScaledInterMatrix(ds, norms[Slice.GW_SCALE_INDEX]), 0.5f)
        );
    }

    abstract float[][] makeCleanScaledInterMatrix(Dataset ds, NormalizationType interNorm);


    public void cleanUpMatricesBySparsity() {

        projectedData = SimilarityMatrixTools.getCosinePearsonCorrMatrix(gwCleanMatrix, 50, generator.nextLong());
        File file1 = new File(outputDirectory, "zscore_matrix.npy");
        MatrixTools.saveMatrixTextNumpy(file1.getAbsolutePath(), gwCleanMatrix);

        File file2 = new File(outputDirectory, "corr_matrix.npy");
        MatrixTools.saveMatrixTextNumpy(file2.getAbsolutePath(), projectedData);

        runUmapAndSaveMatrices(projectedData, outputDirectory,
                rowIndexToIntervalMap);

        /*
        SliceMatrixCleaner matrixCleanupReduction = new SliceMatrixCleaner(gwCleanMatrix,
                generator.nextLong(), outputDirectory, metric);
        gwCleanMatrix = matrixCleanupReduction.getCleanFilteredZscoredMatrix(rowIndexToIntervalMap, weights);
        projectedData = MatrixImputer.imputeUntilNoNansOnlyNN(gwCleanMatrix);
        //corrData = SimilarityMatrixTools.getCosinePearsonCorrMatrix(gwCleanMatrix, 50, generator.nextLong());
        */
    }

    private void runUmapAndSaveMatrices(float[][] data, File outputDirectory,
                                        Map<Integer, SubcompartmentInterval> rowIndexToIntervalMap) {
        int[] indices = new int[data.length];
        for (int i = 0; i < indices.length; i++) {
            SubcompartmentInterval interval = rowIndexToIntervalMap.get(i);
            indices[i] = interval.getChrIndex();
        }

        SimpleScatterPlot plotter = new SimpleScatterPlot(data);
        File outfile = new File(outputDirectory, "umap");
        plotter.plot(indices, outfile.getAbsolutePath());
    }

    public synchronized double processKMeansClusteringResult(Cluster[] clusters,
                                                             GenomeWideList<SubcompartmentInterval> subcompartments,
                                                             boolean useCorr) {

        Set<SubcompartmentInterval> subcompartmentIntervals = new HashSet<>();
        if (MixerGlobals.printVerboseComments) {
            System.out.println("GW Composite data vs clustered into " + clusters.length + " clusters");
        }

        double withinClusterSumOfSquares = 0;
        int numGoodClusters = 0;
        int genomewideCompartmentID = 0;

        Set<Integer> badIndices = new HashSet<>();

        for (int z = 0; z < clusters.length; z++) {
            Cluster cluster = clusters[z];
            int currentClusterID = ++genomewideCompartmentID;

            if (MixerGlobals.printVerboseComments) {
                System.out.println("Size of cluster " + currentClusterID + " - " + cluster.getMemberIndexes().length);
            }

            if (cluster.getMemberIndexes().length < 5) {
                for (int badIndex : cluster.getMemberIndexes()) {
                    badIndices.add(badIndex);
                }
                withinClusterSumOfSquares += Float.MAX_VALUE;
            } else {
                numGoodClusters++;
            }

            for (int i : cluster.getMemberIndexes()) {
                if (useCorr) {
                    withinClusterSumOfSquares +=
                            RobustEuclideanDistance.getNonNanMeanSquaredError(cluster.getCenter(), projectedData[i]);
                } else {
                    withinClusterSumOfSquares +=
                            RobustEuclideanDistance.getNonNanMeanSquaredError(cluster.getCenter(), gwCleanMatrix[i]);
                }

                if (rowIndexToIntervalMap.containsKey(i)) {
                    SubcompartmentInterval interv = rowIndexToIntervalMap.get(i);
                    if (interv != null) {
                        subcompartmentIntervals.add(generateNewSubcompartment(interv, currentClusterID));
                    }
                }
            }
        }

        withinClusterSumOfSquares = withinClusterSumOfSquares / numGoodClusters;
        if (MixerGlobals.printVerboseComments) {
            System.out.println("Final WCSS " + withinClusterSumOfSquares);
        }

        subcompartments.addAll(new ArrayList<>(subcompartmentIntervals));
        SliceUtils.reSort(subcompartments);

        return withinClusterSumOfSquares;
    }

    public synchronized void processGMMClusteringResult(int[] clusterID,
                                                        GenomeWideList<SubcompartmentInterval> subcompartments) {

        Set<SubcompartmentInterval> subcompartmentIntervals = new HashSet<>();

        for (int i = 0; i < projectedData.length; i++) {
            if (rowIndexToIntervalMap.containsKey(i)) {
                SubcompartmentInterval interv = rowIndexToIntervalMap.get(i);
                if (interv != null) {
                    subcompartmentIntervals.add(generateNewSubcompartment(interv, clusterID[i]));
                }
            }
        }

        subcompartments.addAll(new ArrayList<>(subcompartmentIntervals));
        SliceUtils.reSort(subcompartments);
    }

    protected SubcompartmentInterval generateNewSubcompartment(SubcompartmentInterval interv, int currentClusterID) {
        SubcompartmentInterval newInterv = (SubcompartmentInterval) interv.deepClone();
        newInterv.setClusterID(currentClusterID);
        return newInterv;
    }

    public float[][] getData(boolean getCorrMatrix) {
        if (getCorrMatrix) {
            return projectedData;
        }
        return gwCleanMatrix;
    }

    public void appendDataAlongExistingRows(CompositeGenomeWideMatrix additionalData) {
        if (gwCleanMatrix.length != additionalData.gwCleanMatrix.length) {
            System.err.println("***************************************\n" +
                    "Dimension mismatch: " + gwCleanMatrix.length + " != " + additionalData.gwCleanMatrix.length);
        } else {
            gwCleanMatrix = FloatMatrixTools.concatenate(gwCleanMatrix, additionalData.gwCleanMatrix);
        }
    }

    public GWBadIndexFinder getBadIndices() {
        return badIndexLocations;
    }
}

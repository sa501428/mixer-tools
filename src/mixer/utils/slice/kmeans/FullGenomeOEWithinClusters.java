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

package mixer.utils.slice.kmeans;

import javastraw.feature1D.GenomeWideList;
import javastraw.reader.Dataset;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.type.NormalizationType;
import javastraw.tools.MatrixTools;
import mixer.MixerGlobals;
import mixer.utils.similaritymeasures.SimilarityMetric;
import mixer.utils.slice.cleaning.GWBadIndexFinder;
import mixer.utils.slice.gmm.GenomeWideGMMRunner;
import mixer.utils.slice.matrices.SliceMatrix;
import mixer.utils.slice.structures.SliceUtils;
import mixer.utils.slice.structures.SubcompartmentInterval;

import java.io.File;
import java.util.*;

public class FullGenomeOEWithinClusters {
    public static int startingClusterSizeK = 2;
    public static int numClusterSizeKValsUsed = 10;
    public static int numAttemptsForKMeans = 3;
    protected final File outputDirectory;
    private final ChromosomeHandler chromosomeHandler;
    private final SliceMatrix sliceMatrix;
    private final int maxIters = 200;

    private final Random generator = new Random(0);

    public FullGenomeOEWithinClusters(List<Dataset> datasets, ChromosomeHandler chromosomeHandler, int resolution,
                                      NormalizationType[] intraNorms, NormalizationType[] interNorms,
                                      File outputDirectory, long seed, SimilarityMetric metric) {
        this.chromosomeHandler = chromosomeHandler;
        this.outputDirectory = outputDirectory;
        generator.setSeed(seed);

        GWBadIndexFinder badIndexFinder = new GWBadIndexFinder(chromosomeHandler.getAutosomalChromosomesArray(),
                resolution, interNorms);
        badIndexFinder.createInternalBadList(datasets, chromosomeHandler.getAutosomalChromosomesArray());

        sliceMatrix = new SliceMatrix(
                chromosomeHandler, datasets.get(0), intraNorms[0], interNorms[0], resolution, outputDirectory,
                generator.nextLong(), badIndexFinder, metric);

        if (datasets.size() > 1) {
            sliceMatrix.normalizePerDataset();
        }

        for (int dI = 1; dI < datasets.size(); dI++) {
            SliceMatrix additionalData = new SliceMatrix(chromosomeHandler, datasets.get(dI),
                    intraNorms[dI], interNorms[dI], resolution, outputDirectory,
                    generator.nextLong(), badIndexFinder, metric);
            additionalData.normalizePerDataset();
            sliceMatrix.appendDataAlongExistingRows(additionalData);
            // todo append weights
        }

        sliceMatrix.cleanUpMatricesBySparsity();
    }

    public void extractFinalGWSubcompartments(String prefix) {

        Map<Integer, GenomeWideList<SubcompartmentInterval>> kmeansClustersToResults = new HashMap<>();
        Map<Integer, GenomeWideList<SubcompartmentInterval>> gmmClustersToResults = new HashMap<>();
        Map<Integer, List<List<Integer>>> kmeansIndicesMap = new HashMap<>();

        if (MixerGlobals.printVerboseComments) {
            sliceMatrix.exportData();
        }

        GenomeWideKmeansRunner kmeansRunner = new GenomeWideKmeansRunner(chromosomeHandler, sliceMatrix);
        double[][] iterToWcssAicBic = new double[4][numClusterSizeKValsUsed];
        Arrays.fill(iterToWcssAicBic[1], Double.MAX_VALUE);
        Arrays.fill(iterToWcssAicBic[2], Double.MAX_VALUE);
        Arrays.fill(iterToWcssAicBic[3], Double.MAX_VALUE);

        System.out.println("Genomewide clustering");
        for (int numMaxIters : new int[]{100}) {
            for (int z = 0; z < numClusterSizeKValsUsed; z++) {
                runRepeatedKMeansClusteringLoop(1, kmeansRunner, iterToWcssAicBic, z,
                        numMaxIters, kmeansClustersToResults, kmeansIndicesMap);
            }
        }

        for (int z = 0; z < numClusterSizeKValsUsed; z++) {
            runRepeatedKMeansClusteringLoop(numAttemptsForKMeans, kmeansRunner, iterToWcssAicBic, z,
                    maxIters, kmeansClustersToResults, kmeansIndicesMap);

            exportKMeansClusteringResults(z, iterToWcssAicBic, kmeansClustersToResults, prefix);

            runGMMClusteringLoop(z, 20, kmeansIndicesMap.get(z), gmmClustersToResults);
            exportGMMClusteringResults(z, gmmClustersToResults, prefix);
        }
        System.out.println(".");
    }

    private void exportKMeansClusteringResults(int z, double[][] iterToWcssAicBic,
                                               Map<Integer, GenomeWideList<SubcompartmentInterval>> numClustersToResults,
                                               String prefix) {
        int k = z + startingClusterSizeK;
        String outIterPath = new File(outputDirectory, "clusterSize_WCSS_AIC_BIC.npy").getAbsolutePath();
        MatrixTools.saveMatrixTextNumpy(outIterPath, iterToWcssAicBic);
        GenomeWideList<SubcompartmentInterval> gwList = numClustersToResults.get(k);
        SliceUtils.collapseGWList(gwList);
        File outBedFile = new File(outputDirectory, prefix + "_" + k + "_kmeans_clusters.bed");
        gwList.simpleExport(outBedFile);
    }

    private void exportGMMClusteringResults(int z, Map<Integer, GenomeWideList<SubcompartmentInterval>> numClustersToResults,
                                            String prefix) {
        int k = z + startingClusterSizeK;
        if (numClustersToResults.containsKey(k)) {
            GenomeWideList<SubcompartmentInterval> gwList = numClustersToResults.get(k);
            SliceUtils.collapseGWList(gwList);
            File outBedFile = new File(outputDirectory, prefix + "_" + k + "_gmm_clusters.bed");
            gwList.simpleExport(outBedFile);
        }
    }


    private void runRepeatedKMeansClusteringLoop(int attemptsForKMeans, GenomeWideKmeansRunner kmeansRunner,
                                                 double[][] iterToWcssAicBic, int z, int maxIters,
                                                 Map<Integer, GenomeWideList<SubcompartmentInterval>> numClustersToResults,
                                                 Map<Integer, List<List<Integer>>> indicesMap) {
        int numClusters = z + startingClusterSizeK;
        for (int p = 0; p < attemptsForKMeans; p++) {
            kmeansRunner.prepareForNewRun(numClusters);
            kmeansRunner.launchKmeansGWMatrix(generator.nextLong(), maxIters);

            int numActualClustersThisAttempt = kmeansRunner.getNumActualClusters();
            if (numActualClustersThisAttempt == numClusters) {
                double wcss = kmeansRunner.getWithinClusterSumOfSquares();
                if (wcss < iterToWcssAicBic[1][z]) {
                    setMseAicBicValues(z, iterToWcssAicBic, numClusters, wcss);
                    indicesMap.put(z, kmeansRunner.getIndicesMapCopy());
                    numClustersToResults.put(numClusters, kmeansRunner.getFinalCompartments());
                }
            }
            System.out.print(".");
        }
    }

    private void runGMMClusteringLoop(int z, int maxIters, List<List<Integer>> startingIndices,
                                      Map<Integer, GenomeWideList<SubcompartmentInterval>> numClustersToResults) {
        int numClusters = z + startingClusterSizeK;
        GenomeWideGMMRunner gmmRunner = new GenomeWideGMMRunner(chromosomeHandler, sliceMatrix);
        gmmRunner.launch(numClusters, maxIters, numClustersToResults, startingIndices);
        System.out.print("*");
    }

    private void setMseAicBicValues(int z, double[][] iterToWcssAicBic, int numClusters, double sumOfSquares) {
        iterToWcssAicBic[0][z] = numClusters;
        iterToWcssAicBic[1][z] = sumOfSquares;
        // AIC
        iterToWcssAicBic[2][z] = sumOfSquares + 2 * sliceMatrix.getNumColumns() * numClusters;
        // BIC .5*k*d*log(n)
        iterToWcssAicBic[3][z] = sumOfSquares + 0.5 * sliceMatrix.getNumColumns() * numClusters * Math.log(sliceMatrix.getNumRows());
    }
}

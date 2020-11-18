/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2011-2020 Rice University, Baylor College of Medicine, Aiden Lab
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

package mixer.utils.slice;

import javastraw.featurelist.GenomeWideList;
import javastraw.reader.ChromosomeHandler;
import javastraw.reader.Dataset;
import javastraw.type.NormalizationType;
import mixer.MixerGlobals;
import mixer.utils.common.DoubleMatrixTools;
import mixer.utils.common.IntMatrixTools;
import mixer.utils.similaritymeasures.SimilarityMetric;
import mixer.utils.slice.cleaning.BadIndexFinder;
import mixer.utils.slice.cleaning.LeftOverClusterIdentifier;
import mixer.utils.slice.kmeansfloat.Cluster;
import mixer.utils.slice.kmeansfloat.ClusterTools;
import mixer.utils.slice.matrices.CompositeGenomeWideDensityMatrix;
import mixer.utils.slice.matrices.SliceMatrix;
import mixer.utils.slice.structures.SliceUtils;
import mixer.utils.slice.structures.SubcompartmentInterval;

import java.io.File;
import java.util.*;

public class FullGenomeOEWithinClusters {
    public static int startingClusterSizeK = 5;
    public static int numClusterSizeKValsUsed = 2;//10
    public static int numAttemptsForKMeans = 10;
    private final List<Dataset> datasets;
    protected final File outputDirectory;
    private final ChromosomeHandler chromosomeHandler;
    private final int resolution;
    private final NormalizationType norm;
    private final CompositeGenomeWideDensityMatrix interMatrix;
    private final SimilarityMetric metric;

    public FullGenomeOEWithinClusters(List<Dataset> datasets, ChromosomeHandler chromosomeHandler, int resolution, NormalizationType norm,
                                      File outputDirectory, Random generator, String[] referenceBedFiles,
                                      SimilarityMetric metric) {
        this.chromosomeHandler = chromosomeHandler;
        this.resolution = resolution;
        this.norm = norm;
        this.outputDirectory = outputDirectory;
        this.datasets = datasets;
        this.metric = metric;

        BadIndexFinder badIndexFinder = new BadIndexFinder(datasets,
                chromosomeHandler.getAutosomalChromosomesArray(), resolution, norm);

        interMatrix = new SliceMatrix(
                chromosomeHandler, datasets.get(0), norm, resolution, outputDirectory, generator, referenceBedFiles,
                badIndexFinder, 0, metric);

        for (int dI = 1; dI < datasets.size(); dI++) {
            SliceMatrix additionalData = new SliceMatrix(chromosomeHandler, datasets.get(dI),
                    norm, resolution, outputDirectory, generator, new String[]{}, badIndexFinder, dI, metric);
            interMatrix.appendDataAlongExistingRows(additionalData);
        }
    }

    public void extractFinalGWSubcompartments(Random generator, List<String> inputHicFilePaths,
                                              String prefix, int index, boolean compareMaps) {

        Map<Integer, GenomeWideList<SubcompartmentInterval>> numItersToResults = new HashMap<>();

        if (MixerGlobals.printVerboseComments) {
            interMatrix.exportData();
        }

        GenomeWideKmeansRunner kmeansRunner = new GenomeWideKmeansRunner(chromosomeHandler, interMatrix);

        double[][] iterToWcssAicBic = new double[4][numClusterSizeKValsUsed];
        Arrays.fill(iterToWcssAicBic[1], Double.MAX_VALUE);
        Arrays.fill(iterToWcssAicBic[2], Double.MAX_VALUE);
        Arrays.fill(iterToWcssAicBic[3], Double.MAX_VALUE);

        System.out.println("Genomewide clustering");
        for (int z = 0; z < numClusterSizeKValsUsed; z++) {

            int k = z + startingClusterSizeK;
            Cluster[] bestClusters = null;
            int[] bestIDs = null;
            int[][] novelIDsForIndx = null;

            for (int p = 0; p < numAttemptsForKMeans; p++) {

                kmeansRunner.prepareForNewRun(k);
                kmeansRunner.launchKmeansGWMatrix(generator.nextLong());

                int numActualClustersThisAttempt = kmeansRunner.getNumActualClusters();
                double wcss = kmeansRunner.getWithinClusterSumOfSquares();

                if (wcss < iterToWcssAicBic[1][z]) {
                    setMseAicBicValues(z, iterToWcssAicBic, numActualClustersThisAttempt, wcss);
                    numItersToResults.put(k, kmeansRunner.getFinalCompartments());
                    bestClusters = kmeansRunner.getRecentClustersClone();
                    bestIDs = kmeansRunner.getRecentIDsClone();
                    novelIDsForIndx = kmeansRunner.getRecentIDsForIndex();
                }
                System.out.print(".");
            }

            ClusterTools.performStatisticalAnalysisBetweenClusters(outputDirectory, "final_gw_" + k, bestClusters, bestIDs, metric);
            IntMatrixTools.saveMatrixTextNumpy((new File(outputDirectory, "novel_ids_for_index_" + k + ".npy")).getAbsolutePath(), novelIDsForIndx);
        }
        System.out.println(".");

        if (!compareMaps) {
            System.out.println("Post processing");
            LeftOverClusterIdentifier identifier = new LeftOverClusterIdentifier(chromosomeHandler, datasets.get(0), norm, resolution);
            identifier.identify(numItersToResults, interMatrix.getBadIndices());
        }

        DoubleMatrixTools.saveMatrixTextNumpy(new File(outputDirectory, "clusterSize_WCSS_AIC_BIC.npy").getAbsolutePath(), iterToWcssAicBic);

        String hicFileName = SliceUtils.cleanUpPath(inputHicFilePaths.get(index));
        for (Integer key : numItersToResults.keySet()) {
            GenomeWideList<SubcompartmentInterval> gwList = numItersToResults.get(key);
            SliceUtils.collapseGWList(gwList);
            File outBedFile = new File(outputDirectory, prefix + key + "_clusters_" + hicFileName + ".subcompartment.bed");
            gwList.simpleExport(outBedFile);
        }
    }

    private void setMseAicBicValues(int z, double[][] iterToWcssAicBic, int numClusters, double sumOfSquares) {
        iterToWcssAicBic[0][z] = numClusters;
        iterToWcssAicBic[1][z] = sumOfSquares;
        // AIC
        iterToWcssAicBic[2][z] = sumOfSquares + 2 * interMatrix.getWidth() * numClusters;
        // BIC .5*k*d*log(n)
        iterToWcssAicBic[3][z] = sumOfSquares + 0.5 * interMatrix.getWidth() * numClusters * Math.log(interMatrix.getLength());
    }
}

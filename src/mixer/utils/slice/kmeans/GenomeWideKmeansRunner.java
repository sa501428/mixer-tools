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

package mixer.utils.slice.kmeans;

import javastraw.reader.basics.ChromosomeHandler;
import mixer.MixerGlobals;
import mixer.utils.drive.DriveMatrix;
import robust.concurrent.kmeans.clustering.Cluster;
import robust.concurrent.kmeans.clustering.KMeansListener;
import robust.concurrent.kmeans.clustering.RobustConcurrentKMeans;
import robust.concurrent.kmeans.clustering.RobustConcurrentKMedians;

import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicBoolean;

public class GenomeWideKmeansRunner {

    private final float[][] matrix;
    private final ChromosomeHandler chromosomeHandler;
    private final DriveMatrix interMatrix;
    private final AtomicBoolean thisRunIsNotDone = new AtomicBoolean(true);
    private KmeansResult result = null;

    private final boolean useCorrMatrix;
    private final boolean useKMedians;

    public GenomeWideKmeansRunner(ChromosomeHandler chromosomeHandler,
                                  DriveMatrix interMatrix,
                                  boolean useCorrMatrix, boolean useKmedians) {
        this.useCorrMatrix = useCorrMatrix;
        this.interMatrix = interMatrix;
        matrix = interMatrix.getData(useCorrMatrix);
        this.chromosomeHandler = chromosomeHandler;
        this.useKMedians = useKmedians;
    }

    public void prepareForNewRun(int numClusters) {
        result = new KmeansResult(numClusters, chromosomeHandler);
        thisRunIsNotDone.set(true);
    }

    public void launchKmeansGWMatrix(long seed, int maxIters) {

        if (matrix.length > 0 && matrix[0].length > 0) {
            if (MixerGlobals.printVerboseComments) {
                System.out.println("Using seed " + seed);
            }

            int numClusters = result.getNumClustersDesired();
            RobustConcurrentKMeans kMeans;
            if (useKMedians) {
                kMeans = new RobustConcurrentKMedians(matrix,
                        numClusters, maxIters, seed);
            } else {
                kMeans = new RobustConcurrentKMeans(matrix,
                        numClusters, maxIters, seed);
            }

            KMeansListener kMeansListener = new KMeansListener() {
                @Override
                public void kmeansMessage(String s) {
                    if (MixerGlobals.printVerboseComments) {
                        System.out.println(s);
                    }
                }

                @Override
                public void kmeansComplete(Cluster[] preSortedClusters) {
                    Cluster[] clusters = ClusterTools.getSortedClusters(preSortedClusters);
                    System.out.print(".");
                    result.processResultAndUpdateScoringMetrics(clusters, interMatrix, useKMedians,
                            useCorrMatrix, 2 * seed);
                    thisRunIsNotDone.set(false);
                }

                @Override
                public void kmeansError(Throwable throwable) {
                    System.err.println("Slice Error - " + throwable.getLocalizedMessage());
                    System.exit(98);
                }
            };
            kMeans.addKMeansListener(kMeansListener);
            kMeans.run();
        }

        waitUntilDone();
    }

    private void waitUntilDone() {
        while (thisRunIsNotDone.get()) {
            System.out.print("*");
            try {
                TimeUnit.SECONDS.sleep(5);
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
        }
    }

    public int getNumColumns() {
        return matrix[0].length;
    }

    public int getNumRows() {
        return matrix.length;
    }

    public KmeansResult getResult() {
        return result;
    }
}

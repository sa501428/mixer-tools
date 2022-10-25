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

package mixer.utils.slice.cleaning.archive;

import javastraw.feature1D.GenomeWide1DList;
import javastraw.reader.Dataset;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.type.NormalizationType;
import javastraw.tools.ExtractingOEDataUtils;
import javastraw.tools.HiCFileTools;
import mixer.utils.similaritymeasures.RobustEuclideanDistance;
import mixer.utils.similaritymeasures.SimilarityMetric;
import mixer.utils.slice.cleaning.BadIndexFinder;
import mixer.utils.slice.cleaning.NearDiagonalTrim;
import mixer.utils.slice.kmeans.ClusterTools;
import mixer.utils.slice.structures.SubcompartmentInterval;

import java.util.*;

public class LeftOverClusterIdentifier {
    private final float threshold = 3f;
    private final ChromosomeHandler chromosomeHandler;
    private final Dataset dataset;
    private final NormalizationType norm;
    private final int resolution;
    private static final SimilarityMetric euclidean = RobustEuclideanDistance.SINGLETON;

    public LeftOverClusterIdentifier(ChromosomeHandler chromosomeHandler, Dataset dataset, NormalizationType norm, int resolution) {
        this.chromosomeHandler = chromosomeHandler;
        this.dataset = dataset;
        this.norm = norm;
        this.resolution = resolution;
    }

    public void identify(Map<Integer, GenomeWide1DList<SubcompartmentInterval>> results,
                         BadIndexFinder badIndexFinder) {

        for (Chromosome chr1 : chromosomeHandler.getAutosomalChromosomesArray()) {
            final MatrixZoomData zd = HiCFileTools.getMatrixZoomData(dataset, chr1, chr1, resolution);
            if (zd == null) continue;

            float[][] allDataForRegion = null;
            try {
                allDataForRegion = HiCFileTools.getOEMatrixForChromosome(dataset, zd, chr1, resolution,
                        norm, threshold, ExtractingOEDataUtils.ThresholdType.TRUE_OE,
                        true, true, 1, 0, true);
                NearDiagonalTrim.nanFill(chr1, allDataForRegion, resolution);

            } catch (Exception e) {
                e.printStackTrace();
                System.exit(99);
            }

            if (allDataForRegion == null) {
                System.err.println("Missing Data " + zd.getKey());
                return;
            }

            Set<Integer> worstIndices = badIndexFinder.getEmptyIndices(chr1);
            Set<Integer> indicesMissing = new HashSet<>();

            for (int k = 0; k < chr1.getLength() / resolution + 1; k++) {
                if (worstIndices.contains(k)) {
                    continue;
                }

                indicesMissing.add(k);
            }

            // do it this way because of additional internal filter
            Integer firstEntryKey = (Integer) results.keySet().toArray()[0];
            for (SubcompartmentInterval handledInterval : results.get(firstEntryKey).getFeatures("" + chr1.getIndex())) {
                int binXStart = handledInterval.getX1() / resolution;
                int binXEnd = handledInterval.getX2() / resolution;

                for (int j = binXStart; j < binXEnd; j++) {
                    indicesMissing.remove(j);
                }
            }

            if (indicesMissing.size() > 0) {
                for (Integer key : results.keySet()) {
                    GenomeWide1DList<SubcompartmentInterval> listForKey = results.get(key);
                    Map<Integer, float[]> cIDToCenter = getClusterCenters(allDataForRegion, listForKey.getFeatures("" + chr1.getIndex()), resolution);

                    List<SubcompartmentInterval> newlyAssignedSubcompartments = getNewlyAssignedCompartments(chr1, cIDToCenter, indicesMissing, allDataForRegion, resolution);

                    listForKey.addAll(newlyAssignedSubcompartments);

                    results.put(key, listForKey);
                }
            }
            System.out.print(".");
        }
        System.out.println(".");
    }

    private Map<Integer, float[]> getClusterCenters(float[][] allDataForRegion, List<SubcompartmentInterval> intervals, int resolution) {

        Map<Integer, float[]> cIDToCenter = new HashMap<>();
        Map<Integer, int[]> cIDToSize = new HashMap<>();

        for (SubcompartmentInterval interval : intervals) {
            int binXStart = interval.getX1() / resolution;
            int binXEnd = interval.getX2() / resolution;
            int cID = interval.getClusterID();
            float[] totalSum = new float[allDataForRegion[0].length];
            int[] totalCounts = new int[allDataForRegion[0].length];

            for (int i = binXStart; i < binXEnd; i++) {
                for (int j = 0; j < allDataForRegion[i].length; j++) {
                    if (!Float.isNaN(allDataForRegion[i][j])) {
                        totalSum[j] += allDataForRegion[i][j];
                        totalCounts[j] += 1;
                    }
                }
            }

            if (cIDToSize.containsKey(cID)) {
                float[] sumVec = cIDToCenter.get(cID);
                for (int j = 0; j < sumVec.length; j++) {
                    totalSum[j] += sumVec[j];
                }

                int[] countVec = cIDToSize.get(cID);
                for (int j = 0; j < countVec.length; j++) {
                    totalCounts[j] += countVec[j];
                }
            }
            cIDToCenter.put(cID, totalSum);
            cIDToSize.put(cID, totalCounts);
        }

        cIDToCenter.replaceAll((k, v) -> ClusterTools.normalize(cIDToCenter.get(k), cIDToSize.get(k)));

        return cIDToCenter;
    }


    private List<SubcompartmentInterval> getNewlyAssignedCompartments(Chromosome chromosome, Map<Integer, float[]> cIDToCenter,
                                                                      Set<Integer> indicesMissing, float[][] allDataForRegion, int resolution) {

        List<SubcompartmentInterval> intervals = new ArrayList<>();

        for (int indx : indicesMissing) {

            int metaID = getClosestClusterID(allDataForRegion[indx], cIDToCenter);

            int x1 = indx * resolution;
            int x2 = x1 + resolution;

            intervals.add(new SubcompartmentInterval(chromosome, x1, x2, metaID));
        }
        return intervals;
    }

    private int getClosestClusterID(float[] vector, Map<Integer, float[]> cIDToCenter) {
        int currID = Integer.MAX_VALUE;
        double overallDistance = Double.MAX_VALUE;
        boolean nothingChanged = true;

        for (Integer key : cIDToCenter.keySet()) {
            double newDistance = euclidean.distance(cIDToCenter.get(key), vector);

            if (newDistance < overallDistance) {
                overallDistance = newDistance;
                currID = key;
                nothingChanged = false;
            }
        }
        if (nothingChanged) {
            System.err.println("Error 787 " + overallDistance + " - " + cIDToCenter.keySet());
        }
        return currID;
    }
}

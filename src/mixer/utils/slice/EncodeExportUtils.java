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

package mixer.utils.slice;

import mixer.MixerGlobals;
import mixer.utils.drive.DriveMatrix;
import mixer.utils.slice.kmeans.KmeansEvaluator;
import mixer.utils.slice.structures.ENCODESubcompartmentInterval;
import mixer.utils.slice.structures.SliceUtils;
import mixer.utils.slice.structures.SubcompartmentInterval;

import java.io.File;
import java.io.FileWriter;
import java.util.*;
import java.util.stream.IntStream;

public class EncodeExportUtils {
    public static void exportSubcompartments(DriveMatrix sliceMatrix, Map<Integer, List<List<Integer>>> kmeansIndicesMap,
                                             KmeansEvaluator evaluator, String prefix, File outputDirectory,
                                             int startingClusterSizeK) {

        List<Integer> keys = new ArrayList<>(kmeansIndicesMap.keySet());
        Collections.sort(keys);

        int[] clusterSizes = convertToIntArray(keys, startingClusterSizeK);
        int bestIndex = getBestIndex(keys, evaluator);
        bestIndex = Math.max(bestIndex, getIndexOf(clusterSizes, 5));

        int[][] ids = new int[sliceMatrix.getData(false).length][keys.size()];
        for (int[] row : ids) {
            Arrays.fill(row, -1);
        }

        for (int index = 0; index < keys.size(); index++) {
            List<List<Integer>> colorList = kmeansIndicesMap.get(keys.get(index));
            for (int i = 0; i < colorList.size(); i++) {
                for (Integer val : colorList.get(i)) {
                    ids[val][index] = i;
                }
            }
        }

        Map<Integer, SubcompartmentInterval> map = sliceMatrix.getRowIndexToIntervalMap();
        List<SubcompartmentInterval> encodeSubcompartmentIntervals = new ArrayList<>();
        for (int i = 0; i < ids.length; i++) {
            if (map.containsKey(i)) {
                SubcompartmentInterval interv = map.get(i);
                if (interv != null) {
                    encodeSubcompartmentIntervals.add(
                            generateENCODESubcompartment(interv, clusterSizes, ids[i], bestIndex));
                }
            }
        }

        encodeSubcompartmentIntervals = SliceUtils.collapseSubcompartmentIntervals(encodeSubcompartmentIntervals);

        File outBedFile = new File(outputDirectory, "slice_subcompartment_clusters.bed");

        try {
            FileWriter filewriter = new FileWriter(outBedFile);
            filewriter.write(ENCODESubcompartmentInterval.getHeader() + "\n");
            filewriter.write("#subcompartments called by mixer-tools/slice version " + MixerGlobals.versionNum + "\n");
            for (SubcompartmentInterval interval : encodeSubcompartmentIntervals) {
                filewriter.write(interval.toString() + "\n");
            }
            filewriter.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    private static ENCODESubcompartmentInterval generateENCODESubcompartment(SubcompartmentInterval interv,
                                                                             int[] clusterSizes, int[] id,
                                                                             int indx) {
        return new ENCODESubcompartmentInterval(interv.getChrIndex(), interv.getChrName(),
                interv.getX1(), interv.getX2(), id[indx], clusterSizes, id);
    }

    private static int[] convertToIntArray(List<Integer> keys, int offset) {
        int[] arr = new int[keys.size()];
        for (int i = 0; i < arr.length; i++) {
            arr[i] = keys.get(i) + offset;
        }
        return arr;
    }

    private static int getBestIndex(List<Integer> keys, KmeansEvaluator evaluator) {
        for (int k = keys.size() - 1; k > -1; k--) {
            Integer key = keys.get(k);
            double corr = evaluator.getCorr(key);
            if (corr > 0 && corr < 0.8) {
                return key;
            }
        }
        return -1;
    }

    private static int getIndexOf(int[] array, int target) {
        return IntStream.range(0, array.length)
                .filter(i -> target == array[i])
                .findFirst() // first occurrence
                .orElse(-1); // No element found
    }
}

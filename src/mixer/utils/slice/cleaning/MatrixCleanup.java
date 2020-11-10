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

package mixer.utils.slice.cleaning;

import javastraw.reader.ExtractingOEDataUtils;
import mixer.utils.common.FloatMatrixTools;
import mixer.utils.common.IntMatrixTools;
import mixer.utils.shuffle.Metrics;
import mixer.utils.slice.structures.SubcompartmentInterval;

import java.io.File;
import java.util.*;

public class MatrixCleanup {
    public static final int BATCHED_NUM_ROWS = 1; // 10
    protected static final float zScoreThreshold = 3f;
    private final static float PERCENT_NAN_ALLOWED = .5f;
    public static boolean USE_COSINE = false;
    public static boolean USE_CORRELATION = false;
    protected final File outputDirectory;
    protected final Random generator;
    protected float[][] data;
    private static final boolean saveIntermediateData = false;

    public MatrixCleanup(float[][] interMatrix, long seed, File outputDirectory) {
        data = ExtractingOEDataUtils.simpleLogWithCleanup(interMatrix, 1);

        System.out.println("matrix size " + data.length + " x " + data[0].length);
        generator = new Random(seed);
        this.outputDirectory = outputDirectory;
    }

    public static Set<Integer> getBadIndices(float[][] matrix) {

        Set<Integer> badIndices = new HashSet<>();

        int[] numNans = getNumberOfNansInRow(matrix);
        int[] numZerosNotNans = getNumberZerosNotNansInRow(matrix);

        float n = matrix[0].length;

        for (int i = 0; i < numNans.length; i++) {
            if ((float) (numNans[i] + numZerosNotNans[i]) / n > PERCENT_NAN_ALLOWED) {
                badIndices.add(i);
            }
        }

        return badIndices;
    }

    public static float[][] filterOutColumnsAndRowsNonSymmetricMatrix(float[][] interMatrix, Map<Integer, SubcompartmentInterval> original) {
        Set<Integer> badIndices = getBadIndices(interMatrix);
        if (badIndices.size() == 0) {
            return interMatrix;
        }

        System.out.println("interMatrix.length " + interMatrix.length + " badIndices.size() " + badIndices.size());

        int counter = 0;
        int[] newIndexToOrigIndex = new int[interMatrix.length - badIndices.size()];
        for (int i = 0; i < interMatrix.length; i++) {
            if (!badIndices.contains(i)) {
                newIndexToOrigIndex[counter++] = i;
            }
        }

        float[][] newMatrix = new float[newIndexToOrigIndex.length][interMatrix[0].length];
        Map<Integer, SubcompartmentInterval> newRowIndexToIntervalMap = new HashMap<>();
        for (int i = 0; i < newMatrix.length; i++) {
            int tempI = newIndexToOrigIndex[i];
            System.arraycopy(interMatrix[tempI], 0, newMatrix[i], 0, newMatrix[0].length);
            newRowIndexToIntervalMap.put(i, (SubcompartmentInterval) original.get(newIndexToOrigIndex[i]).deepClone());
        }

        original.clear();
        original.putAll(newRowIndexToIntervalMap);

        return newMatrix;
    }

    private static int[] getNumberZerosNotNansInRow(float[][] matrix) {
        int[] numZeros = new int[matrix.length];
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                if (!Float.isNaN(matrix[i][j]) && matrix[i][j] < 1e-5) {
                    numZeros[i]++;
                }
            }
        }
        return numZeros;
    }

    private static int[] getNumberOfNansInRow(float[][] matrix) {
        int[] numNans = new int[matrix.length];
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                if (Float.isNaN(matrix[i][j])) {
                    numNans[i]++;
                }
            }
        }
        return numNans;
    }

    public float[][] getSimpleCleaningOfMatrixAppendCorr(Map<Integer, SubcompartmentInterval> rowIndexToIntervalMap) {
        FloatMatrixTools.thresholdNonZerosByZscoreToNanDownColumn(data, zScoreThreshold, BATCHED_NUM_ROWS);

        data = filterOutColumnsAndRowsNonSymmetricMatrix(data, rowIndexToIntervalMap);
        System.out.println("matrix size " + data.length + " x " + data[0].length);

        if (saveIntermediateData) {
            int numCentroids = Math.max(data.length / 50, 200);
            saveMatricesLocally(data, numCentroids, outputDirectory, rowIndexToIntervalMap);
        }

        FloatMatrixTools.inPlaceZscoreDownColsNoNan(data, BATCHED_NUM_ROWS);

        if (USE_COSINE) {
            int numCentroids = Math.max(data.length / 50, 200);
            data = CorrelationTools.getMinimallySufficientNonNanSimilarityMatrix(data,
                    numCentroids, Metrics.Type.COSINE);
            File temp = new File(outputDirectory, "cosine.npy");
            FloatMatrixTools.saveMatrixTextNumpy(temp.getAbsolutePath(), data);
        }
        return data;
    }

    private void saveMatricesLocally(float[][] data, int numCentroids, File outputDirectory, Map<Integer, SubcompartmentInterval> rowIndexToIntervalMap) {
        //0 - cosine, 1 - corr, 2 - mse, 3 - mae, 4 - EMD
        String[] names = {"pre_cosine", "correlation", "pre_mse", "mae", "pre_emd", "pre_kl", "pre_js"};
        String[] suffixes = {"_vs_centroids.npy", "_all.npy"};
        int[] numToUse = {numCentroids, data.length};

        for (int z = 0; z < 2; z++) {
            Metrics.Type k = Metrics.Type.JS;
            //if (k == 1 || k == 3 || k == 4) continue;
            float[][] result = CorrelationTools.getMinimallySufficientNonNanSimilarityMatrix(data, numToUse[z], k);
            File temp = new File(outputDirectory, k.toString() + suffixes[z]);
            FloatMatrixTools.saveMatrixTextNumpy(temp.getAbsolutePath(), result);
            result = null;
            System.gc();
        }

        // indices
        int[][] indices = new int[data.length][2];
        for (int i = 0; i < indices.length; i++) {
            SubcompartmentInterval interval = rowIndexToIntervalMap.get(i);
            indices[i][0] = interval.getChrIndex();
            indices[i][1] = interval.getX1();
        }
        File temp = new File(outputDirectory, "gw_indices.npy");
        IntMatrixTools.saveMatrixTextNumpy(temp.getAbsolutePath(), indices);
    }
}
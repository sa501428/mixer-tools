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

package mixer.utils.common;

import javastraw.tools.MatrixTools;

import javax.imageio.ImageIO;
import java.awt.*;
import java.awt.image.BufferedImage;
import java.io.File;
import java.util.Arrays;
import java.util.List;


/**
 * Helper methods to handle matrix operations
 */
@SuppressWarnings("ALL")
public class FloatMatrixTools {

    public static void thresholdNonZerosByZscoreToNanDownColumn(float[][] matrix, float threshold, int batchSize) {
        float[] colMeans = getColNonZeroMeansNonNan(matrix, batchSize);
        float[] colStdDevs = getColNonZeroStdDevNonNans(matrix, colMeans, batchSize);

        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                float val = matrix[i][j];
                if (!Float.isNaN(val) && val > 1e-10) {
                    int newJ = j / batchSize;
                    float newVal = (val - colMeans[newJ]) / colStdDevs[newJ];
                    if (newVal > threshold) { // || newVal < -threshold || val < 1e-10
                        matrix[i][j] = Float.NaN;
                    }
                }
            }
        }
    }

    public static void inPlaceZscoreDownColsNoNan(float[][] matrix, int batchSize) {
        float[] colMeans = getColNonZeroMeansNonNan(matrix, batchSize);
        float[] colStdDevs = getColNonZeroStdDevNonNans(matrix, colMeans, batchSize);

        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                float val = matrix[i][j];
                if (!Float.isNaN(val)) {
                    int newJ = j / batchSize;
                    matrix[i][j] = (val - colMeans[newJ]) / colStdDevs[newJ];
                }
            }
        }
    }

    public static float[] getColNonZeroStdDevNonNans(float[][] matrix, float[] means, int batchSize) {

        float[] stdDevs = new float[means.length];
        int[] colNonNans = new int[means.length];

        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                float val = matrix[i][j];
                if (!Float.isNaN(val) && val > 1e-10) {
                    int newJ = j / batchSize;
                    float diff = val - means[newJ];
                    stdDevs[newJ] += diff * diff;
                    colNonNans[newJ] += 1;
                }
            }
        }

        for (int k = 0; k < stdDevs.length; k++) {
            stdDevs[k] = (float) Math.sqrt(stdDevs[k] / Math.max(colNonNans[k], 1));
        }

        return stdDevs;
    }

    public static float[] getColNonZeroMeansNonNan(float[][] matrix, int batchSize) {
        float[] colMeans = new float[matrix[0].length / batchSize + 1];
        int[] colNonNans = new int[matrix[0].length / batchSize + 1];
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                float val = matrix[i][j];
                if (!Float.isNaN(val) && val > 1e-10) {
                    int newJ = j / batchSize;
                    colMeans[newJ] += val;
                    colNonNans[newJ] += 1;
                }
            }
        }

        for (int k = 0; k < colMeans.length; k++) {
            colMeans[k] = colMeans[k] / Math.max(colNonNans[k], 1);
        }

        return colMeans;
    }

    public static void fill(float[][] allDataForRegion, float val) {
        for (int i = 0; i < allDataForRegion.length; i++) {
            Arrays.fill(allDataForRegion[i], val);
        }
    }

    public static void saveMatrixTextNumpy(String filename, float[][] matrix) {
        MatrixTools.saveMatrixTextNumpy(filename, matrix);
    }

    public static float[] flattenedRowMajorOrderMatrix(float[][] matrix) {
        int m = matrix.length;
        int n = matrix[0].length;

        int numElements = m * n;
        float[] flattenedMatrix = new float[numElements];

        int index = 0;
        for (int i = 0; i < m; i++) {
            System.arraycopy(matrix[i], 0, flattenedMatrix, index, n);
            index += n;
        }
        return flattenedMatrix;
    }

    public static float[] getRowMajorOrderFlattendedSectionFromMatrix(float[][] matrix, int numCols) {
        int numRows = matrix.length - numCols;

        int numElements = numRows * numCols;
        float[] flattenedMatrix = new float[numElements];

        int index = 0;
        for (int i = numCols; i < numRows; i++) {
            System.arraycopy(matrix[i], 0, flattenedMatrix, index, numCols);
            index += numCols;
        }
        return flattenedMatrix;
    }

    public static float[][] concatenate(float[][] matrix1, float[][] matrix2) {
        float[][] combo = new float[matrix1.length][matrix1[0].length + matrix2[0].length];
        for (int i = 0; i < matrix1.length; i++) {
            System.arraycopy(matrix1[i], 0, combo[i], 0, matrix1[i].length);
            System.arraycopy(matrix2[i], 0, combo[i], matrix1[i].length, matrix2[i].length);
        }
        return combo;
    }

    public static float[][] transpose(float[][] matrix) {
        float[][] result = new float[matrix[0].length][matrix.length];
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                result[j][i] = matrix[i][j];
            }
        }
        return result;
    }

    // column length assumed identical and kept the same
    public static float[][] stitchMultipleMatricesTogetherByRowDim(List<float[][]> data) {
        if (data.size() == 1) return data.get(0);

        int colNums = data.get(0)[0].length;
        int rowNums = 0;
        for (float[][] mtrx : data) {
            rowNums += mtrx.length;
        }

        float[][] aggregate = new float[rowNums][colNums];

        int rowOffSet = 0;
        for (float[][] region : data) {
            copyFromAToBRegion(region, aggregate, rowOffSet, 0);
            rowOffSet += region.length;
        }

        return aggregate;
    }

    public static void copyFromAToBRegion(float[][] source, float[][] destination, int rowOffSet, int colOffSet) {
        for (int i = 0; i < source.length; i++) {
            System.arraycopy(source[i], 0, destination[i + rowOffSet], colOffSet, source[0].length);
        }
    }

    // column length assumed identical and kept the same
    public static float[][] stitchMultipleMatricesTogetherByColDim(List<float[][]> data) {
        if (data.size() == 1) return data.get(0);

        int rowNums = data.get(0).length;
        int colNums = 0;
        for (float[][] mtrx : data) {
            colNums += mtrx[0].length;
        }

        float[][] aggregate = new float[rowNums][colNums];

        int colOffSet = 0;
        for (float[][] region : data) {
            copyFromAToBRegion(region, aggregate, 0, colOffSet);
            colOffSet += region[0].length;
        }

        return aggregate;
    }

    public static float[][] cleanUpMatrix(float[][] matrix) {
        for (int r = 0; r < matrix.length; r++) {
            for (int c = 0; c < matrix[r].length; c++) {
                if (Float.isNaN(matrix[r][c]) || Float.isInfinite(matrix[r][c]) || Math.abs(matrix[r][c]) < 1E-10) {
                    matrix[r][c] = 0;
                }
            }
        }
        return matrix;
    }

    public static void saveMatrixToPNG(File file, float[][] matrix, boolean useLog) {
        double range = getMaxVal(matrix);
        double minVal = 0;
        if (useLog) {
            minVal = Math.log(1 + getMinVal(matrix));
            range = Math.log(1 + range) - minVal;
        }

        BufferedImage image = new BufferedImage(matrix.length, matrix[0].length, BufferedImage.TYPE_INT_RGB);
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[0].length; j++) {
                if (useLog) {
                    image.setRGB(i, j, mixColors((Math.log(1 + matrix[i][j]) - minVal) / range));
                } else {
                    image.setRGB(i, j, mixColors(matrix[i][j] / range));
                }
            }
        }

        try {
            ImageIO.write(image, "png", file);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public static void saveOEMatrixToPNG(File file, double[][] matrix) {
        double max = getMaxAbsLogVal(matrix);
        int zoom = 100;

        BufferedImage image = new BufferedImage(zoom * matrix.length, zoom * matrix[0].length, BufferedImage.TYPE_INT_RGB);
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[0].length; j++) {
                int color = mixOEColors(Math.log(matrix[i][j]), max);
                for (int zi = zoom * i; zi < zoom * (i + 1); zi++) {
                    for (int zj = zoom * j; zj < zoom * (j + 1); zj++) {
                        image.setRGB(zi, zj, color);
                    }
                }
            }
        }

        try {
            ImageIO.write(image, "png", file);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    private static float getMaxAbsLogVal(double[][] matrix) {
        double maxVal = Math.abs(Math.log(matrix[0][0]));
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                double temp = Math.abs(Math.log(matrix[i][j]));
                if (temp > maxVal) {
                    maxVal = temp;
                }
            }
        }
        return (float) maxVal;
    }

    private static float getMaxVal(float[][] matrix) {
        float maxVal = matrix[0][0];
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                if (matrix[i][j] > maxVal) {
                    maxVal = matrix[i][j];
                }
            }
        }
        return maxVal;
    }

    private static float getMinVal(float[][] matrix) {
        float minVal = matrix[0][0];
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                if (matrix[i][j] < minVal) {
                    minVal = matrix[i][j];
                }
            }
        }
        return minVal;
    }

    public static int mixColors(double ratio) {
        int color1 = Color.WHITE.getRGB();
        int color2 = Color.RED.getRGB();

        int mask1 = 0x00ff00ff;
        int mask2 = 0xff00ff00;

        int f2 = (int) (256 * ratio);
        int f1 = 256 - f2;

        return (((((color1 & mask1) * f1) + ((color2 & mask1) * f2)) >> 8) & mask1)
                | (((((color1 & mask2) * f1) + ((color2 & mask2) * f2)) >> 8) & mask2);
    }

    private static int mixOEColors(double val0, double max) {
        int color1 = Color.WHITE.getRGB();
        int color2 = Color.RED.getRGB();
        double val = Math.abs(val0);
        if (val0 < 0) {
            color2 = Color.BLUE.getRGB();
        }

        double ratio = val / max;

        int mask1 = 0x00ff00ff;
        int mask2 = 0xff00ff00;

        int f2 = (int) (256 * ratio);
        int f1 = 256 - f2;

        return (((((color1 & mask1) * f1) + ((color2 & mask1) * f2)) >> 8) & mask1)
                | (((((color1 & mask2) * f1) + ((color2 & mask2) * f2)) >> 8) & mask2);
    }

    public static void log(float[][] matrix, int pseudocount) {
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                matrix[i][j] = (float) Math.log(matrix[i][j] + pseudocount);
            }
        }
    }
}

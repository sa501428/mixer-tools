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

package mixer.utils.slice.gmm;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.RealMatrix;

import java.util.Arrays;
import java.util.List;

public class GMMTools {

    public static float[][] getInitialMeans(List<List<Integer>> groupsOfIndices, int numClusters, float[][] data) {
        float[][] meanVectors = new float[numClusters][data[0].length];
        int[][] counts = new int[numClusters][data[0].length];

        for (int k = 0; k < numClusters; k++) {
            for (int i : groupsOfIndices.get(k)) {
                for (int j = 0; j < data[i].length; j++) {
                    if (!Float.isNaN(data[i][j])) {
                        meanVectors[k][j] += data[i][j];
                        counts[k][j]++;
                    }
                }
            }
        }

        // todo is this a problem???
        for (int k = 0; k < numClusters; k++) {
            for (int j = 0; j < data[0].length; j++) {
                if (counts[k][j] < 1e-10) {
                    System.err.println("Something is going wrong in these calcs");
                }
                meanVectors[k][j] /= Math.max(counts[k][j], 1);
            }
            System.out.println("mu at " + k + " is " + Arrays.toString(meanVectors[k]));
        }

        return meanVectors;
    }

    public static float[][] getNewWeightedAverage(int numClusters, float[][] data, float[][] r) {
        float[][] meanVectors = new float[numClusters][data[0].length];
        float[][] counts = new float[numClusters][data[0].length];
        //float[] n = GMMTools.addUpAllRows(r);

        for (int k = 0; k < numClusters; k++) {
            for (int i = 0; i < data.length; i++) {
                for (int j = 0; j < data[i].length; j++) {
                    if (!Float.isNaN(data[i][j])) {
                        meanVectors[k][j] += r[i][k] * data[i][j];
                        counts[k][j] += r[i][k];
                    }
                }
            }
        }
        // todo is this a problem?
        for (int k = 0; k < numClusters; k++) {
            for (int j = 0; j < data[0].length; j++) {
                if (counts[k][j] < 1e-10) {
                    System.err.println("Something is going wrong in these calcs");
                }
                meanVectors[k][j] /= counts[k][j];
            }
            System.out.println("mu2 at " + k + " is " + Arrays.toString(meanVectors[k]));
            System.out.println("counts2 at " + k + " is " + Arrays.toString(counts[k]));
        }
        System.out.println("---");

        return meanVectors;
    }

    public static float multivariateNormal(float[] x, float[] meanVector, float[][] covarianceMatrix) {
        int[] status = new int[x.length];
        Arrays.fill(status, -1);
        int n = getStatus(x, meanVector, status);
        if (n == 0) {
            System.err.println("x is " + Arrays.toString(x));
            System.err.println("mu is " + Arrays.toString(meanVector));
            System.err.println("cov is " + Arrays.toString(covarianceMatrix));
        }
        RealMatrix cov = getSubsetCovMatrix(covarianceMatrix, status, n);
        LUDecomposition lu = new LUDecomposition(cov);
        RealMatrix diff = validSubtract(x, meanVector, status, n);

        double dval = multivariateNormalCalc(n, lu, diff);
        float val = (float) dval;
        if (Float.isInfinite(val) || Float.isNaN(val)) {
            System.out.println("mnc " + dval + " - " + val);
        }
        return (float) dval;
    }

    public static double multivariateNormalCalc2(int n, LUDecomposition lu, RealMatrix diff) {
        return Math.pow(2 * Math.PI, -n / 2.0) *
                Math.pow(determinant(lu), -0.5) *
                Math.exp(chainMultiply(diff, inverse(lu)) / 2.0);
    }

    public static double multivariateNormalCalc(int n, LUDecomposition lu, RealMatrix diff) {

        double denom = Math.sqrt(Math.pow(2 * Math.PI, n) * determinant(lu));
        if (!Double.isFinite(denom) || denom < 1e-3) {
            System.out.println("denom " + denom);

        }
        double num = Math.exp(chainMultiply(diff, inverse(lu)) / 2.0);
        if (!Double.isFinite(num)) {
            System.out.println("num " + num);

        }
        double result = num / denom;

        if (Double.isNaN(result) || Double.isInfinite(result)) {
            System.out.println("result is " + result);
        }

        return result;
    }

    public static double chainMultiply(RealMatrix diff, RealMatrix inverseCov) {
        RealMatrix matrix = diff.transpose().multiply(inverseCov).multiply(diff);
        if (matrix.getRowDimension() != 1 || matrix.getColumnDimension() != 1) {
            System.err.println("Matrix has wrong dimensions; is " +
                    matrix.getRowDimension() + "x" + matrix.getColumnDimension() +
                    " should be 1x1");
        }
        return matrix.getEntry(0, 0);
    }

    public static RealMatrix inverse(LUDecomposition lu) {
        return lu.getSolver().getInverse();
    }

    public static double determinant(LUDecomposition lu) {
        return lu.getDeterminant();
    }

    public static RealMatrix validSubtract(float[] x, float[] mu, int[] status, int n) {
        double[][] result = new double[n][1];
        for (int k = 0; k < status.length; k++) {
            if (status[k] > -1) {
                result[status[k]][0] = x[k] - mu[k];
            }
        }
        return new Array2DRowRealMatrix(result);
    }

    public static RealMatrix getSubsetCovMatrix(float[][] covMatrix, int[] status, int n) {
        double[][] subset = new double[n][n];
        for (int i = 0; i < covMatrix.length; i++) {
            if (status[i] > -1) {
                for (int j = i; j < covMatrix.length; j++) {
                    if (status[j] > -1) {
                        subset[status[i]][status[j]] = covMatrix[i][j];
                    }
                }
            }
        }
        return new Array2DRowRealMatrix(subset);
    }

    public static int getStatus(float[] x, float[] meanVector, int[] status) {
        int counter = 0;
        for (int i = 0; i < x.length; i++) {
            if (!Float.isNaN(x[i] - meanVector[i])) {
                status[i] = counter;
                counter++;
            }
        }
        if (counter == 0) {
            System.err.println("something is wrong; counter 0");
        }
        return counter;
    }

    public static float[] addUpAllRows(float[][] matrix) {
        float[] result = new float[matrix[0].length];
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                result[j] += matrix[i][j];
            }
        }
        return result;
    }

    public static float[][] getProbabilityOfClusterForRow(int numClusters, float[][] data, float[] pi, float[][] meanVectors, float[][][] covMatrices) {
        float[][] r = new float[data.length][numClusters];
        for (int n = 0; n < data.length; n++) {
            float localSum = 0;
            for (int k = 0; k < numClusters; k++) {
                r[n][k] = pi[k] * multivariateNormal(data[n], meanVectors[k], covMatrices[k]);

                localSum += r[n][k];

                if (Double.isNaN(r[n][k]) || Double.isInfinite(r[n][k])) {
                    printMatrix(covMatrices[k], k);
                    System.err.println("mu is " + Arrays.toString(meanVectors[k]));
                    System.err.println("Pi is " + pi[k]);
                    System.err.println("R is " + r[n][k] + " Local sum " + localSum);
                    System.exit(8);
                }
            }
            for (int k = 0; k < numClusters; k++) {
                r[n][k] /= localSum;
                if (Double.isNaN(r[n][k])) {
                    System.err.println("Local sum " + localSum);
                    System.exit(9);
                }
            }
        }
        return r;
    }

    public static void printMatrix(float[][] matrix, int k) {
        System.out.println("Printing matrix " + k);
        for (float[] c : matrix) {
            System.out.println(Arrays.toString(c));
        }
    }
}

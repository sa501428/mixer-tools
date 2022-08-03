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

package mixer.utils.magic;

import mixer.MixerGlobals;

import java.util.Arrays;

public class FinalScale {

    private final static float tol = .0005f;
    private final static float percentLowRowSumExcluded = 0.0001f;
    private final static float dp = percentLowRowSumExcluded / 2;
    private final static float percentZValuesToIgnore = 0;//0.0025f;
    private final static float dp1 = 0;
    private final static double tolerance = 5e-4;
    private final static int maxIter = 100;
    private final static int totalIterations = 3 * maxIter;
    private final static float minErrorThreshold = .02f;
    private static final float OFFSET = .5f;

    public static float[][] scaleMatrix(SymmLLInterMatrix ic, int[] targetVectorInitial) {

        double low, zHigh, zLow;
        int rLowIndex, zLowIndex, zHighIndex;
        float localPercentLowRowSumExcluded = percentLowRowSumExcluded;
        float localPercentZValuesToIgnore = percentZValuesToIgnore;

        //	find the matrix dimensions
        int k = targetVectorInitial.length;
        float[] current, row, col, rowBackup, dc;
        float[] dr = new float[k];
        int[] bad = new int[k];
        int[] bad1 = new int[k];
        float[] s = new float[k];
        double[] zz = new double[Math.min(k, Integer.MAX_VALUE - 1)];
        double[] r0 = new double[Math.min(k, Integer.MAX_VALUE - 1)];

        int[] zTargetVector = copy(targetVectorInitial);
        float[] calculatedVectorB = new float[k];
        float[] one = new float[k];
        Arrays.fill(one, 1);


        double[] reportErrorForIteration = new double[totalIterations + 3];
        int[] numRoundsForAllIterations = new int[totalIterations + 3];

        int l = 0;
        for (int p = 0; p < k; p++) {
            if (zTargetVector[p] > 0) {
                zz[l++] = zTargetVector[p];
            }
        }
        zz = dealWithSorting(zz, l);

        // unlikely to exceed max int for LowIndex; HighIndex
        // for now we will only sort one vector and hope that suffices
        zLowIndex = (int) Math.max(0, l * localPercentZValuesToIgnore + OFFSET);
        zHighIndex = (int) Math.min(l - 1, l * (1.0 - localPercentZValuesToIgnore) + OFFSET);
        zLow = zz[zLowIndex];
        zHigh = zz[zHighIndex];

        for (int p = 0; p < k; p++) {
            double valZ = zTargetVector[p];
            if (valZ > 0 && (valZ < zLow || valZ > zHigh)) {
                zTargetVector[p] = 0;
            }
        }


        for (int p = 0; p < k; p++) {
            if (zTargetVector[p] == 0) {
                one[p] = 0;
            }
        }

        int[] numNonZero = ic.getNumNonZerosInRow();

        //	find relevant percentiles
        int n0 = 0;
        for (int p = 0; p < k; p++) {
            int valP = numNonZero[p];
            if (valP > 0) {
                r0[n0++] = valP;
            }
        }
        r0 = dealWithSorting(r0, n0);

        rLowIndex = (int) Math.max(0, n0 * localPercentLowRowSumExcluded + OFFSET);
        low = r0[rLowIndex];


        //	find the "bad" rows and exclude them
        for (int p = 0; p < k; p++) {
            if ((numNonZero[p] < low && zTargetVector[p] > 0)) {
                bad[p] = 1;
                zTargetVector[p] = 0;
            }
        }

        row = sparseMultiplyGetRowSums(ic, one);
        rowBackup = copy(row);

        for (int p = 0; p < k; p++) {
            dr[p] = 1 - bad[p];
        }
        dc = copy(dr);
        one = copy(dr);

        // treat separately rows for which z[p] = 0
        for (int p = 0; p < k; p++) {
            if (zTargetVector[p] == 0) {
                one[p] = 0;
            }
        }
        for (int p = 0; p < k; p++) {
            bad1[p] = (int) (1 - one[p]);
        }

        current = copy(dr);
        //	start iterations
        //	row is the current rows sum; dr and dc are the current rows and columns scaling vectors
        double ber = 10.0 * (1.0 + tolerance);
        double err = ber;
        int iter = 0;
        int fail;
        int nErr = 0;
        double[] errors = new double[10000];
        int allIterationsI = 0;

        // if perc or perc1 reached upper bound or the total number of iteration is too high, exit
        while ((ber > tolerance || err > 5.0 * tolerance) && iter < maxIter && allIterationsI < totalIterations
                && localPercentLowRowSumExcluded <= 0.2 && localPercentZValuesToIgnore <= 0.1) {

            iter++;
            allIterationsI++;
            fail = 1;

            col = scaleUpdateSums(ic, bad1, zTargetVector, s, row, dr, dc);
            row = scaleUpdateSums(ic, bad1, zTargetVector, s, col, dc, dr);

            // calculate current scaling vector
            for (int p = 0; p < k; p++) {
                calculatedVectorB[p] = (float) Math.sqrt(dr[p] * dc[p]);
            }

            //	calculate the current error
            ber = 0;
            for (int p = 0; p < k; p++) {
                if (bad1[p] == 1) continue;
                double tempErr = Math.abs(calculatedVectorB[p] - current[p]);
                if (tempErr > ber) {
                    ber = tempErr;
                }
            }

            reportErrorForIteration[allIterationsI - 1] = ber;
            numRoundsForAllIterations[allIterationsI - 1] = iter;

            //	since calculating the error in row sums requires matrix-vector multiplication we are doing this every 10
            //	iterations
            if (iter % 10 == 0) {
                col = sparseMultiplyGetRowSums(ic, calculatedVectorB);
                err = 0;
                for (int p = 0; p < k; p++) {
                    if (bad1[p] == 1) continue;
                    double tempErr = Math.abs((col[p] * calculatedVectorB[p] - zTargetVector[p]));
                    if (err < tempErr) {
                        err = tempErr;
                    }
                }
                errors[nErr++] = err;
            }

            current = copy(calculatedVectorB);

            // check whether convergence rate is satisfactory
            // if less than 5 iterations (so less than 5 errors) and less than 2 row sums errors, there is nothing to check

            if (ber < tolerance && (nErr < 2 || errors[nErr - 1] < 0.5 * errors[nErr - 2])) continue;

            if (iter > 5) {
                for (int q = 1; q <= 5; q++) {
                    if (reportErrorForIteration[allIterationsI - q] * (1.0 + minErrorThreshold) < reportErrorForIteration[allIterationsI - q - 1]) {
                        fail = 0;
                        break;
                    }
                }

                if (nErr >= 2 && errors[nErr - 1] > 0.75 * errors[nErr - 2]) {
                    fail = 1;
                }

                if (iter >= maxIter) {
                    fail = 1;
                }

                if (fail == 1) {
                    localPercentLowRowSumExcluded += dp;
                    localPercentZValuesToIgnore += dp1;
                    nErr = 0;
                    rLowIndex = (int) Math.max(0, n0 * localPercentLowRowSumExcluded + OFFSET);
                    low = r0[rLowIndex];
                    zLowIndex = (int) Math.max(0, l * localPercentZValuesToIgnore + OFFSET);
                    zHighIndex = (int) Math.min(l - 1, l * (1.0 - localPercentZValuesToIgnore) + OFFSET);
                    zLow = zz[zLowIndex];
                    zHigh = zz[zHighIndex];
                    for (int p = 0; p < k; p++) {
                        if (zTargetVector[p] > 0 && (zTargetVector[p] < zLow || zTargetVector[p] > zHigh)) {
                            zTargetVector[p] = 0;
                        }
                    }
                    for (int p = 0; p < k; p++) {
                        if ((numNonZero[p] < low && zTargetVector[p] > 0)) {
                            bad[p] = 1;
                            bad1[p] = 1;
                            one[p] = 0;
                            zTargetVector[p] = 0;
                        }
                    }


                    ber = 10.0 * (1.0 + tol);
                    err = 10.0 * (1.0 + tol);

                    //	if the current error is larger than 5 iteration ago start from scratch,
                    //	otherwise continue from the current position
                    if (reportErrorForIteration[allIterationsI - 1] > reportErrorForIteration[allIterationsI - 6]) {
                        for (int p = 0; p < k; p++) {
                            dr[p] = 1 - bad[p];
                        }
                        dc = copy(dr);
                        one = copy(dr);
                        current = copy(dr);
                        row = copy(rowBackup);
                    } else {
                        for (int p = 0; p < k; p++) {
                            dr[p] *= (1 - bad[p]);
                        }
                        for (int p = 0; p < k; p++) {
                            dc[p] *= (1 - bad[p]);
                        }
                    }
                    iter = 0;
                }
            }
        }

        //	find the final error in row sums
        if (iter % 10 == 0) {
            col = sparseMultiplyGetRowSums(ic, calculatedVectorB);
            err = 0;
            for (int p = 0; p < k; p++) {
                if (bad1[p] == 1) continue;
                double tempErr = Math.abs(col[p] * calculatedVectorB[p] - zTargetVector[p]);
                if (err < tempErr)
                    err = tempErr;
            }
        }

        reportErrorForIteration[allIterationsI + 1] = ber;
        reportErrorForIteration[allIterationsI + 2] = err;

        for (int p = 0; p < k; p++) {
            if (bad[p] == 1) {
                calculatedVectorB[p] = Float.NaN;
            }
        }

        if (MixerGlobals.printVerboseComments) {
            System.out.println(allIterationsI);
            System.out.println(localPercentLowRowSumExcluded);
            System.out.println(localPercentZValuesToIgnore);
            System.out.println(Arrays.toString(reportErrorForIteration));
        }

        return ic.scaleBalance(calculatedVectorB);
    }

    private static float[] scaleUpdateSums(SymmLLInterMatrix ic, int[] bad1, int[] zTargetVector,
                                           float[] s, float[] dim, float[] dDim, float[] dnDim) {
        int k = zTargetVector.length;
        for (int p = 0; p < k; p++) if (bad1[p] == 1) dim[p] = 1.0f;
        for (int p = 0; p < k; p++) s[p] = zTargetVector[p] / dim[p];
        for (int p = 0; p < k; p++) dDim[p] *= s[p];
        float[] newDim = sparseMultiplyGetRowSums(ic, dDim);
        for (int p = 0; p < k; p++) newDim[p] *= dnDim[p];
        return newDim;
    }

    private static float[] copy(float[] original) {
        float[] copy = new float[original.length];
        System.arraycopy(original, 0, copy, 0, copy.length);
        return copy;
    }

    private static int[] copy(int[] original) {
        int[] copy = new int[original.length];
        System.arraycopy(original, 0, copy, 0, copy.length);
        return copy;
    }


    private static double[] dealWithSorting(double[] vector, int length) {
        double[] realVector = new double[length];
        System.arraycopy(vector, 0, realVector, 0, length);
        Arrays.sort(realVector);
        return realVector;
    }

    private static float[] sparseMultiplyGetRowSums(SymmLLInterMatrix ic,
                                                    float[] vector) {
        double[] sumVector = ic.sparseMultiply(vector);
        return convertToFloats(sumVector);
    }

    private static float[] convertToFloats(double[] original) {
        float[] copy = new float[original.length];
        for (int k = 0; k < copy.length; k++) {
            copy[k] = (float) original[k];
        }
        return copy;
    }


}

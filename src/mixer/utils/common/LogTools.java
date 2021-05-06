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

public class LogTools {
    public static void simpleLogWithCleanup(float[][] matrix, float badVal) {
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                float val = matrix[i][j];
                if (!Float.isNaN(val)) {
                    val = (float) Math.log(val + 1);
                    if (Float.isInfinite(val)) {
                        matrix[i][j] = badVal;
                    } else {
                        matrix[i][j] = val;
                    }
                }
            }
        }
    }

    public static float getMaxAbsLogVal(float[][] matrix) {
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

    public static void applySimpleLog(float[][] interMatrix) {
        for (int i = 0; i < interMatrix.length; i++) {
            for (int j = 0; j < interMatrix[i].length; j++) {
                float val = interMatrix[i][j];
                if (Float.isNaN(val)) {
                    interMatrix[i][j] = 0;
                } else {
                    val = (float) Math.log(val + 1);
                    if (Float.isInfinite(val)) {
                        interMatrix[i][j] = 0;
                    } else {
                        interMatrix[i][j] = val;
                    }
                }
            }
        }
    }
}
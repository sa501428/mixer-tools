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

package mixer.utils.shuffle.scoring;

public class DerivativeScoring extends ShuffleScore {

    public DerivativeScoring(float[][] matrix, Integer[] rBounds, Integer[] cBounds) {
        super(matrix, rBounds, cBounds);
    }

    @Override
    public double score() {
        double diff = 0;
        int numElements = 0;
        for (int rI = 0; rI < rBounds.length - 1; rI++) {
            for (int i = rBounds[rI]; i < rBounds[rI + 1] - 1; i++) {

                for (int cI = 0; cI < cBounds.length - 1; cI++) {
                    for (int j = cBounds[cI]; j < cBounds[cI + 1] - 1; j++) {
                        diff += Math.abs(matrix[i][j] - matrix[i + 1][j]);
                        diff += Math.abs(matrix[i][j] - matrix[i][j + 1]);
                        diff += Math.abs(matrix[i][j] - matrix[i + 1][j + 1]);
                        diff += Math.abs(matrix[i + 1][j] - matrix[i][j + 1]);
                        numElements++;
                    }
                }
            }
        }
        return diff / numElements;
    }
}
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

package mixer.utils.shuffle.scoring;

import mixer.utils.matrix.ShuffledIndices;

import java.util.HashMap;
import java.util.Map;

public class VarianceScoring extends ShuffleScore {
    public VarianceScoring(float[][] matrix, ShuffledIndices rBounds, ShuffledIndices cBounds, boolean useSymmetry) {
        super(matrix, rBounds, cBounds, useSymmetry);
    }

    @Override
    protected float score(Integer[] rBounds, Integer[] cBounds, Integer[] rIDs, Integer[] cIDs) {
        double sumOfSquareErr = 0;
        Map<String, Double> sumMap = new HashMap<>();
        Map<String, Long> numRegionMap = new HashMap<>();
        long numElements = populateMeanMap(sumMap, numRegionMap);


        for (int rI = 0; rI < rBounds.length - 1; rI++) {
            for (int cI = 0; cI < cBounds.length - 1; cI++) {
                String key = getKey(rI, cI);
                if (numRegionMap.get(key) > 0) {
                    double mu = sumMap.get(key) / numRegionMap.get(key);
                    for (int i = rBounds[rI]; i < rBounds[rI + 1]; i++) {
                        for (int j = cBounds[cI]; j < cBounds[cI + 1]; j++) {
                            double v = matrix[i][j] - mu;
                            sumOfSquareErr += v * v;
                        }
                    }
                }
            }
        }
        return (float) (sumOfSquareErr / numElements);
    }
}

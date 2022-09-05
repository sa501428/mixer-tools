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

import java.util.Map;

public abstract class ShuffleScore {
    protected final float[][] matrix;
    protected final Integer[] rBounds;
    protected final Integer[] cBounds;

    public ShuffleScore(float[][] matrix, Integer[] rBounds, Integer[] cBounds) {
        this.matrix = matrix;
        this.rBounds = rBounds;
        this.cBounds = cBounds;
    }

    public float score(boolean isBaseline) {
        if (isBaseline) {
            return baselineScore();
        }
        return score(rBounds, cBounds);
    }

    private float baselineScore() {
        return score(new Integer[]{0, matrix.length},
                new Integer[]{0, matrix[0].length});
    }

    protected abstract float score(Integer[] rBounds, Integer[] cBounds);

    protected String getKey(int rI, int cI) {
        if (rI <= cI) {
            return rI + "_" + cI;
        }
        //return cI+"_"+rI;
        return rI + "_" + cI;
    }

    protected long populateMeanMap(Map<String, Double> sumMap, Map<String, Long> numRegionMap) {
        long numElements = 0;
        for (int rI = 0; rI < rBounds.length - 1; rI++) {
            for (int cI = 0; cI < cBounds.length - 1; cI++) {
                double sum = 0;
                long numInRegion = 0;
                for (int i = rBounds[rI]; i < rBounds[rI + 1]; i++) {
                    for (int j = cBounds[cI]; j < cBounds[cI + 1]; j++) {
                        sum += matrix[i][j];
                        numInRegion++;
                    }
                }

                String key = getKey(rI, cI);
                if (sumMap.containsKey(key)) {
                    sumMap.put(key, sumMap.get(key) + sum);
                    numRegionMap.put(key, numRegionMap.get(key) + numInRegion);
                } else {
                    sumMap.put(key, sum);
                    numRegionMap.put(key, numInRegion);
                }
                numElements += numInRegion;
            }
        }
        return numElements;
    }
}

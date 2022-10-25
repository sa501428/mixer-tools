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

import java.util.Map;

public abstract class ShuffleScore {
    protected final float[][] matrix;
    protected final ShuffledIndices rBounds;
    protected final ShuffledIndices cBounds;
    private final boolean useSymmetry;

    public ShuffleScore(float[][] matrix, ShuffledIndices rBounds, ShuffledIndices cBounds, boolean useSymmetry) {
        this.matrix = matrix;
        this.rBounds = rBounds;
        this.cBounds = cBounds;
        this.useSymmetry = useSymmetry;
    }

    public float score(boolean isBaseline) {
        if (isBaseline) {
            return baselineScore();
        }
        return score(rBounds.boundaries, cBounds.boundaries,
                rBounds.ids, cBounds.ids);
    }

    private float baselineScore() {
        return score(new Integer[]{0, matrix.length},
                new Integer[]{0, matrix[0].length},
                new Integer[]{0}, new Integer[]{0});
    }

    protected abstract float score(Integer[] rBounds, Integer[] cBounds,
                                   Integer[] rIDs, Integer[] cIDs);

    protected String getKey(int rI, int cI) {
        if (useSymmetry) {
            int id1 = rBounds.ids[rI];
            int id2 = cBounds.ids[cI];
            if (id1 <= id2) return id1 + "_" + id2;
            return id2 + "_" + id1;
        }
        return rI + "_" + cI;
    }

    protected long populateMeanMap(Map<String, Double> sumMap, Map<String, Long> numRegionMap) {
        long numElements = 0;
        for (int rI = 0; rI < rBounds.boundaries.length - 1; rI++) {
            for (int cI = 0; cI < cBounds.boundaries.length - 1; cI++) {
                double sum = 0;
                long numInRegion = 0;
                for (int i = rBounds.boundaries[rI]; i < rBounds.boundaries[rI + 1]; i++) {
                    for (int j = cBounds.boundaries[cI]; j < cBounds.boundaries[cI + 1]; j++) {
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

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

package mixer.utils.slice.cleaning.utils;

import javastraw.tools.ParallelizationTools;
import mixer.MixerGlobals;
import mixer.utils.slice.matrices.MatrixAndWeight;
import mixer.utils.slice.structures.SubcompartmentInterval;

import java.util.Map;
import java.util.Set;
import java.util.concurrent.atomic.AtomicInteger;

public class RowCleaner extends DimensionCleaner {
    protected final static float PERCENT_NAN_ALLOWED = .7f;
    private final Map<Integer, SubcompartmentInterval> rowIndexToIntervalMap;

    public RowCleaner(float[][] data, Map<Integer, SubcompartmentInterval> rowIndexToIntervalMap, int[] weights) {
        super(data, weights);
        this.rowIndexToIntervalMap = rowIndexToIntervalMap;
    }

    @Override
    protected MatrixAndWeight filterOutBadIndices(Set<Integer> badIndices, float[][] matrix, int[] weights) {
        Map<Integer, SubcompartmentInterval> original = rowIndexToIntervalMap;
        if (MixerGlobals.printVerboseComments) {
            System.out.println("interMatrix.length " + matrix.length + " badIndices.size() " + badIndices.size());
        }

        float[][] newMatrix = MatrixRowCleaner.makeNewMatrixAndUpdateIndices(matrix, original, badIndices);


        return new MatrixAndWeight(newMatrix, weights);
    }

    @Override
    protected int[] getNumberOfNansInDimension(float[][] matrix) {
        // invalid rows
        int[] totalNumInvalids = new int[matrix.length];

        AtomicInteger index = new AtomicInteger(0);
        ParallelizationTools.launchParallelizedCode(() -> {
            int i = index.getAndIncrement();
            while (i < data.length) {
                int numInvalids = 0;
                for (int j = 0; j < matrix[i].length; j++) {
                    if (Float.isNaN(matrix[i][j])) {
                        numInvalids++;
                    }
                }
                totalNumInvalids[i] = numInvalids;
                i = index.getAndIncrement();
            }
        });

        return totalNumInvalids;
    }

    @Override
    protected int getLimit() {
        return (int) Math.ceil(PERCENT_NAN_ALLOWED * data[0].length);
    }

    @Override
    protected boolean useOnlyCorr() {
        return false;
    }

    @Override
    protected float[][] getAppropriatelyFlippedMatrix() {
        return data;
    }
}

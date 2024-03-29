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

package mixer.utils.slice.matrices;

import javastraw.reader.basics.Chromosome;

import java.util.Map;

/**
 * Container class for dimensions of a matrix and tracking indices
 */

public class Dimension {
    public int length = 0;
    public final int[] offset;
    public int[] interval;

    // simple binning
    public Dimension(Chromosome[] chromosomes, int resolution) {
        offset = new int[chromosomes.length];
        interval = new int[chromosomes.length];
        for (int i = 0; i < chromosomes.length; i++) {
            length += (int) (chromosomes[i].getLength() / resolution + 1);
            if (i < chromosomes.length - 1) {
                offset[i + 1] = length;
            }
            interval[i] = length;
        }
    }

    // compressed binning
    public Dimension(Chromosome[] chromosomes, Map<Integer, Integer> indexToFilteredLength) {
        offset = new int[chromosomes.length];
        interval = new int[chromosomes.length];
        for (int i = 0; i < chromosomes.length; i++) {
            int val = indexToFilteredLength.get(chromosomes[i].getIndex());
            length += val;
            if (i < chromosomes.length - 1) {
                offset[i + 1] = length;
            }
            interval[i] = val;
        }
    }
}
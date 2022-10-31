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

package mixer.utils.shuffle;

import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ChromosomeArrayPair;
import javastraw.reader.basics.ChromosomeHandler;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class Partition {

    // batch 2s / split 1
    private final static int[] chroms1A = new int[]{1, 4, 5, 8, 9, 12, 13, 16, 17, 20, 21};
    private final static int[] chroms1B = new int[]{2, 3, 6, 7, 10, 11, 14, 15, 18, 19, 22};

    // split 2
    private final static int[] chroms2A = new int[]{1, 3, 4, 5, 8, 11, 13, 14, 15, 16, 17, 18};
    private final static int[] chroms2B = new int[]{2, 6, 7, 9, 10, 12, 19, 20, 21, 22};

    // split 3
    private final static int[] chroms3A = new int[]{1, 3, 6, 8, 10, 11, 12, 17, 21};
    private final static int[] chroms3B = new int[]{2, 4, 5, 7, 9, 13, 14, 15, 16, 18, 19, 20, 22};

    // split 4
    private final static int[] chroms4A = new int[]{1, 4, 6, 7, 9, 11, 12, 15, 16, 17, 18, 22};
    private final static int[] chroms4B = new int[]{2, 3, 5, 8, 10, 13, 14, 19, 20, 21};

    // split 5
    private final static int[] chroms5A = new int[]{1, 2, 3, 5, 6, 9, 15, 16, 21, 22};
    private final static int[] chroms5B = new int[]{4, 7, 8, 10, 11, 12, 13, 14, 17, 18, 19, 20};

    public static ChromosomeArrayPair getChromosomePartition(ChromosomeHandler chromosomeHandler, Type mapType) {
        if (mapType == Type.SKIP_BY_TWOS) {
            return chromosomeHandler.splitAutosomesAndSkipByTwos();
        } else if (mapType == Type.FIRST_HALF_VS_SECOND_HALF) {
            return chromosomeHandler.splitAutosomesIntoHalves();
        } else if (mapType == Type.ODDS_VS_EVENS) {
            return new ChromosomeArrayPair(chromosomeHandler.extractOddOrEvenAutosomes(true),
                    chromosomeHandler.extractOddOrEvenAutosomes(false));
        } else if (mapType == Type.SPLIT1) {
            return createPartition(chromosomeHandler, chroms1A);
        } else if (mapType == Type.SPLIT2) {
            return createPartition(chromosomeHandler, chroms2A);
        } else if (mapType == Type.SPLIT3) {
            return createPartition(chromosomeHandler, chroms3A);
        } else if (mapType == Type.SPLIT4) {
            return createPartition(chromosomeHandler, chroms4A);
        } else if (mapType == Type.SPLIT5) {
            return createPartition(chromosomeHandler, chroms5A);
        }
        return null;
    }

    private static ChromosomeArrayPair createPartition(ChromosomeHandler handler, int[] chroms) {
        Set<Integer> partA = convertToSet(chroms);
        List<Chromosome> cA = new ArrayList<>(20);
        List<Chromosome> cB = new ArrayList<>(20);
        for (Chromosome chromosome : handler.getAutosomalChromosomesArray()) {
            if (partA.contains(chromosome.getIndex())) {
                cA.add(chromosome);
            } else {
                cB.add(chromosome);
            }
        }

        return new ChromosomeArrayPair(cA.toArray(new Chromosome[0]),
                cB.toArray(new Chromosome[0]));
    }

    private static Set<Integer> convertToSet(int[] chroms) {
        Set<Integer> chromSet = new HashSet<>();
        for (int val : chroms) {
            chromSet.add(val);
        }
        return chromSet;
    }

    public enum Type {
        ODDS_VS_EVENS, FIRST_HALF_VS_SECOND_HALF, SKIP_BY_TWOS,
        SPLIT1, SPLIT2, SPLIT3, SPLIT4, SPLIT5
    }
}

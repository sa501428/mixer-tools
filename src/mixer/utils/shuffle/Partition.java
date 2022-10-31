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

import javastraw.reader.basics.ChromosomeArrayPair;
import javastraw.reader.basics.ChromosomeHandler;

public class Partition {
    public static ChromosomeArrayPair getChromosomePartition(ChromosomeHandler chromosomeHandler, Type mapType) {
        if (mapType == Type.SKIP_BY_TWOS) {
            return chromosomeHandler.splitAutosomesAndSkipByTwos();
        } else if (mapType == Type.FIRST_HALF_VS_SECOND_HALF) {
            return chromosomeHandler.splitAutosomesIntoHalves();
        } else if (mapType == Type.ODDS_VS_EVENS) {
            return new ChromosomeArrayPair(chromosomeHandler.extractOddOrEvenAutosomes(true),
                    chromosomeHandler.extractOddOrEvenAutosomes(false));
        }
        return null;
    }

    public enum Type {ODDS_VS_EVENS, FIRST_HALF_VS_SECOND_HALF, SKIP_BY_TWOS}
}

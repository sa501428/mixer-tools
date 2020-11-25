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

package mixer.utils.shuffle;

import mixer.utils.common.FloatMatrixTools;

import java.io.File;

public class AggregateMatrix {

    private final String stem;
    private float[][] aggregate = null;

    public AggregateMatrix(boolean isBaseline) {
        if (isBaseline) {
            stem = "baseline";
        } else {
            stem = "shuffled";
        }
    }

    private static void addBToA(float[][] a, float[][] b) {
        if (a.length == b.length && a[0].length == b[0].length) {
            for (int i = 0; i < a.length; i++) {
                for (int j = 0; j < a[i].length; j++) {
                    a[i][j] += b[i][j];
                }
            }
        } else {
            System.err.println("dimensions incorrect " + a.length + "==" + b.length
                    + "; " + a[0].length + "==" + b[0].length);
        }
    }

    private static void divideBy(float[][] matrix, float scalar) {
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                matrix[i][j] /= scalar;
            }
        }
    }

    public synchronized void addBToA(float[][] matrix) {
        if (aggregate == null) {
            aggregate = new float[matrix.length][matrix[0].length];
        }
        addBToA(aggregate, matrix);
    }

    public void scaleForNumberOfRounds(int numRounds) {
        divideBy(aggregate, numRounds);
    }

    public void saveToPNG(File outfolder, InterOnlyMatrix.InterMapType mapType) {
        File mapFile = new File(outfolder, mapType.toString() + "_" + stem + ".png");
        File mapLogFile = new File(outfolder, mapType.toString() + "_log_" + stem + ".png");
        FloatMatrixTools.saveMatrixToPNG(mapFile, aggregate, false);
        FloatMatrixTools.saveMatrixToPNG(mapLogFile, aggregate, true);

        //FloatMatrixTools.saveMatrixTextNumpy(new File(outfolder, mapTypes[y].toString() + "_baseline.npy").getAbsolutePath(), baselineM);
        //FloatMatrixTools.saveMatrixTextNumpy(new File(outfolder, mapTypes[y].toString() + "_shuffle.npy").getAbsolutePath(), shuffleM);
        //FloatMatrixTools.saveMatrixTextNumpy(new File(outfolder, mapTypes[y].toString() + "_matrix.npy").getAbsolutePath(), interMatrix.getMatrix());
        //EntropyCalculations ec = new EntropyCalculations(shuffleFile, baselineFile, shuffleLogFile, baselineLogFile, shuffleM);
    }
}

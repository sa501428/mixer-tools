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

package mixer.utils.custom;

import javastraw.feature1D.GenomeWide1DList;
import javastraw.reader.Dataset;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.type.NormalizationType;
import javastraw.tools.ExtractingOEDataUtils;
import javastraw.tools.HiCFileTools;
import mixer.utils.slice.cleaning.NearDiagonalTrim;
import mixer.utils.slice.structures.SubcompartmentInterval;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class ExtractionB4 {

    private static final int RESOLUTION = 50000;

    public static void extract() {
        Dataset ds = HiCFileTools.extractDatasetForCLT(
                "/Users/mshamim/Desktop/hicfiles/gm12878_rh14_30.hic",
                false, false, false);
        Chromosome chrom19 = ds.getChromosomeHandler().getChromosomeFromName("19");

        final MatrixZoomData zd = HiCFileTools.getMatrixZoomData(ds, chrom19, chrom19, RESOLUTION);
        NormalizationType type = ds.getNormalizationHandler().getNormTypeFromString("KR");

        try {
            float[][] matrix = HiCFileTools.getOEMatrixForChromosome(ds, zd, chrom19, RESOLUTION,
                    type, 10f,
                    ExtractingOEDataUtils.ThresholdType.TRUE_OE_LOG,
                    //ExtractingOEDataUtils.ThresholdType.TRUE_OE,
                    true, true, 1, 0, true);
            removeEmptyValues(matrix);
            NearDiagonalTrim.nanFill(chrom19, matrix, RESOLUTION);
            float[] sums = getAbsRowSums(matrix);
            fillEmptyRows(matrix, sums);

            for (int numCentroids = 5; numCentroids < 12; numCentroids++) {
                QuickClustering clustering = new QuickClustering(matrix, numCentroids, 13768L, 1000);
                int[] result = clustering.cluster();

                List<SubcompartmentInterval> intervalList = new ArrayList<>();
                for (int i = 0; i < result.length; i++) {
                    SubcompartmentInterval interval = new SubcompartmentInterval(chrom19,
                            i * RESOLUTION, (i + 1) * RESOLUTION, result[i]);
                    intervalList.add(interval);
                }
                GenomeWide1DList<SubcompartmentInterval> finalCompartments = new GenomeWide1DList<>(ds.getChromosomeHandler());
                finalCompartments.addAll(intervalList);

                File outBedFile = new File("/Users/mshamim/Desktop/B4_Gold",
                        "B4_LOG_" + numCentroids + "_clusters.bed"
                        //"B4_TOE_"+numCentroids + "_clusters.bed"
                );
                finalCompartments.simpleExport(outBedFile);
            }

        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    private static void removeEmptyValues(float[][] matrix) {
        for (int i = 0; i < matrix.length; i++) {
            for (int j = i; j < matrix[i].length; j++) {
                if (Float.isNaN(matrix[i][j])) continue;
                if (Math.abs(matrix[i][j]) < 1e-10) {
                    // todo
                }
            }
        }
    }

    private static void fillEmptyRows(float[][] matrix, float[] sums) {
        for (int i = 0; i < sums.length; i++) {
            if (sums[i] < 1e-5) {
                Arrays.fill(matrix[i], -10000);
            }
        }
    }

    private static float[] getAbsRowSums(float[][] matrix) {
        float[] sums = new float[matrix[0].length];
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                float val = matrix[i][j];
                if (!Float.isNaN(val)) {
                    sums[i] += Math.abs(val);
                }
            }
        }
        return sums;
    }
}

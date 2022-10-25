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

package mixer.utils.matrix;

import javastraw.reader.Dataset;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.type.NormalizationType;
import javastraw.tools.ExtractingOEDataUtils;
import javastraw.tools.HiCFileTools;
import mixer.utils.common.LogTools;
import mixer.utils.similaritymeasures.SimilarityMetric;
import mixer.utils.slice.cleaning.SimilarityMatrixTools;

public class IntraOnlyMatrix extends HiCMatrix {

    public IntraOnlyMatrix(Dataset ds, NormalizationType norm, int resolution, Chromosome chromosome,
                           INTRA_TYPE intra_type, SimilarityMetric metric, int compressionFactor) {
        super(ds, norm, resolution, new Chromosome[]{chromosome}, new Chromosome[]{chromosome},
                metric, true, intra_type, true, compressionFactor);
    }

    protected void fillInChromosomeRegion(Dataset ds, float[][] matrix, MatrixZoomData zd, Chromosome chr1, int offsetIndex1,
                                          Chromosome chr2, int offsetIndex2, boolean needToFlip) {

        int lengthChr1 = (int) Math.ceil((float) chr1.getLength() / resolution);
        int lengthChr2 = (int) Math.ceil((float) chr2.getLength() / resolution);

        float[][] allDataForRegion = null;
        try {
            if (intraType == INTRA_TYPE.DEFAULT) {
                allDataForRegion = HiCFileTools.extractLocalBoundedRegionFloatMatrix(zd, 0, lengthChr1,
                        0, lengthChr2, lengthChr1, lengthChr2, norm, true);
            } else if (intraType == INTRA_TYPE.LOG_BASE_EXPECTED) {
                allDataForRegion = HiCFileTools.getOEMatrixForChromosome(ds, zd, chr1, resolution,
                        norm, 10, ExtractingOEDataUtils.ThresholdType.LOG_BASE_EXP_OF_OBS,
                        true, true, 1, 0, true);
            } else if (intraType == INTRA_TYPE.TRUE_OE) {
                allDataForRegion = HiCFileTools.getOEMatrixForChromosome(ds, zd, chr1, resolution,
                        norm, 10, ExtractingOEDataUtils.ThresholdType.TRUE_OE,
                        true, true, 1, 0, true);
            } else if (intraType == INTRA_TYPE.JUST_LOG) {
                allDataForRegion = HiCFileTools.extractLocalBoundedRegionFloatMatrix(zd, 0, lengthChr1,
                        0, lengthChr2, lengthChr1, lengthChr2, norm, true);
                LogTools.simpleLogWithCleanup(allDataForRegion, 0);
            } else {
                System.err.println("Invalid Matrix type " + intraType);
                System.exit(9);
            }
            if (metric != null) {
                allDataForRegion = SimilarityMatrixTools.getNonNanSimilarityMatrix(allDataForRegion,
                        metric, 1, 283746L);
            }
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(99);
        }

        for (int i = 0; i < allDataForRegion.length; i++) {
            System.arraycopy(allDataForRegion[i], 0,
                    matrix[offsetIndex1 + i], offsetIndex2, allDataForRegion[i].length);
        }
    }
}

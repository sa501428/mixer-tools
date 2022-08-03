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

package mixer.utils.common;

import javastraw.reader.Dataset;
import javastraw.reader.norm.NormalizationVector;
import javastraw.reader.type.HiCZoom;
import javastraw.reader.type.NormalizationType;
import javastraw.tools.HiCFileTools;
import javastraw.tools.MatrixTools;

import java.io.File;

public class Unused {

    public static void run() {
        Dataset ds = HiCFileTools.extractDatasetForCLT(
                "/Users/mshamim/Desktop/hicfiles/SCALE/hap1_SCALE_30.hic",
                false, false, false);
        File outfolder = new File("/Users/mshamim/Desktop/out_norms");
        NormalizationType gw_norm = ds.getNormalizationTypesMap().get("GW_SCALE");
        NormalizationType inter_norm = ds.getNormalizationTypesMap().get("INTER_SCALE");
        for (int i = 1; i < 23; i++) {
            NormalizationVector nv_gw = ds.getNormalizationVector(i, new HiCZoom(HiCZoom.HiCUnit.BP, 25000), gw_norm);
            NormalizationVector nv_inter = ds.getNormalizationVector(i, new HiCZoom(HiCZoom.HiCUnit.BP, 25000), inter_norm);

            float[][] result = new float[2][(int) nv_gw.getData().getLength()];
            for (int k = 0; k < nv_gw.getData().getLength(); k++) {
                result[0][k] = (float) nv_gw.getData().get(k);
                result[1][k] = (float) nv_inter.getData().get(k);
            }

            File outfile = new File(outfolder, i + ".npy");
            MatrixTools.saveMatrixTextNumpy(outfile.getAbsolutePath(), result);
        }
    }
}

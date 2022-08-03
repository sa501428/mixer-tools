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

package mixer.utils.estimate;

import javastraw.reader.Dataset;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.block.Block;
import javastraw.reader.block.ContactRecord;
import javastraw.reader.mzd.Matrix;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.type.HiCZoom;
import javastraw.reader.type.NormalizationType;

import java.util.List;

public class EstimateTarget {
    private static long getGenomeSize(Chromosome[] chromosomes) {
        long genomeLength = 0L;
        for (Chromosome chromosome : chromosomes) {
            genomeLength += chromosome.getLength();
        }
        return genomeLength;
    }

    public void estimate(Dataset ds) {
        int resolution = 500000;
        HiCZoom zoom = ds.getZoomForBPResolution(resolution);
        NormalizationType norm = ds.getNormalizationHandler().getNormTypeFromString("NONE");
        Chromosome[] chromosomes = ds.getChromosomeHandler().getAutosomalChromosomesArray();

        long totalNumberOfInterContacts = 0L;

        for (int i = 0; i < chromosomes.length; i++) {
            for (int j = i + 1; j < chromosomes.length; j++) {

                Matrix matrix = ds.getMatrix(chromosomes[i], chromosomes[j]);
                if (matrix == null) continue;

                MatrixZoomData zd = matrix.getZoomData(zoom);
                if (zd == null) continue;

                long rowEnd = chromosomes[i].getLength() / resolution;
                long columnEnd = chromosomes[j].getLength() / resolution;
                List<Block> blocks = zd.getNormalizedBlocksOverlapping(0, 0, rowEnd, columnEnd, norm,
                        false);

                for (Block block : blocks) {
                    for (ContactRecord cr : block.getContactRecords()) {
                        float counts = cr.getCounts();
                        if (Float.isNaN(counts) || Float.isInfinite(counts)) continue;
                        totalNumberOfInterContacts += counts;
                    }
                }
            }
        }

        totalNumberOfInterContacts *= 2;
        System.out.println("Total number of Inter-chromosomal contacts is " + totalNumberOfInterContacts);
        double estimate = (double) (getGenomeSize(chromosomes) / (totalNumberOfInterContacts / (200 * 5)));
        System.out.println("Estimate of SLICE resolution " + estimate);
    }
}

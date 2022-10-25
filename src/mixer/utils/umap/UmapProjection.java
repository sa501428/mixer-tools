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

package mixer.utils.umap;

import mixer.utils.drive.DriveMatrix;
import mixer.utils.slice.gmm.SimpleScatterPlot;
import mixer.utils.slice.structures.SubcompartmentInterval;
import tagbio.umap.Umap;

import java.io.File;
import java.util.Arrays;
import java.util.List;
import java.util.Map;

public class UmapProjection {

    private final float[][] umapProjection;
    private final Map<Integer, SubcompartmentInterval> rowIndexToIntervalMap;

    public UmapProjection(DriveMatrix matrix, boolean useCorrelation) {
        umapProjection = get2DProjection(matrix.getData(useCorrelation));
        rowIndexToIntervalMap = matrix.getRowIndexToIntervalMap();
    }

    public static float[][] get2DProjection(float[][] points) {
        if (points[0].length > 2) {
            final Umap umap = new Umap();
            umap.setNumberComponents(2);
            umap.setNumberNearestNeighbours(50); // 50 // 15 -> 50 for more global picture
            umap.setThreads(10);
            umap.setMinDist(0.5f); // 0.2 ->0.8 -> 0.5  //0.1f -> 0.2f for more general features
            umap.setVerbose(false);
            umap.setSeed(0L);
            return umap.fitTransform(points);
        }
        return points;
    }

    public void plotProjection(File outputDirectory, List<List<Integer>> colorList, String stem) {
        int[] colors = new int[umapProjection.length];
        Arrays.fill(colors, -1);
        for (int index = 0; index < colorList.size(); index++) {
            for (Integer val : colorList.get(index)) {
                colors[val] = index;
            }
        }

        plotProjection(outputDirectory, colors, stem);
    }


    public void plotProjection(File outputDirectory, int[] colors, String stem) {
        SimpleScatterPlot plotter = new SimpleScatterPlot(umapProjection);
        File outfile = new File(outputDirectory, "umap_" + stem);
        plotter.plot(colors, outfile.getAbsolutePath());
    }

    public void runUmapAndColorByChromosome(File outputDirectory) {
        int[] indices = new int[umapProjection.length];
        for (int i = 0; i < indices.length; i++) {
            SubcompartmentInterval interval = rowIndexToIntervalMap.get(i);
            indices[i] = interval.getChrIndex();
        }

        plotProjection(outputDirectory, indices, "chrom");
    }
}

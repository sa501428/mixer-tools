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

package mixer.algos;

import javastraw.feature2D.Feature2D;
import javastraw.feature2D.Feature2DList;
import javastraw.feature2D.Feature2DParser;
import javastraw.reader.Dataset;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.mzd.Matrix;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.type.HiCZoom;
import javastraw.reader.type.NormalizationType;
import javastraw.tools.HiCFileTools;
import javastraw.tools.MatrixTools;
import mixer.clt.CommandLineParserForMixer;
import mixer.clt.MixerCLT;
import mixer.utils.network.MapQFilter;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Created by muhammadsaadshamim on 9/14/15.
 */
public class FineTune extends MixerCLT {

    private String inputBedpeFile, outputPath;
    private int resolution = 5000;
    private Dataset ds;

    public FineTune() {
        super("finetune [-r resolution] [-k norm] <hicfile> <bedpe> [output_path]");
    }

    @Override
    protected void readMixerArguments(String[] args, CommandLineParserForMixer mixerParser) {
        if (args.length != 4) {
            printUsageAndExit(7);
        }

        inputBedpeFile = args[2];
        outputPath = args[3];
        ds = HiCFileTools.extractDatasetForCLT(args[1], false, false, false);

        NormalizationType preferredNorm = mixerParser.getNormalizationTypeOption(ds.getNormalizationHandler());
        if (preferredNorm != null)
            norm = preferredNorm;

        try {
            int specifiedResolution = mixerParser.getMultipleResolutionOptions().get(0);
            resolution = specifiedResolution;
        } catch (Exception e) {
        }
    }


    @Override
    public void run() {

        ChromosomeHandler handler = ds.getChromosomeHandler();
        if (givenChromosomes != null)
            handler = HiCFileTools.stringToChromosomes(givenChromosomes, handler);

        Feature2DList listA = Feature2DParser.loadFeatures(inputBedpeFile, handler, false, null, false);

        MapQFilter.removeLowMapQFeatures(listA, 5000, ds, handler, norm);

        chooseMaxValPixel(listA, handler);

    }

    private void chooseMaxValPixel(Feature2DList list, ChromosomeHandler handler) {
        final Map<String, Chromosome> chrNameToIndex = new HashMap<>();
        for (Chromosome chr : handler.getChromosomeArray()) {
            chrNameToIndex.put(Feature2DList.getKey(chr, chr), chr);
        }
        list.filterLists((String chr, List<Feature2D> feature2DList) -> {
            try {
                return locationsOfMaxVals(resolution, chrNameToIndex.get(chr), ds, feature2DList, norm);
            } catch (Exception e) {
                System.err.println("Unable to remove low mapQ entries for " + chr);
            }
            return new ArrayList<>();
        });

    }

    private List<Feature2D> locationsOfMaxVals(int resolution, Chromosome chromosome,
                                               Dataset ds, List<Feature2D> feature2DList,
                                               NormalizationType norm) {
        List<Feature2D> newList = new ArrayList<>();
        Matrix matrix = ds.getMatrix(chromosome, chromosome);
        MatrixZoomData zd = matrix.getZoomData(new HiCZoom(HiCZoom.HiCUnit.BP, resolution));

        for (Feature2D loop : feature2DList) {
            long start1 = loop.getStart1() / resolution - 1;
            long start2 = loop.getStart2() / resolution - 1;
            long end1 = loop.getEnd1() / resolution + 1;
            long end2 = loop.getEnd2() / resolution + 1;

            try {
                double[][] data = HiCFileTools.extractLocalBoundedRegion(zd, start1, end1,
                        start2, end2, (int) (end1 - start1), (int) (end2 - start2),
                        norm, false).getData();
                System.out.println(chromosome.getName());
                System.out.println(start1 * resolution + " " + end1 * resolution);
                System.out.println(start2 * resolution + " " + end2 * resolution);
                MatrixTools.saveMatrixTextNumpy(new File(outputPath, "out.npy").getAbsolutePath(),
                        data);
                System.out.println();

                System.exit(6);
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
        return newList;
    }

}

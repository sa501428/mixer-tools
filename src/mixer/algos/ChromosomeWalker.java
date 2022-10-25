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

import javastraw.reader.Dataset;
import javastraw.reader.type.NormalizationType;
import javastraw.tools.HiCFileTools;
import mixer.clt.CommandLineParserForMixer;
import mixer.clt.MixerCLT;
import mixer.utils.walker.LocalGenomeRegion;

import java.io.File;
import java.io.FileOutputStream;
import java.io.ObjectOutputStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;

public class ChromosomeWalker extends MixerCLT {

    private static int expectedMinSize = 3;
    private final int numPixelOverlapWhileSliding = 5;
    private Dataset ds;
    private File outputDirectory;
    private int resolution = 1000;
    private int neuralNetSize = 500;

    //whole assembly linker
    public ChromosomeWalker() {
        super("walk [-r res1,res2] [-k normalization] [-w minimum size]" +
                "  <hicFile> <assembly_file> <output_directory>");
    }

    public static void writeStrictIntsToFile(File path, List<Integer> positions) {
        try (ObjectOutputStream write = new ObjectOutputStream(new FileOutputStream(path))) {
            for (Integer pos : positions) {
                write.writeObject(pos);
            }
        } catch (Exception eo) {
            eo.printStackTrace();
        }
    }

    public static void writeStrictMapToFile(File path, Map<Integer, LocalGenomeRegion> indexToRegion) {
        List<Integer> keys = new ArrayList<>(indexToRegion.keySet());
        Collections.sort(keys);

        try (ObjectOutputStream write = new ObjectOutputStream(new FileOutputStream(path))) {
            for (Integer key : keys) {
                LocalGenomeRegion region = indexToRegion.get(key);

                write.writeObject(region.toString());
            }
        } catch (Exception eo) {
            eo.printStackTrace();
        }
    }

    @Override
    protected void readMixerArguments(String[] args, CommandLineParserForMixer mixerParser) {
        if (args.length != 4) {
            // 3 - standard, 5 - when list/control provided
            printUsageAndExit(8);  // this will exit
        }

        ds = HiCFileTools.extractDatasetForCLT(args[1], false, false, true);
        outputDirectory = HiCFileTools.createValidDirectory(args[3]);


        NormalizationType preferredNorm = mixerParser.getNormalizationTypeOption(ds.getNormalizationHandler());
        if (preferredNorm != null)
            norm = preferredNorm;

        List<Integer> potentialResolution = mixerParser.getMultipleResolutionOptions();
        if (potentialResolution != null) {
            resolution = potentialResolution.get(0);
        }

        int specifiedCliqueSize = mixerParser.getWindowSizeOption();
        if (specifiedCliqueSize > 1) {
            expectedMinSize = specifiedCliqueSize;
        }

        updateNumberOfCPUThreads(mixerParser);

        int specifiedMatrixSize = mixerParser.getMatrixSizeOption();
        if (specifiedMatrixSize > 10) {
            neuralNetSize = specifiedMatrixSize;
        }
    }


    @Override
    public void run() {


    }
}

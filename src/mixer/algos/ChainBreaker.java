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

import java.io.File;
import java.util.List;

public class ChainBreaker extends MixerCLT {

    public static int expectedCliqueSize = 5;
    private Dataset ds;
    private File outputDirectory;
    private int initialResolution = 25000;

    protected ChainBreaker(String usage) {
        super("chainbreaker  <hicFile> <assembly_file> <output_directory>");
    }

    @Override
    protected void readMixerArguments(String[] args, CommandLineParserForMixer mixerParser) {
        if (args.length != 4) {
            printUsageAndExit(2);  // this will exit
        }

        ds = HiCFileTools.extractDatasetForCLT(args[1], false, false, true);
        outputDirectory = HiCFileTools.createValidDirectory(args[3]);

        NormalizationType preferredNorm = mixerParser.getNormalizationTypeOption(ds.getNormalizationHandler());
        if (preferredNorm != null)
            norm = preferredNorm;

        List<Integer> potentialResolution = mixerParser.getMultipleResolutionOptions();
        if (potentialResolution != null) {
            initialResolution = potentialResolution.get(0);
        }

        updateNumberOfCPUThreads(mixerParser);
    }

    @Override
    public void run() {


    }
}

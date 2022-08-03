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

import com.google.common.primitives.Ints;
import javastraw.reader.Dataset;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.type.HiCZoom;
import javastraw.reader.type.NormalizationType;
import javastraw.tools.HiCFileTools;
import mixer.clt.CommandLineParserForMixer;
import mixer.clt.MixerCLT;
import mixer.utils.aba.ABADataStack;
import mixer.utils.bed.BedFile;

import java.io.File;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.atomic.AtomicInteger;

public class ABA extends MixerCLT {
    private final Object key = new Object();
    private final boolean saveAllData = true;
    private String bedListPath;
    private File outputDirectory;
    private Dataset ds;
    private int window = 20;
    private int[] resolutions = new int[]{100};


    /**
     * Usage for ABA
     */
    public ABA() {
        super("aba [-w window] [-r resolution(s)] [-c chromosomes]" +
                " [-k NONE/VC/VC_SQRT/KR] [--threads num_of_threads]" +
                " <hicFile(s)> <BedFile> <SaveFolder>");
    }

    public static String getBasicUsage() {
        return "aba <hicFile(s)> <BedFile> <SaveFolder>";
    }

    @Override
    protected void readMixerArguments(String[] args, CommandLineParserForMixer parser) {
        if (args.length != 4) {
            printUsageAndExit(8);
        }

        bedListPath = args[2];
        outputDirectory = HiCFileTools.createValidDirectory(args[3]);

        List<String> summedHiCFiles = Arrays.asList(args[1].split("\\+"));
        ds = HiCFileTools.extractDatasetForCLT(summedHiCFiles.get(0), false, false, false);

        NormalizationType preferredNorm = parser.getNormalizationTypeOption(ds.getNormalizationHandler());
        if (preferredNorm != null)
            norm = preferredNorm;

        int potentialWindow = parser.getWindowSizeOption();
        if (potentialWindow > 0)
            window = potentialWindow;

        List<Integer> possibleResolutions = parser.getMultipleResolutionOptions();
        if (possibleResolutions != null && possibleResolutions.size() > 0) {
            resolutions = Ints.toArray(possibleResolutions);
        }

        updateNumberOfCPUThreads(parser);
    }

    @Override
    public void run() {

        int L = 2 * window + 1;
        for (final int resolution : HiCFileTools.filterResolutions(ds.getBpZooms(), resolutions)) {

            AtomicInteger numOfRegions = new AtomicInteger(0);

            System.out.println("Processing for resolution " + resolution);
            HiCZoom zoom = new HiCZoom(HiCZoom.HiCUnit.BP, resolution);

            ChromosomeHandler handler = ds.getChromosomeHandler();
            if (givenChromosomes != null)
                handler = HiCFileTools.stringToChromosomes(givenChromosomes, handler);

            // Metrics resulting from apa filtering
            final Map<String, Integer[]> filterMetrics = new HashMap<>();

            BedFile bedFile = new BedFile(bedListPath, handler);
            if (bedFile.getNumTotalFeatures() < 1) {
                System.err.println("bed file is empty or incorrect path provided.");
                System.exit(3);
            }


            final AtomicInteger currentProgressStatus = new AtomicInteger(0);
            Chromosome[] chromosomes = handler.getChromosomeArrayWithoutAllByAll();
            double maxProgressStatus = chromosomes.length;

            ExecutorService executor = Executors.newFixedThreadPool(numCPUThreads);
            ABADataStack.initializeDataSaveFolder(outputDirectory, "" + resolution);

            AtomicInteger counter = new AtomicInteger(0);
            for (int l = 0; l < numCPUThreads; l++) {
                executor.execute(() -> {
                    int index = counter.getAndIncrement();
                    while (index < chromosomes.length) {
                        Chromosome chrom = chromosomes[index];
                        MatrixZoomData zd = HiCFileTools.getMatrixZoomData(ds, chrom, chrom, zoom);
                        if (zd != null) {

                            bedFile.doABAOnChrom(chrom, zd, outputDirectory, resolution, numOfRegions, L, window, bedFile, norm);
                            System.out.print(((int) Math.floor((100.0 * currentProgressStatus.incrementAndGet()) / maxProgressStatus)) + "% ");
                        }
                        index = counter.getAndIncrement();
                    }
                });
            }

            executor.shutdown();

            // Wait until all threads finish
            while (!executor.isTerminated()) {
            }

            System.out.println("Exporting ABA results...");
            //save data as int array
            ABADataStack.exportGenomeWideData(numOfRegions.get(), saveAllData);
            ABADataStack.clearAllData();
        }
        System.out.println("ABA complete");
    }
}

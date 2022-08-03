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
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.norm.NormalizationPicker;
import javastraw.reader.type.NormalizationHandler;
import javastraw.reader.type.NormalizationType;
import javastraw.tools.HiCFileTools;
import mixer.clt.CommandLineParserForMixer;
import mixer.clt.MixerCLT;
import mixer.utils.InterChromosomeRegion;
import mixer.utils.bed.BedFileMappings;
import mixer.utils.magic.ClusteringMagic;
import mixer.utils.magic.MagicMatrix;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public class Magic extends MixerCLT {

    public static final boolean PROJECT_TO_UMAP = true;
    private final Random generator = new Random(22871L);
    private Dataset ds;
    private int resolution = 1000;
    private File outputDirectory;
    private String prefix, bedpath;
    private NormalizationType norm = NormalizationHandler.NONE;
    private boolean doScale = false;
    private boolean useZScore = false;

    // subcompartment lanscape identification via clustering enrichment
    public Magic(String command) {
        super("magic [-r resolution] [--verbose] [-k norm] " +
                //"<-k NONE/VC/VC_SQRT/KR/SCALE> [--compare reference.bed] [--has-translocation] " +
                "<file.hic> <K0,KF,nK> <initial_clusters.bed> <outfolder> <prefix_>\n" +
                "   K0 - minimum number of clusters\n" +
                "   KF - maximum number of clusters\n" +
                "   nK - number of times to rerun k-medians");
    }

    @Override
    protected void readMixerArguments(String[] args, CommandLineParserForMixer mixerParser) {
        if (args.length != 6) {
            printUsageAndExit(5);
        }

        ds = HiCFileTools.extractDatasetForCLT(args[1], false, false, false);
        try {
            String[] valString = args[2].split(",");
            ClusteringMagic.startingClusterSizeK = Integer.parseInt(valString[0]);
            ClusteringMagic.numClusterSizeKValsUsed = Integer.parseInt(valString[1]) - ClusteringMagic.startingClusterSizeK;
            ClusteringMagic.numAttemptsForKMeans = Integer.parseInt(valString[2]);
        } catch (Exception e) {
            printUsageAndExit(5);
        }

        bedpath = args[3];
        outputDirectory = HiCFileTools.createValidDirectory(args[4]);
        prefix = args[5];

        NormalizationType potentialNorm = mixerParser.getNormalizationTypeOption(ds.getNormalizationHandler());
        if (potentialNorm != null) {
            norm = potentialNorm;
        }

        List<Integer> possibleResolutions = mixerParser.getMultipleResolutionOptions();
        if (possibleResolutions != null) {
            if (possibleResolutions.size() > 1)
                System.err.println("Only one resolution can be specified\nUsing " + possibleResolutions.get(0));
            resolution = possibleResolutions.get(0);
        }

        doScale = mixerParser.getScaleOption();
        useZScore = mixerParser.getZScoreOption();

        long[] possibleSeeds = mixerParser.getMultipleSeedsOption();
        if (possibleSeeds != null && possibleSeeds.length > 0) {
            for (long seed : possibleSeeds) {
                generator.setSeed(seed);
            }
        }
    }

    @Override
    public void run() {

        ChromosomeHandler chromosomeHandler = ds.getChromosomeHandler();
        if (givenChromosomes != null)
            chromosomeHandler = HiCFileTools.stringToChromosomes(givenChromosomes, chromosomeHandler);

        List<InterChromosomeRegion> regionsToIgnore = new ArrayList<>();

        NormalizationType validNormForFiltering = NormalizationPicker.getFirstValidNormInThisOrder(ds,
                new String[]{"SCALE", "KR", "VC", "VC_SQRT"});

        BedFileMappings mappings = new BedFileMappings(bedpath, chromosomeHandler, resolution, ds, validNormForFiltering);
        MagicMatrix matrix = new MagicMatrix(ds, chromosomeHandler, resolution, norm,
                outputDirectory, generator.nextLong(), mappings, regionsToIgnore, doScale, useZScore);
        matrix.export(new File(outputDirectory, "magic.npy").getAbsolutePath());

        ClusteringMagic clustering = new ClusteringMagic(matrix, outputDirectory, chromosomeHandler, 10L);
        clustering.extractFinalGWSubcompartments(prefix);
        System.out.println("\nClustering complete");
    }
}

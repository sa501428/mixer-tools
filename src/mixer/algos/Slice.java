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
import javastraw.reader.type.NormalizationType;
import javastraw.tools.HiCFileTools;
import mixer.clt.CommandLineParserForMixer;
import mixer.clt.MixerCLT;
import mixer.utils.slice.cleaning.SliceMatrixCleaner;
import mixer.utils.slice.kmeans.FullGenomeOEWithinClusters;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

/**
 * experimental code
 * <p>
 * Created by muhammadsaadshamim on 9/14/15.
 */
public class Slice extends MixerCLT {

    public static final int INTRA_SCALE_INDEX = 0;
    public static final int INTER_SCALE_INDEX = 1;
    public static final int GW_SCALE_INDEX = 2;
    public static final boolean USE_INTER_CORR_CLUSTERING = false;
    public static final boolean PROJECT_TO_UMAP = true;
    public static final boolean USE_WEIGHTED_MEAN = false;
    public static boolean FILTER_OUTLIERS = true;
    private final List<Dataset> datasetList = new ArrayList<>();
    private final List<String> inputHicFilePaths = new ArrayList<>();
    private final Random generator = new Random(22871L);
    private int resolution = 100000;
    private Dataset ds;
    private File outputDirectory;
    private List<NormalizationType[]> normsList;
    private String prefix = "";
    public static boolean USE_KMEANS = false, USE_KMEDIANS = true;
    public static boolean USE_ENCODE_MODE = false;

    // subcompartment lanscape identification via clustering enrichment
    public Slice(String command) {
        super("slice [-r resolution] [--verbose] " +
                //"<-k NONE/VC/VC_SQRT/KR/SCALE> [--compare reference.bed] [--has-translocation] " +
                "<file.hic> <K0,KF,nK> <outfolder> <prefix_>\n" +
                "   K0 - minimum number of clusters\n" +
                "   KF - maximum number of clusters\n" +
                "   nK - number of times to rerun kmeans");
    }

    @Override
    protected void readMixerArguments(String[] args, CommandLineParserForMixer mixerParser) {
        if (args.length != 5) {
            printUsageAndExit(5);
        }

        for (String path : args[1].split(",")) {
            System.out.println("Extracting " + path);
            inputHicFilePaths.add(path);
            datasetList.add(HiCFileTools.extractDatasetForCLT(path, false, false, true));
        }

        try {
            String[] valString = args[2].split(",");
            FullGenomeOEWithinClusters.startingClusterSizeK = Integer.parseInt(valString[0]);
            FullGenomeOEWithinClusters.numClusterSizeKValsUsed = Integer.parseInt(valString[1]) -
                    FullGenomeOEWithinClusters.startingClusterSizeK;
            FullGenomeOEWithinClusters.numAttemptsForKMeans = Integer.parseInt(valString[2]);
        } catch (Exception e) {
            printUsageAndExit(5);
        }

        ds = datasetList.get(0);
        outputDirectory = HiCFileTools.createValidDirectory(args[3]);
        prefix = args[4];
        normsList = populateNormalizations(datasetList);

        List<Integer> possibleResolutions = mixerParser.getMultipleResolutionOptions();
        if (possibleResolutions != null) {
            if (possibleResolutions.size() > 1)
                System.err.println("Only one resolution can be specified\nUsing " + possibleResolutions.get(0));
            resolution = possibleResolutions.get(0);
        }

        long[] possibleSeeds = mixerParser.getMultipleSeedsOption();
        if (possibleSeeds != null && possibleSeeds.length > 0) {
            for (long seed : possibleSeeds) {
                generator.setSeed(seed);
            }
        }

        int subsampling = mixerParser.getSubsamplingOption();
        if (subsampling > 0) {
            SliceMatrixCleaner.NUM_PER_CENTROID = subsampling;
        }

        USE_ENCODE_MODE = mixerParser.getENCODEOption();
    }


    private List<NormalizationType[]> populateNormalizations(List<Dataset> datasetList) {
        List<NormalizationType[]> normsList = new ArrayList<>();
        for (Dataset ds : datasetList) {
            NormalizationType[] norms = new NormalizationType[3];
            norms[INTRA_SCALE_INDEX] = NormalizationPicker.getFirstValidNormInThisOrder(ds, new String[]{"SCALE", "KR"});
            norms[INTER_SCALE_INDEX] = NormalizationPicker.getFirstValidNormInThisOrder(ds, new String[]{"INTER_SCALE", "INTER_KR"});
            norms[GW_SCALE_INDEX] = NormalizationPicker.getFirstValidNormInThisOrder(ds, new String[]{"GW_SCALE", "GW_KR"});
            normsList.add(norms);
        }
        return normsList;
    }

    @Override
    public void run() {

        ChromosomeHandler chromosomeHandler = ds.getChromosomeHandler();
        if (givenChromosomes != null)
            chromosomeHandler = HiCFileTools.stringToChromosomes(givenChromosomes, chromosomeHandler);

        if (datasetList.size() < 1) return;

        FullGenomeOEWithinClusters withinClusters = new FullGenomeOEWithinClusters(datasetList,
                chromosomeHandler, resolution, normsList, outputDirectory, generator.nextLong());
        withinClusters.extractFinalGWSubcompartments(prefix);

        System.out.println("\nClustering complete");
    }
}

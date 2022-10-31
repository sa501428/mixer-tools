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

import javastraw.feature1D.GenomeWide1DList;
import javastraw.reader.Dataset;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.type.NormalizationType;
import javastraw.tools.HiCFileTools;
import javastraw.tools.UNIXTools;
import mixer.clt.CommandLineParserForMixer;
import mixer.clt.MixerCLT;
import mixer.utils.matrix.InterOnlyMatrix;
import mixer.utils.shuffle.Partition;
import mixer.utils.shuffle.ShuffleAction;
import mixer.utils.similaritymeasures.SimilarityMetric;
import mixer.utils.slice.structures.SliceUtils;
import mixer.utils.slice.structures.SubcompartmentInterval;
import mixer.utils.umap.UMAPAction;

import java.io.File;
import java.util.List;
import java.util.Random;

/**
 * experimental code
 * <p>
 * Created by muhammadsaadshamim on 9/14/15.
 */
public class ShuffleAndUMAP extends MixerCLT {

    //private final Random generator = new Random(22871L);
    private long randomSeed = 22871L;
    private Dataset ds;
    private int resolution = 100000;
    private int compressionFactor = 8;
    private File outputDirectory;
    private String[] prefix;
    private String[] referenceBedFiles;
    private final boolean useIntraMap, useGWMap;
    private InterOnlyMatrix.INTRA_TYPE intraType;
    private SimilarityMetric metric = null;
    private final boolean doUMAP, doShuffle, useSymmetry, useOriginal;

    // subcompartment lanscape identification via clustering enrichment
    public ShuffleAndUMAP(String name) {
        super("[intra-]shuffle-umap [-r resolution] [-k NONE/VC/VC_SQRT/KR/SCALE] [-w window] [--verbose] " +
                "<file.hic> <subcompartment.bed(s)> <outfolder> <prefix>");
        useGWMap = name.contains("gw");
        useIntraMap = name.contains("intra");
        doShuffle = name.contains("shuffle");
        useOriginal = name.contains("original");
        doUMAP = name.contains("umap");
        useSymmetry = name.contains("symm");
    }

    @Override
    protected void readMixerArguments(String[] args, CommandLineParserForMixer mixerParser) {
        if (args.length != 5) {
            printUsageAndExit(51);
        }

        ds = HiCFileTools.extractDatasetForCLT(args[1], false, false, true);

        referenceBedFiles = args[2].split(",");
        outputDirectory = HiCFileTools.createValidDirectory(args[3]);
        prefix = args[4].split(",");

        NormalizationType preferredNorm = mixerParser.getNormalizationTypeOption(ds.getNormalizationHandler());
        if (preferredNorm != null) norm = preferredNorm;

        List<Integer> possibleResolutions = mixerParser.getMultipleResolutionOptions();
        if (possibleResolutions != null) {
            if (possibleResolutions.size() > 1)
                System.err.println("Only one resolution can be specified\nUsing " + possibleResolutions.get(0));
            resolution = possibleResolutions.get(0);
        }

        intraType = InterOnlyMatrix.getIntraType(mixerParser.getMapTypeOption());
        metric = SimilarityMetric.getMetric(mixerParser.getCorrelationTypeOption());

        randomSeed = mixerParser.getSeedOption();

        int minSize = mixerParser.getWindowSizeOption();
        if (minSize > 1) {
            compressionFactor = minSize;
        }
    }

    @Override
    public void run() {

        ChromosomeHandler chromosomeHandler = ds.getChromosomeHandler();
        if (givenChromosomes != null)
            chromosomeHandler = HiCFileTools.stringToChromosomes(givenChromosomes, chromosomeHandler);

        UNIXTools.makeDir(outputDirectory);

        Partition.Type[] mapTypes = {Partition.Type.SPLIT1, Partition.Type.SPLIT2,
                Partition.Type.SPLIT3, Partition.Type.SPLIT4, Partition.Type.SPLIT5};
        if (useOriginal) {
            mapTypes = new Partition.Type[]{Partition.Type.ODDS_VS_EVENS,
                    Partition.Type.SKIP_BY_TWOS, Partition.Type.FIRST_HALF_VS_SECOND_HALF};
        }

        if (doShuffle) {
            for (int i = 0; i < referenceBedFiles.length; i++) {
                GenomeWide1DList<SubcompartmentInterval> subcompartments =
                        SliceUtils.loadFromSubcompartmentBEDFile(chromosomeHandler, referenceBedFiles[i]);
                System.out.println("Processing " + prefix[i]);
                File newFolder = new File(outputDirectory, "shuffle_" + prefix[i]);
                UNIXTools.makeDir(newFolder);
                Random generator = new Random(randomSeed);
                ShuffleAction matrix;
                if (useIntraMap) {
                    matrix = new ShuffleAction(ds, norm, resolution, compressionFactor,
                            metric, intraType, useSymmetry, mapTypes);
                    matrix.runIntraAnalysis(subcompartments, newFolder, generator);
                } else {
                    matrix = new ShuffleAction(ds, norm, resolution, compressionFactor, metric, useSymmetry, mapTypes);
                    matrix.runInterAnalysis(subcompartments, newFolder, generator);
                }
                matrix.savePlotsAndResults(newFolder, prefix[i]);
            }
            System.out.println("Shuffle complete");
        }

        if (doUMAP) {
            UMAPAction umap;
            if (useIntraMap) {
                umap = new UMAPAction(ds, norm, resolution, compressionFactor, metric, intraType, mapTypes);
            } else {
                umap = new UMAPAction(ds, norm, resolution, compressionFactor, metric, useGWMap, mapTypes);
            }
            umap.runAnalysis(referenceBedFiles, outputDirectory, chromosomeHandler);
        }
    }


    private Chromosome[] getStringToChromosomes(List<String> givenChromosomes, ChromosomeHandler chromosomeHandler) {
        Chromosome[] result = new Chromosome[givenChromosomes.size()];
        for (int i = 0; i < result.length; i++) {
            String chrom = givenChromosomes.get(i);
            result[i] = chromosomeHandler.getChromosomeFromName(chrom);
        }
        return result;
    }
}

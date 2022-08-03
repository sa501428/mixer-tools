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
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.tools.HiCFileTools;
import javastraw.tools.MatrixTools;
import mixer.clt.CommandLineParserForMixer;
import mixer.clt.MixerCLT;
import mixer.utils.network.ConnectedComponentBFS;
import mixer.utils.network.LoopAnchor;

import java.io.File;
import java.util.*;

/**
 * Created by muhammadsaadshamim on 9/14/15.
 */
public class Network extends MixerCLT {

    private final int bufferWidth = 3;
    private String inputBedpeFile, outputPath;
    private int resolution = 5000;
    private Dataset ds;

    public Network() {
        super("network [-r resolution] <hicfile> <bedpe> <output_stem>");
    }

    @Override
    protected void readMixerArguments(String[] args, CommandLineParserForMixer mixerParser) {
        if (args.length != 4) {
            printUsageAndExit(6);
        }

        inputBedpeFile = args[2];
        outputPath = args[3];
        ds = HiCFileTools.extractDatasetForCLT(args[1], false, false, false);

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

        List<Integer> networkSizes = identifyNetworks(listA);
        int[] sizes = networkSizes.stream().mapToInt(i -> i).toArray();

        MatrixTools.saveMatrixTextNumpy(outputPath + "loop_network_sizes.npy", sizes);
        listA.exportFeatureList(new File(outputPath + "_large_networks.bedpe"), false, Feature2DList.ListFormat.NA);
    }

    private List<Integer> identifyNetworks(Feature2DList list) {
        List<Integer> loopNetworkSizes = Collections.synchronizedList(new ArrayList<>());

        list.filterLists((chr, feature2DList) -> {
            ConnectedComponentBFS cc = processNetworksForChromosome(resolution, feature2DList);

            loopNetworkSizes.addAll(cc.getClusterSizes());

            return extractBigNetworks(cc, feature2DList);
        });

        return loopNetworkSizes;
    }

    private List<Feature2D> extractBigNetworks(ConnectedComponentBFS cc, List<Feature2D> feature2DList) {
        List<Feature2D> bigNetworks = new ArrayList<>();
        for (Set<Integer> network : cc.getLargeNetworks()) {
            for (int id : network) {
                bigNetworks.add(feature2DList.get(id));
            }
        }
        return bigNetworks;
    }

    private ConnectedComponentBFS processNetworksForChromosome(int resolution, List<Feature2D> feature2DList) {
        List<LoopAnchor> upStreamList = new ArrayList<>();
        List<LoopAnchor> downStreamList = new ArrayList<>();
        //Matrix matrix = ds.getMatrix(chromosome, chromosome);
        //MatrixZoomData zd = matrix.getZoomData(new HiCZoom(HiC.Unit.BP, resolution));

        for (int i = 0; i < feature2DList.size(); i++) {
            Feature2D loop = feature2DList.get(i);
            long mid1 = loop.getMidPt1() / resolution;
            long mid2 = loop.getMidPt2() / resolution;

            handleAnchor(mid1, upStreamList, i);
            handleAnchor(mid2, downStreamList, i);
        }
        Map<Integer, Set<Integer>> adjacencyMatrix = generateAdjacencies(upStreamList, downStreamList);
        ConnectedComponentBFS cc = new ConnectedComponentBFS(adjacencyMatrix, feature2DList.size());
        return cc;
    }

    private Map<Integer, Set<Integer>> generateAdjacencies(List<LoopAnchor> upStreamList, List<LoopAnchor> downStreamList) {
        Map<Integer, Set<Integer>> adjacencyMatrix = new HashMap<>();
        addContactsToMap(adjacencyMatrix, upStreamList);
        addContactsToMap(adjacencyMatrix, downStreamList);
        return adjacencyMatrix;
    }

    private void addContactsToMap(Map<Integer, Set<Integer>> adjacencyMatrix, List<LoopAnchor> anchors) {
        for (LoopAnchor anchor : anchors) {
            if (anchor.size() > 1) {
                anchor.populateAdjacencyMatrix(adjacencyMatrix);
            }
        }
    }

    private void handleAnchor(long midBin, List<LoopAnchor> anchors, Integer id) {
        boolean notHandled = true;
        for (LoopAnchor anchor : anchors) {
            if (anchor.overlaps(midBin)) {
                notHandled = false;
                anchor.add(midBin, id);
            }
        }
        if (notHandled) {
            anchors.add(new LoopAnchor(midBin, id, bufferWidth));
        }
    }
}

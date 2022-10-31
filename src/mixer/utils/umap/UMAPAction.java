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

import javastraw.feature1D.GenomeWide1DList;
import javastraw.reader.Dataset;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.type.NormalizationType;
import javastraw.tools.MatrixTools;
import javastraw.tools.UNIXTools;
import mixer.utils.common.FloatMatrixTools;
import mixer.utils.matrix.GWInterOnlyMatrix;
import mixer.utils.matrix.HiCMatrix;
import mixer.utils.matrix.InterOnlyMatrix;
import mixer.utils.matrix.IntraOnlyMatrix;
import mixer.utils.shuffle.Partition;
import mixer.utils.similaritymeasures.SimilarityMetric;
import mixer.utils.slice.cleaning.IntraMatrixCleaner;
import mixer.utils.slice.structures.SliceUtils;
import mixer.utils.slice.structures.SubcompartmentInterval;
import tagbio.umap.Umap;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public class UMAPAction {

    private final Random generator = new Random(0);
    private final Dataset ds;
    private final NormalizationType norm;
    private final int compressionFactor;
    private final int resolution;
    private final SimilarityMetric metric;
    private final boolean useGWMap;
    private final Partition.Type[] mapTypes;
    private HiCMatrix.INTRA_TYPE intraType = null;
    private boolean isIntra = false;

    public UMAPAction(Dataset ds, NormalizationType norm, int resolution, int compressionFactor,
                      SimilarityMetric metric, boolean useGWMap, Partition.Type[] mapTypes) {
        this.resolution = resolution;
        this.compressionFactor = compressionFactor;
        this.ds = ds;
        this.norm = norm;
        this.metric = metric;
        this.useGWMap = useGWMap;
        this.mapTypes = mapTypes;
    }

    public UMAPAction(Dataset ds, NormalizationType norm, int resolution, int compressionFactor,
                      SimilarityMetric metric, HiCMatrix.INTRA_TYPE intraType, Partition.Type[] mapTypes) {
        this(ds, norm, resolution, compressionFactor, metric, false, mapTypes);
        this.intraType = intraType;
        isIntra = true;
    }

    public void runAnalysis(String[] bedFiles, File outputDirectory,
                            ChromosomeHandler chromosomeHandler) {

        List<GenomeWide1DList<SubcompartmentInterval>> allSubcompartments = new ArrayList<>();
        for (int i = 0; i < bedFiles.length; i++) {
            GenomeWide1DList<SubcompartmentInterval> subcompartments =
                    SliceUtils.loadFromSubcompartmentBEDFile(chromosomeHandler, bedFiles[i]);
            SliceUtils.collapseGWList(subcompartments);
            allSubcompartments.add(subcompartments);
        }

        try {
            writeToFile(outputDirectory, bedFiles);
        } catch (Exception e) {
            e.printStackTrace();
        }

        int numBedFiles = allSubcompartments.size();
        // todo, using only big size? todo sorting picture

        if (useGWMap) {
            final HiCMatrix interMatrix = GWInterOnlyMatrix.getMatrix(ds, norm, resolution, metric, chromosomeHandler, compressionFactor);
            float[][] matrix = interMatrix.getMatrix();

            FloatMatrixTools.cleanUpMatrix(matrix, true);
            System.out.println("Saving files as png; dimensions " + matrix.length + " " + matrix[0].length);
            File file = new File(outputDirectory, "matrix.png");
            FloatMatrixTools.saveMatrixToPNG(file, matrix, false);
            //file = new File(outputDirectory, "matrix_log.png");
            //FloatMatrixTools.saveMatrixToPNG(file, matrix, true);

            System.out.println("GW UMAP Run");

            processMatrixForUMAP(matrix, numBedFiles,
                    interMatrix.getRowChromosomes(), interMatrix.getRowOffsets(),
                    allSubcompartments, outputDirectory, "genomewide", "rows");
        } else if (isIntra) {
            Chromosome[] chromosomes = chromosomeHandler.getAutosomalChromosomesArray();
            for (int y = 0; y < chromosomes.length; y++) {
                final HiCMatrix interMatrix = new IntraOnlyMatrix(ds, norm, resolution, chromosomes[y],
                        intraType, metric, compressionFactor);

                float[][] matrix = interMatrix.getMatrix();

                processMatrixForUMAP(matrix, numBedFiles,
                        interMatrix.getRowChromosomes(), interMatrix.getRowOffsets(),
                        allSubcompartments, outputDirectory, chromosomes[y].getName(), "rows");
            }
        } else {
            for (int y = 0; y < mapTypes.length; y++) {

                final HiCMatrix interMatrix = InterOnlyMatrix.getMatrix(ds, norm, resolution, mapTypes[y], metric);
                interMatrix.applySimpleLog();

                float[][] matrix = interMatrix.getMatrix();
                float[][] matrixT = FloatMatrixTools.transpose(matrix);

                IntraMatrixCleaner.rollingAverage(matrix, compressionFactor);
                IntraMatrixCleaner.rollingAverage(matrixT, compressionFactor);

                processMatrixForUMAP(matrix, numBedFiles,
                        interMatrix.getRowChromosomes(), interMatrix.getRowOffsets(),
                        allSubcompartments, outputDirectory, mapTypes[y].toString(), "rows");

                processMatrixForUMAP(matrixT, numBedFiles,
                        interMatrix.getColChromosomes(), interMatrix.getColOffsets(),
                        allSubcompartments, outputDirectory, mapTypes[y].toString(), "cols");
            }
        }
    }

    private void processMatrixForUMAP(float[][] matrix, int numBedFiles, Chromosome[] chromosomes, int[] offsets,
                                      List<GenomeWide1DList<SubcompartmentInterval>> allSubcompartments,
                                      File outputDirectory, String folderName, String filePrefix) {


        int[][] indicesToClusterIDs = new int[matrix.length][numBedFiles];

        for (int z = 0; z < numBedFiles; z++) {
            populateIndexToClusterIDMap(chromosomes, offsets, allSubcompartments.get(z), indicesToClusterIDs, z);
        }

        File outfolder = new File(outputDirectory, folderName);
        UNIXTools.makeDir(outfolder);
        runUmapAndSaveMatrices(matrix, outfolder, filePrefix, indicesToClusterIDs);
    }

    private void runUmapAndSaveMatrices(float[][] initialData, File outputDirectory, String outstem,
                                        int[][] initialIndexToIDs) {

        ZeroRowRemover rowRemover = new ZeroRowRemover(initialData, initialIndexToIDs);

        int numThreads = Math.max(1, Runtime.getRuntime().availableProcessors());
        System.out.println("Running UMAP with " + outstem);
        final Umap umap = new Umap();
        umap.setNumberComponents(2);
        umap.setNumberNearestNeighbours(100); // 50 // 15 -> 50 for more global picture
        umap.setThreads(numThreads);
        umap.setMinDist(0.5f); // 0.2 ->0.8 -> 0.5  //0.1f -> 0.2f for more general features
        umap.setVerbose(false);
        umap.setSeed(generator.nextLong());
        final float[][] result = umap.fitTransform(rowRemover.getCleanData());

        File temp = new File(outputDirectory, "umap_" + outstem + "_embedding.npy");
        FloatMatrixTools.saveMatrixTextNumpy(temp.getAbsolutePath(), result);
        System.gc();

        temp = new File(outputDirectory, outstem + "genome_indices.npy");
        MatrixTools.saveMatrixTextNumpy(temp.getAbsolutePath(), rowRemover.getCleanIndices());

        ScatterPlot sc = new ScatterPlot(result, rowRemover.getCleanIndices());
        temp = new File(outputDirectory, outstem + "_plots");
        sc.plot(temp.getAbsolutePath());
    }

    private void populateIndexToClusterIDMap(Chromosome[] chromosomes, int[] offsets,
                                             GenomeWide1DList<SubcompartmentInterval> subcompartments,
                                             int[][] indicesToClusterIDs, int colIndex) {
        for (int x = 0; x < chromosomes.length; x++) {
            Chromosome chrom = chromosomes[x];
            List<SubcompartmentInterval> intervalList = subcompartments.getFeatures("" + chrom.getIndex());
            for (SubcompartmentInterval interval : intervalList) {
                int xStart = interval.getX1() / resolution;
                int xEnd = interval.getX2() / resolution;
                int clusterID = interval.getClusterID();

                for (int k = xStart; k < xEnd; k++) {
                    final int actualPosition = k + offsets[x];
                    indicesToClusterIDs[actualPosition][colIndex] = clusterID;
                }
            }
        }
    }

    private void writeToFile(File outfolder, String[] filenames) throws IOException {
        FileWriter myWriter = new FileWriter(new File(outfolder, "bedfile_order.txt"));
        for (String filename : filenames) {
            myWriter.write(filename + "\n");
        }
        myWriter.close();
    }
}

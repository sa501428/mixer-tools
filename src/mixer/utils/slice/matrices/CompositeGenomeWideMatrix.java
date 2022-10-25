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

package mixer.utils.slice.matrices;

import javastraw.feature1D.GenomeWide1DList;
import javastraw.reader.Dataset;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.type.NormalizationType;
import javastraw.tools.MatrixTools;
import mixer.MixerGlobals;
import mixer.algos.Slice;
import mixer.utils.common.ArrayTools;
import mixer.utils.common.FloatMatrixTools;
import mixer.utils.common.ZScoreTools;
import mixer.utils.drive.DriveMatrix;
import mixer.utils.slice.cleaning.BadIndexFinder;
import mixer.utils.slice.cleaning.SimilarityMatrixTools;
import mixer.utils.slice.cleaning.SliceMatrixCleaner;
import mixer.utils.slice.structures.SliceUtils;
import mixer.utils.slice.structures.SubcompartmentInterval;

import java.io.File;
import java.util.*;

public abstract class CompositeGenomeWideMatrix extends DriveMatrix {
    protected final NormalizationType[] norms;
    protected final int resolution;
    protected final Map<Integer, SubcompartmentInterval> rowIndexToIntervalMap = new HashMap<>();
    protected final Chromosome[] chromosomes;
    protected final Random generator = new Random(2436);
    protected final File outputDirectory;
    private MatrixAndWeight gwCleanMatrix, projectedData = null;
    protected final BadIndexFinder badIndexLocations;
    protected final int maxClusterSizeExpected;

    public CompositeGenomeWideMatrix(ChromosomeHandler chromosomeHandler, Dataset ds,
                                     NormalizationType[] norms,
                                     int resolution,
                                     File outputDirectory, long seed,
                                     BadIndexFinder badIndexLocations,
                                     int maxClusterSizeExpected) {
        this.maxClusterSizeExpected = maxClusterSizeExpected;
        this.norms = norms;
        if (MixerGlobals.printVerboseComments) {
            System.out.println("Norms: " + Arrays.toString(norms));
        }
        this.resolution = resolution;
        this.outputDirectory = outputDirectory;
        generator.setSeed(seed);
        this.badIndexLocations = badIndexLocations;

        chromosomes = chromosomeHandler.getAutosomalChromosomesArray();
        gwCleanMatrix = makeCleanScaledInterMatrix(ds, norms[Slice.INTER_SCALE_INDEX]);
        //gwCleanMatrix = makeComboNormMatrix(ds, norms);
    }

    private MatrixAndWeight makeComboNormMatrix(Dataset ds, NormalizationType[] norms) {
        return concatenate(
                makeCleanScaledInterMatrix(ds, norms[Slice.INTER_SCALE_INDEX]),
                makeCleanScaledInterMatrix(ds, norms[Slice.GW_SCALE_INDEX])
        );
    }

    abstract MatrixAndWeight makeCleanScaledInterMatrix(Dataset ds, NormalizationType interNorm);

    public void cleanUpMatricesBySparsity() {

        SliceMatrixCleaner matrixCleanupReduction = new SliceMatrixCleaner(gwCleanMatrix.matrix,
                generator.nextLong(), outputDirectory, resolution);
        gwCleanMatrix = matrixCleanupReduction.getCleanFilteredZscoredMatrix(rowIndexToIntervalMap,
                gwCleanMatrix.weights);

        inPlaceScaleSqrtWeightCol();

        File file0 = new File(outputDirectory, "genome_indices.npy");
        MatrixTools.saveMatrixTextNumpy(file0.getAbsolutePath(), getGenomeIndices());

        File file1 = new File(outputDirectory, "slice_matrix.npy");
        MatrixTools.saveMatrixTextNumpy(file1.getAbsolutePath(), gwCleanMatrix.matrix);

        if (Slice.USE_INTER_CORR_CLUSTERING || Slice.PROJECT_TO_UMAP) {
            projectedData = new MatrixAndWeight(SimilarityMatrixTools.getCosinePearsonCorrMatrix(gwCleanMatrix.matrix,
                    50, generator.nextLong()), gwCleanMatrix.weights);


            File file2 = new File(outputDirectory, "corr_slice_matrix.npy");
            MatrixTools.saveMatrixTextNumpy(file2.getAbsolutePath(), projectedData.matrix);
        }
    }

    protected int[][] getGenomeIndices() {
        int n = gwCleanMatrix.matrix.length;
        int[][] coordinates = new int[n][3];
        for (int i = 0; i < n; i++) {
            SubcompartmentInterval interval = rowIndexToIntervalMap.get(i);
            coordinates[i][0] = interval.getChrIndex();
            coordinates[i][1] = interval.getX1();
            coordinates[i][2] = interval.getX2();
        }
        return coordinates;
    }

    @Override
    public void inPlaceScaleSqrtWeightCol() {
        ZScoreTools.inPlaceScaleSqrtWeightCol(gwCleanMatrix.matrix, gwCleanMatrix.weights);
    }

    public synchronized void processGMMClusteringResult(int[] clusterID,
                                                        GenomeWide1DList<SubcompartmentInterval> subcompartments) {

        Set<SubcompartmentInterval> subcompartmentIntervals = new HashSet<>();

        for (int i = 0; i < projectedData.matrix.length; i++) {
            if (rowIndexToIntervalMap.containsKey(i)) {
                SubcompartmentInterval interv = rowIndexToIntervalMap.get(i);
                if (interv != null) {
                    subcompartmentIntervals.add(generateNewSubcompartment(interv, clusterID[i]));
                }
            }
        }

        subcompartments.addAll(new ArrayList<>(subcompartmentIntervals));
        SliceUtils.reSort(subcompartments);
    }

    public float[][] getData(boolean getCorrMatrix) {
        if (getCorrMatrix) {
            return projectedData.matrix;
        }
        return gwCleanMatrix.matrix;
    }

    public void appendDataAlongExistingRows(CompositeGenomeWideMatrix additionalData) {
        if (gwCleanMatrix.matrix.length != additionalData.gwCleanMatrix.matrix.length) {
            System.err.println("***************************************\n" +
                    "Dimension mismatch: " + gwCleanMatrix.matrix.length + " != " + additionalData.gwCleanMatrix.matrix.length);
        } else {
            gwCleanMatrix = concatenate(gwCleanMatrix, additionalData.gwCleanMatrix);
        }
    }

    private MatrixAndWeight concatenate(MatrixAndWeight mw1, MatrixAndWeight mw2) {
        return new MatrixAndWeight(
                FloatMatrixTools.concatenate(mw1.matrix, mw2.matrix),
                ArrayTools.concatenate(mw1.weights, mw2.weights)
        );
    }

    public BadIndexFinder getBadIndices() {
        return badIndexLocations;
    }

    public Map<Integer, SubcompartmentInterval> getRowIndexToIntervalMap() {
        return rowIndexToIntervalMap;
    }
}

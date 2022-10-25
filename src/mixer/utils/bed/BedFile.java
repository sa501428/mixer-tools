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

package mixer.utils.bed;

import javastraw.feature1D.GenomeWide1DList;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.type.NormalizationType;
import javastraw.tools.HiCFileTools;
import mixer.MixerGlobals;
import mixer.utils.aba.ABADataStack;
import mixer.utils.slice.structures.SimpleInterval;
import org.apache.commons.math3.linear.RealMatrix;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.regex.Pattern;

public class BedFile {
    private final GenomeWide1DList<SimpleInterval> featureGenomeWide1DList;

    public BedFile(String bedListPath, ChromosomeHandler handler) {
        featureGenomeWide1DList = populateBedFile(handler, bedListPath);
    }

    public static RealMatrix extractLocalizedData(MatrixZoomData zd, SimpleInterval region, int overallWidth,
                                                  int resolution, int window, NormalizationType norm) throws IOException {
        long loopX = region.getX1() / resolution;
        long binXStart = loopX - window;
        long binXEnd = loopX + (window + 1);

        return HiCFileTools.extractLocalBoundedRegion(zd, binXStart, binXEnd, binXStart, binXEnd, overallWidth, overallWidth,
                norm, true);
    }

    public static GenomeWide1DList<SimpleInterval> populateBedFile(ChromosomeHandler handler, String bedFilePath) {
        List<SimpleInterval> anchors = new ArrayList<>();

        try {
            //BufferedReader br = ParsingUtils.openBufferedReader(bedFilePath);
            anchors.addAll(parseBEDFile(bedFilePath, handler));
        } catch (IOException ec) {
            ec.printStackTrace();
        }

        return new GenomeWide1DList<>(handler, anchors);
    }

    private static List<SimpleInterval> parseBEDFile(String bedFilePath, ChromosomeHandler handler) throws IOException {
        BufferedReader bufferedReader = new BufferedReader(new FileReader(bedFilePath), MixerGlobals.bufferSize);

        Set<SimpleInterval> positions = new HashSet<>();
        String nextLine;

        int errorCount = 0;
        while ((nextLine = bufferedReader.readLine()) != null) {
            String[] tokens = Pattern.compile("\t").split(nextLine);

            if (nextLine.startsWith("#")) {
                continue;
            }

            if (tokens.length > 2) {
                // valid line
                String chr1Name = tokens[0];
                int start1 = Integer.parseInt(tokens[1]);
                int end1 = Integer.parseInt(tokens[2]);

                Chromosome chr = handler.getChromosomeFromName(chr1Name);
                if (chr == null) {
                    if (errorCount < 3) {
                        System.err.println("Skipping line: " + nextLine);
                    } else if (errorCount == 3) {
                        System.err.println("Maximum error count exceeded.  Further errors will not be logged");
                        System.err.println("If chromosomes were intentionally filtered, ignore this error.");
                    }

                    errorCount++;
                    continue;
                }

                positions.add(new SimpleInterval(chr.getIndex(), chr.getName(), start1, end1));
            }
        }
        if (positions.size() < 1)
            System.err.println("BED File empty - file may have problems or error was encountered");
        bufferedReader.close();
        return new ArrayList<>(positions);
    }

    public void doABAOnChrom(Chromosome chrom, MatrixZoomData zd, File outputDirectory, int resolution,
                             AtomicInteger numOfRegions, int overallWidth, int window, BedFile bedFile,
                             NormalizationType norm) {
        ABADataStack dataStack = new ABADataStack(overallWidth, outputDirectory, "" + resolution);
        if (MixerGlobals.printVerboseComments) {
            System.out.println("CHR " + chrom.getName() + " " + chrom.getIndex());
        }

        List<SimpleInterval> regions = bedFile.get(chrom.getIndex());
        if (regions == null || regions.size() == 0) {
            if (MixerGlobals.printVerboseComments) {
                System.out.println("Chromosome " + chrom.getName() + " - no regions found");
            }
            return;
        }

        numOfRegions.addAndGet(regions.size());

        for (SimpleInterval region : regions) {
            try {
                RealMatrix newData = extractLocalizedData(zd, region, overallWidth, resolution, window, norm);
                dataStack.addData(newData);
                //dataStack.addData(APAUtils.extractLocalizedData(zd, loop, L, resolution, window, norm));
            } catch (Exception e) {
                System.err.println(e.getMessage());
                System.err.println("Unable to find data for loop: " + region);
            }
        }

        dataStack.updateGenomeWideData();
        //dataStack.exportDataSet(chrom.getName() + 'v' + chr2.getName(), peakNumbers);

    }

    private List<SimpleInterval> get(int index) {
        return featureGenomeWide1DList.getFeatures("" + index);
    }

    public int getNumTotalFeatures() {
        return featureGenomeWide1DList.size();
    }
}

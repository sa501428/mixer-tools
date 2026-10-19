/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2011-2021 Rice University, Baylor College of Medicine, Aiden Lab
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

package mixer.utils.slice.structures;

import javastraw.feature1D.GenomeWide1DList;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ChromosomeHandler;
import mixer.MixerGlobals;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.regex.Pattern;

public class SliceUtils {

    private static final AtomicInteger idCounter = new AtomicInteger(0);

    public static void reSort(GenomeWide1DList<SubcompartmentInterval> subcompartments) {
        subcompartments.filterLists((chr, featureList) -> {
            Collections.sort(featureList);
            return featureList;
        });
    }

    public static void collapseGWList(GenomeWide1DList<SubcompartmentInterval> intraSubcompartments) {
        intraSubcompartments.filterLists((chr, featureList) -> collapseSubcompartmentIntervals(featureList));
    }

    public static List<SubcompartmentInterval> collapseSubcompartmentIntervals(List<SubcompartmentInterval> intervals) {
        if (intervals.size() > 0) {

            Collections.sort(intervals);
            SubcompartmentInterval collapsedInterval = (SubcompartmentInterval) intervals.get(0).deepClone();

            Set<SubcompartmentInterval> newIntervals = new HashSet<>();
            for (int i = 1; i < intervals.size(); i++) {
                SubcompartmentInterval nextInterval = intervals.get(i);
                if (collapsedInterval.overlapsWith(nextInterval)) {
                    collapsedInterval = collapsedInterval.absorbAndReturnNewInterval(nextInterval);
                } else {
                    newIntervals.add(collapsedInterval);
                    collapsedInterval = (SubcompartmentInterval) nextInterval.deepClone();
                }
            }
            newIntervals.add(collapsedInterval);

            List<SubcompartmentInterval> newIntervalsSorted = new ArrayList<>(newIntervals);
            Collections.sort(newIntervalsSorted);

            return newIntervalsSorted;
        }
        return intervals;
    }

    public static void splitGWList(GenomeWide1DList<SubcompartmentInterval> intraSubcompartments, int width) {
        intraSubcompartments.filterLists((chr, featureList) -> splitSubcompartmentIntervals(featureList, width));
    }

    private static List<SubcompartmentInterval> splitSubcompartmentIntervals(List<SubcompartmentInterval> intervals, int width) {
        if (intervals.size() > 0) {

            Collections.sort(intervals);

            Set<SubcompartmentInterval> newIntervals = new HashSet<>();
            for (SubcompartmentInterval currInterval : intervals) {
                newIntervals.addAll(currInterval.splitByWidth(width));
            }

            List<SubcompartmentInterval> newIntervalsSorted = new ArrayList<>(newIntervals);
            Collections.sort(newIntervalsSorted);

            return newIntervalsSorted;
        }
        return intervals;
    }


    public static String cleanUpPath(String filePath) {
        String[] breakUpFileName = filePath.split("/");
        return breakUpFileName[breakUpFileName.length - 1].replaceAll(".hic", "");
    }

    public static void readInFileAndCollapse(String location, ChromosomeHandler handler) {
        GenomeWide1DList<SubcompartmentInterval> subcompartments = loadFromSubcompartmentBEDFile(handler, location);
        SliceUtils.collapseGWList(subcompartments);
        subcompartments.simpleExport(new File(location + "_collapsed.bed"));
    }

    public static void readInFileAndSplitByResolutionLevel(String location, ChromosomeHandler handler) {

        GenomeWide1DList<SubcompartmentInterval> subcompartments = loadFromSubcompartmentBEDFile(handler, location);
        SliceUtils.splitGWList(subcompartments, 100000);
        subcompartments.simpleExport(new File(location + "_split.bed"));

    }

    /**
     * @return List of motif anchors from the provided bed file
     */
    public static GenomeWide1DList<SubcompartmentInterval> loadFromSubcompartmentBEDFile(ChromosomeHandler handler, String bedFilePath) {
        List<SubcompartmentInterval> anchors = new ArrayList<>();

        try {
            //BufferedReader br = ParsingUtils.openBufferedReader(bedFilePath);
            anchors.addAll(parseSubcompartmentBEDFile(bedFilePath, handler));
        } catch (IOException ec) {
            ec.printStackTrace();
        }

        return new GenomeWide1DList<>(handler, anchors);
    }

    /*
     * Methods for handling BED Files
     */

    /**
     * Helper function for actually parsing BED file
     * Ignores any attributes beyond the 5th column (i.e. just chr, positions, and subcompartment are read)
     *
     * @return list of motifs
     */
    private static List<SubcompartmentInterval> parseSubcompartmentBEDFile(String bedFilePath, ChromosomeHandler handler) throws IOException {
        BufferedReader bufferedReader = new BufferedReader(new FileReader(bedFilePath), MixerGlobals.bufferSize);

        Set<SubcompartmentInterval> anchors = new HashSet<>();
        String nextLine;
        Map<String, Integer> idToVal = new HashMap<>();

        int errorCount = 0;
        while ((nextLine = bufferedReader.readLine()) != null) {
            String[] tokens = Pattern.compile("\t").split(nextLine);

            if (tokens.length > 3 && tokens[3].equalsIgnoreCase("NA")) {
                continue;
            }

            if (tokens[0].startsWith("chr") && tokens.length > 3) {
                // valid line
                String chr1Name = tokens[0];
                int start1 = Integer.parseInt(tokens[1]);
                int end1 = Integer.parseInt(tokens[2]);
                String id = tokens[3].toUpperCase();

                int val;
                try {
                    val = Integer.parseInt(tokens[3]);
                } catch (Exception e) {
                    val = 0;
                }
                val = updateMapAndConfirmVal(id, val, idToVal);

                Chromosome chrom = handler.getChromosomeFromName(chr1Name);
                if (chrom == null) {
                    if (errorCount < 3) {
                        System.out.println("Skipping line: " + nextLine);
                    } else if (errorCount == 3) {
                        System.err.println("Maximum error count exceeded.  Further errors will not be logged");
                        System.err.println("If chromosomes were intentionally filtered, ignore this error.");
                    }
                    errorCount++;
                    continue;
                }

                anchors.add(new SubcompartmentInterval(chrom, start1, end1, val, id));
            } else {
                System.out.println("Skipping line: " + nextLine);
            }
        }
        if (anchors.size() < 1) System.err.println("BED File empty - file may have problems or error was encountered");
        bufferedReader.close();
        return new ArrayList<>(anchors);
    }

    private static int updateMapAndConfirmVal(String id, int val, Map<String, Integer> idToVal) {
        if (idToVal.containsKey(id)) {
            return idToVal.get(id);
        } else {
            int newVal = val;
            boolean valueIsAlreadyPresent = false;
            for (Integer valInMap : idToVal.values()) {
                valueIsAlreadyPresent |= valInMap.intValue() == val;
            }

            if (valueIsAlreadyPresent) {
                newVal = idCounter.incrementAndGet();
            } else {
                idCounter.set(val);
                idCounter.incrementAndGet();
            }
            idToVal.put(id, newVal);
            return newVal;
            //int val2 = idCounter.incrementAndGet();
            //idToVal.put(id, val2);
            //return val2;
        }
    }
}

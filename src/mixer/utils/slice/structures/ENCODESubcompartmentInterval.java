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

package mixer.utils.slice.structures;

import javastraw.feature1D.Feature1D;

public class ENCODESubcompartmentInterval extends SubcompartmentInterval {

    private final int[] otherIDs;
    private final int[] clusterSizes;

    public ENCODESubcompartmentInterval(int chrIndex, String chrName, int x1, int x2,
                                        Integer clusterID, int[] clusterSizes,
                                        int[] otherIDs) {
        super(chrIndex, chrName, x1, x2, clusterID);
        this.clusterSizes = clusterSizes;
        this.otherIDs = otherIDs;
    }

    @Override
    public Feature1D deepClone() {
        return new ENCODESubcompartmentInterval(getChrIndex(), getChrName(), getX1(), getX2(),
                getClusterID(), clusterSizes, otherIDs);
    }

    @Override
    public SubcompartmentInterval absorbAndReturnNewInterval(SubcompartmentInterval interval) {
        return new ENCODESubcompartmentInterval(getChrIndex(), getChrName(), getX1(), interval.getX2(),
                getClusterID(), clusterSizes, otherIDs);
    }

    @Override
    public boolean overlapsWith(SubcompartmentInterval o) {
        if (getChrIndex().equals(o.getChrIndex()) && getX2().equals(o.getX1())) {
            if (o instanceof ENCODESubcompartmentInterval) {
                ENCODESubcompartmentInterval o2 = (ENCODESubcompartmentInterval) o;
                return confirmMatchOfAllIDs(o2.otherIDs);
            }
        }
        return false;
    }

    private boolean confirmMatchOfAllIDs(int[] otherIDs2) {
        boolean everythingMatchesSoFar = true;
        for (int z = 0; z < otherIDs.length; z++) {
            everythingMatchesSoFar &= (otherIDs[z] == otherIDs2[z]);
        }
        return everythingMatchesSoFar;
    }

    @Override
    public String toString() {
        return "chr" + clean(getChrName()) + "\t" + getX1() + "\t" + getX2() + "\t" + (getClusterID() + 1) +
                "\t0\t.\t" + getX1() + "\t" + getX2() +
                "\t" + SubcompartmentColors.getColorString((getClusterID() + 1)) + listOutAllIDs();
    }

    private String clean(String chrName) {
        return chrName.replaceFirst("chr", "");
    }

    public static String getHeader() {
        return "#chrom\tchromStart\tchromEnd\tname\tscore\tstrand" +
                "\tthickStart\tthickEnd\titemRGB\tnumAltClusterings" +
                "\taltClusterNum\taltClusterAssignment";
    }

    private String listOutAllIDs() {
        StringBuilder allIDs = new StringBuilder("\t").append(otherIDs.length);
        allIDs.append("\t").append((clusterSizes[0] + 1));
        for (int i = 1; i < otherIDs.length; i++) {
            allIDs.append(",").append(clusterSizes[i]);
        }
        allIDs.append(",").append(getWidthForResolution(1));
        allIDs.append("\t").append(otherIDs[0]);
        for (int i = 1; i < otherIDs.length; i++) {
            allIDs.append(",").append(otherIDs[i]);
        }
        allIDs.append(",").append(0);
        return allIDs.toString();
    }
}

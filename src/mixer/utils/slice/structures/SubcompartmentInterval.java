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

import javastraw.feature1D.Feature1D;
import javastraw.reader.basics.Chromosome;

import java.util.ArrayList;
import java.util.List;

public class SubcompartmentInterval extends SimpleInterval {

    private Integer clusterID;
    private String clusterName;

    public SubcompartmentInterval(int chrIndex, String chrName, int x1, int x2, Integer clusterID, String clusterName) {
        super(chrIndex, chrName, x1, x2);
        this.clusterName = clusterName;
        this.clusterID = clusterID;
    }

    public SubcompartmentInterval(int chrIndex, String chrName, int x1, int x2, Integer clusterID) {
        this(chrIndex, chrName, x1, x2, clusterID, "" + clusterID);
    }

    public SubcompartmentInterval(Chromosome chromosome, int x1, int x2, Integer clusterID) {
        this(chromosome.getIndex(), chromosome.getName(), x1, x2, clusterID, "" + clusterID);
    }

    public SubcompartmentInterval(Chromosome chromosome, int x1, int x2, Integer clusterID, String clusterName) {
        this(chromosome.getIndex(), chromosome.getName(), x1, x2, clusterID, clusterName);
    }

    public Integer getClusterID() {
        return clusterID;
    }

    public void setClusterID(Integer clusterID) {
        this.clusterID = clusterID;
        this.clusterName = "" + clusterID;
    }

    public void setClusterID(Integer clusterID, String clusterName) {
        this.clusterID = clusterID;
        this.clusterName = clusterName;
    }

    public SubcompartmentInterval absorbAndReturnNewInterval(SubcompartmentInterval interval) {
        return new SubcompartmentInterval(getChrIndex(), getChrName(), getX1(), interval.getX2(), clusterID, clusterName);
    }

    public boolean overlapsWith(SubcompartmentInterval o) {
        return getChrIndex().equals(o.getChrIndex()) && clusterID.equals(o.clusterID) && getX2().equals(o.getX1());
    }

    @Override
    public String toString() {
        return (getChrName().startsWith("chr") ? getChrName() : "chr" + getChrName()) + "\t" + getX1() + "\t" + getX2() + "\t" + clusterName + "\t0"
                + "\t.\t" + getX1() + "\t" + getX2() + "\t" + SubcompartmentColors.getColorString(clusterID);
    }

    @Override
    public Feature1D deepClone() {
        return new SubcompartmentInterval(getChrIndex(), getChrName(), getX1(), getX2(), clusterID, clusterName);
    }

    public List<SubcompartmentInterval> splitByWidth(int width) {
        List<SubcompartmentInterval> splitList = new ArrayList<>();

        for (int i = getX1(); i < getX2(); i += width) {
            splitList.add(new SubcompartmentInterval(getChrIndex(), getChrName(), i, Math.min(i + width, getX2()), clusterID, clusterName));
        }

        return splitList;
    }
}
/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2011-2020 Rice University, Baylor College of Medicine, Aiden Lab
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


package mixer.track.feature;

import mixer.MixerGlobals;
import mixer.data.ChromosomeHandler;
import mixer.data.anchor.MotifAnchor;

import java.awt.*;
import java.util.List;
import java.util.*;


/**
 * @author jrobinso, mshamim, mhoeger
 *         <p/>
 *         reflection only used for plotting, should not be used by CLTs
 */
public class Feature2D implements Comparable<Feature2D> {

    static final String genericHeader = "#chr1\tx1\tx2\tchr2\ty1\ty2\tname\tscore\tstrand1\tstrand2\tcolor";
    private static final String genericLegacyHeader = "#chr1\tx1\tx2\tchr2\ty1\ty2\tcolor";
    private static final String BEDPE_SPACER = "\t.\t.\t.\t.";
    public static int tolerance = 0;
    final FeatureType featureType;
    final Map<String, String> attributes;
    final String chr1;
    final String chr2;
    final int start1;
    final int start2;
    int end1;
    int end2;
    private boolean isSelected = false;
    private Feature2D reflection = null;
    private Color color, translucentColor;

    public Feature2D(FeatureType featureType, String chr1, int start1, int end1, String chr2, int start2, int end2, Color c,
                     Map<String, String> attributes) {
        this.featureType = featureType;
        this.chr1 = chr1;
        this.start1 = start1;
        this.end1 = end1;
        this.chr2 = chr2;
        this.start2 = start2;
        this.end2 = end2;
        this.color = (c == null ? Color.black : c);
        setTranslucentColor();
        this.attributes = attributes;
    }

    public static String getDefaultOutputFileHeader() {
        if (MixerGlobals.isLegacyOutputPrintingEnabled) {
            return genericLegacyHeader;
        } else {
            return genericHeader;
        }
    }

    public FeatureType getFeatureType() {
        return this.featureType;
    }

    private String getFeatureName() {
        switch (featureType) {
            case PEAK:
                return "Peak";
            case DOMAIN:
                return "Contact Domain";
            case GENERIC:
            case NONE:
            default:
                return "Feature";
        }
    }

    public String getChr1() {
        return chr1;
    }

    public String getChr2() {
        return chr2;
    }

    public int getStart1() {
        return start1;
    }

    public int getStart2() {
        return start2;
    }

    public int getEnd1() {
        return end1;
    }

    public void setEnd1(int end1) {
        this.end1 = end1;
        if (reflection != null)
            reflection.end2 = end1;
    }

    public int getEnd2() {
        return end2;
    }

    public void setEnd2(int end2) {
        this.end2 = end2;
        if (reflection != null)
            reflection.end1 = end2;
    }

    public int getWidth1() {
        return end1 - start1;
    }

    public int getWidth2() {
        return end2 - start2;
    }

    public int getMidPt1() {
        return midPoint(start1, end1);
    }

    public int getMidPt2() {
        return midPoint(start2, end2);
    }

    private int midPoint(int start, int end) {
        return (int) (start + (end - start) / 2.0);
    }

    public Color getColor() {
        return color;
    }


    public void setColor(Color color) {
        this.color = color;
        if (reflection != null)
            reflection.color = color;
        setTranslucentColor();
    }

    private void setTranslucentColor() {
        translucentColor = new Color(color.getRed(), color.getGreen(), color.getBlue(), 50);
        if (reflection != null)
            reflection.translucentColor = translucentColor;
    }

    public String getOutputFileHeader() {
        StringBuilder output = new StringBuilder(getDefaultOutputFileHeader());

        ArrayList<String> keys = new ArrayList<>(attributes.keySet());
        Collections.sort(keys);

        for (String key : keys) {
            output.append("\t").append(key);
        }

        return output.toString();
    }

    private String simpleString() {
        return chr1 + "\t" + start1 + "\t" + end1 + "\t" + chr2 + "\t" + start2 + "\t" + end2;
    }

    private String justColorString() {
        return "\t" + color.getRed() + "," + color.getGreen() + "," + color.getBlue();
    }

    public String simpleStringWithColor() {
        if (MixerGlobals.isLegacyOutputPrintingEnabled) {
            return simpleString() + justColorString();
        } else {
            return simpleString() + BEDPE_SPACER + justColorString();
        }
    }

    @Override
    public String toString() {
        StringBuilder output = new StringBuilder(simpleStringWithColor());

        ArrayList<String> keys = new ArrayList<>(attributes.keySet());
        Collections.sort(keys);
        for (String key : keys) {
            output.append("\t").append(attributes.get(key));
        }

        return output.toString();
    }

    public ArrayList<String> getAttributeKeys() {
        ArrayList<String> keys = new ArrayList<>(attributes.keySet());
        Collections.sort(keys);
        return keys;
    }

    public String getAttribute(String key) {
        return attributes.get(key);
    }

    public void setAttribute(String key, String newVal) {
        attributes.put(key, newVal);
        // attribute directly shared between reflections
        if (reflection != null)
            reflection.attributes.put(key, newVal);
    }

    public float getFloatAttribute(String key) {
        return Float.parseFloat(attributes.get(key));
    }

    public void addIntAttribute(String key, int value) {
        attributes.put(key, "" + value);
    }

    public void addFloatAttribute(String key, Float value) {
        attributes.put(key, "" + value);
    }

    public void addStringAttribute(String key, String value) {
        attributes.put(key, value);
    }

    /**
     * @param otherFeature
     * @return
     */
    public boolean overlapsWith(Feature2D otherFeature) {

        float window1 = (otherFeature.getEnd1() - otherFeature.getStart1()) / 2;
        float window2 = (otherFeature.getEnd2() - otherFeature.getStart2()) / 2;

        int midOther1 = otherFeature.getMidPt1();
        int midOther2 = otherFeature.getMidPt2();

        return midOther1 >= (this.start1 - window1) && midOther1 <= (this.end1 + window1) && midOther2 >= (this.start2 - window2) && midOther2 <= (this.end2 + window2);
    }

    @Override
    public int compareTo(Feature2D o) {
        int[] comparisons = new int[]{chr1.compareTo(o.chr1), chr2.compareTo(o.chr2), start1 - o.start1,
                start2 - o.start2, end1 - o.end1, end2 - o.end2};
        for (int i : comparisons) {
            if (i != 0)
                return i;
        }
        return 0;
    }

    public boolean isOnDiagonal() {
        return chr1.equals(chr2) && start1 == start2 && end1 == end2;
    }

    public Feature2D reflectionAcrossDiagonal() {
        if (reflection == null) {
            reflection = new Feature2D(featureType, chr2, start2, end2, chr1, start1, end1, color, attributes);
            reflection.reflection = this;
        }
        return reflection;
    }

    public boolean isInLowerLeft() {
        return chr1.equals(chr2) && start2 > start1;
    }

    public boolean isInUpperRight() {
        return chr1.equals(chr2) && start2 < start1;
    }

    public boolean doesNotContainAttributeKey(String attribute) {
        return !attributes.containsKey(attribute);
    }

    public boolean containsAttributeValue(String attribute) {
        return attributes.containsValue(attribute);
    }

    public String getLocationKey() {
        return start1 + "_" + start2;
    }

    @Override
    public boolean equals(Object obj) {
        if (obj == null) {
            return false;
        }
        if (getClass() != obj.getClass()) {
            return false;
        }
        if (this == obj) {
            return true;
        }

        final Feature2D other = (Feature2D) obj;
        if (chr1.equals(other.chr1)) {
            if (chr2.equals(other.chr2)) {
                if (Math.abs(start1 - other.start1) <= tolerance) {
                    if (Math.abs(start2 - other.start2) <= tolerance) {
                        if (Math.abs(end1 - other.end1) <= tolerance) {
                            return Math.abs(end2 - other.end2) <= tolerance;
                        }
                    }
                }
            }
        }

        return false;
    }

    @Override
    public int hashCode() {
        int hash = 3;
        hash = 53 * hash + chr1.hashCode() + end1 - start1;
        hash = 53 * hash + chr2.hashCode() + end2 - start2;
        return hash;
    }

    public void clearAttributes() {
        attributes.clear();
    }

    public List<MotifAnchor> getAnchors(boolean onlyUninitializedFeatures, ChromosomeHandler handler) {
        List<Feature2D> originalFeatures = new ArrayList<>();
        originalFeatures.add(this);

        List<MotifAnchor> anchors = new ArrayList<>();
        if (isOnDiagonal()) {
            // loops should not be on diagonal
            // anchors.add(new MotifAnchor(chr1, start1, end1, originalFeatures, originalFeatures));
        } else {
            List<Feature2D> emptyList = new ArrayList<>();
            anchors.add(new MotifAnchor(chr1, start1, end1, originalFeatures, emptyList));
            anchors.add(new MotifAnchor(chr2, start2, end2, emptyList, originalFeatures));
        }
        return anchors;
    }

    public Feature2D deepCopy() {
        Map<String, String> attrClone = new HashMap<>();
        for (String key : attributes.keySet()) {
            attrClone.put(key, attributes.get(key));
        }
        return new Feature2D(featureType, chr1, start1, end1, chr2, start2, end2, color, attrClone);
    }

    public Feature2DWithMotif toFeature2DWithMotif() {
        return new Feature2DWithMotif(featureType, chr1, start1, end1, chr2, start2, end2, color, attributes);
    }

    public enum FeatureType {
        NONE, PEAK, DOMAIN, GENERIC, SCAFFOLD, SUPERSCAFFOLD, SELECTED_GROUP
    }
}

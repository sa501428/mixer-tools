/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2011-2017 Broad Institute, Aiden Lab
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

package juicebox.track;

import org.broad.igv.feature.Chromosome;
import org.broad.igv.renderer.DataRange;
import org.broad.igv.track.WindowFunction;

import java.awt.*;
import java.util.Collection;

/**
 * @author jrobinso
 *         Date: 8/1/13
 *         Time: 7:51 PM
 */
public interface HiCDataSource {

    String getName();

    void setName(String text);

    Color getPosColor();

    void setColor(Color selectedColor);

    Color getNegColor();

    void setNegColor(Color selectedColor);

    DataRange getDataRange();

    void setDataRange(DataRange dataRange);

    boolean isLog();

    HiCDataPoint[] getData(Chromosome chr, int startBin, int endBin, HiCGridAxis gridAxis, double scaleFactor, WindowFunction windowFunction);

    Collection<WindowFunction> getAvailableWindowFunctions();
}

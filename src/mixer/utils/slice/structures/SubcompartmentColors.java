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

import java.awt.*;

public class SubcompartmentColors {
    private static final int[][] COLORS = {
        {0,114,178},{230,159,0},{0,158,115},{213,62,0},{170,0,170},
        {0,190,210},{240,80,150},{110,70,30},{190,190,0},{50,40,120},
        {140,210,40},{35,35,35},{255,180,160},{100,0,45},{0,90,70},
        {180,140,220},{240,210,50},{100,160,230},{150,100,100},{0,45,65}
    };
    private static Color color(int id) {
        if (id <= 0) return new Color(128,128,128);
        if (id <= COLORS.length) {
            int[] c = COLORS[id-1]; return new Color(c[0],c[1],c[2]);
        }
        return new Color((int)(((long)id-20)*0x9e3779L) & 0xffffff);
    }
    public static Color getColorWithAlpha(int id, int alpha) {
        Color c = color(id);
        return new Color(c.getRed(),c.getGreen(),c.getBlue(),alpha);
    }
    public static String getColorString(Integer id) {
        Color c = color(id);
        return c.getRed()+","+c.getGreen()+","+c.getBlue();
    }
}

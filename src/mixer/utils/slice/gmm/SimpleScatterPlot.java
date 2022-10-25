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

package mixer.utils.slice.gmm;

import mixer.utils.slice.structures.SubcompartmentColors;
import mixer.utils.umap.UmapProjection;

import javax.imageio.ImageIO;
import java.awt.*;
import java.awt.image.BufferedImage;
import java.io.File;
import java.util.List;

public class SimpleScatterPlot {
    private static final Color BACKGROUND_COLOR = Color.BLACK; //Color.WHITE;
    private final int circleOffset = 3;
    private final int circleWidth = circleOffset * 2;
    private final float[][] points;
    private final int width = 1000;
    private final int height = 1000;
    private float minX = 0, maxX = 0;
    private float minY = 0, maxY = 0;
    private int widthX, heightY;
    private String stem = "";

    public SimpleScatterPlot(float[][] points) {
        this.points = UmapProjection.get2DProjection(points);
        if (points[0].length > 2) stem = "_umap";
        updateBounds();
    }

    private void updateBounds() {
        for (int k = 0; k < points.length; k++) {
            minX = Math.min(points[k][0], minX);
            maxX = Math.max(points[k][0], maxX);
            minY = Math.min(points[k][1], minY);
            maxY = Math.max(points[k][1], maxY);
        }
        minX -= circleWidth;
        minY -= circleWidth;
        maxX += circleWidth;
        maxY += circleWidth;
        widthX = (int) Math.ceil(maxX - minX);
        heightY = (int) Math.ceil(maxY - minY);
    }

    public void plot(List<List<Integer>> idList, String absPath) {
        int[] ids = new int[points.length];
        for (int i = 0; i < idList.size(); i++) {
            for (int k : idList.get(i)) {
                ids[k] = i;
            }
        }
        plot(ids, absPath);
    }

    public void plot(int[] ids, String absPath) {
        plotIndividualMap(ids, absPath + stem + ".png");
    }

    private void plotIndividualMap(int[] ids, String absPath) {
        BufferedImage image = new BufferedImage(width, height, BufferedImage.TYPE_INT_ARGB);
        Graphics2D g = (Graphics2D) image.getGraphics();
        colorBackground(g);

        for (int i = 0; i < points.length; i++) {
            drawCircle(g, points[i], ids[i]);
        }

        save(image, absPath);
    }

    private void colorBackground(Graphics2D g) {
        g.setColor(BACKGROUND_COLOR);
        g.fillRect(0, 0, width, height);
    }

    private void save(BufferedImage image, String absPath) {
        File file = new File(absPath);
        try {
            ImageIO.write(image, "png", file);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    private void drawCircle(Graphics g, float[] point, int id) {
        Color color = SubcompartmentColors.getColorWithAlpha(id + 1, 128);
        g.setColor(color);
        int newX = (int) ((point[0] - circleOffset - minX) * width / widthX);
        int newY = (int) ((point[1] - circleOffset - minY) * height / heightY);
        g.fillOval(newX, newY, circleWidth, circleWidth);
    }
}
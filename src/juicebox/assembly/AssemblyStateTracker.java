/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2011-2018 Broad Institute, Aiden Lab
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

package juicebox.assembly;

import juicebox.gui.SuperAdapter;

import java.util.Stack;

/**
 * Created by nathanielmusial on 7/5/17.
 */
public class AssemblyStateTracker {
    private final Stack<AssemblyScaffoldHandler> undoStack;
    private final Stack<AssemblyScaffoldHandler> redoStack;
    private final AssemblyScaffoldHandler initialAssemblyScaffoldHandler;
    private final SuperAdapter superAdapter;

    public AssemblyStateTracker(AssemblyScaffoldHandler assemblyScaffoldHandler, SuperAdapter superAdapter) {

        undoStack = new Stack<>();
        undoStack.push(assemblyScaffoldHandler);
        redoStack = new Stack<>();
        this.initialAssemblyScaffoldHandler = assemblyScaffoldHandler;
        this.superAdapter = superAdapter;
    }

    public AssemblyScaffoldHandler getAssemblyHandler() {
        return undoStack.peek();
    }

    public AssemblyScaffoldHandler getNewAssemblyHandler() {
        return new AssemblyScaffoldHandler(undoStack.peek());
    }

    public void resetState() {
        undoStack.clear();
        assemblyActionPerformed(initialAssemblyScaffoldHandler, true);
        superAdapter.safeClearAllMZDCache();
    }

    public AssemblyScaffoldHandler getInitialAssemblyScaffoldHandler() {
        return initialAssemblyScaffoldHandler;
    }

    public void assemblyActionPerformed(AssemblyScaffoldHandler assemblyScaffoldHandler, boolean refreshMap) {
        redoStack.clear();
        undoStack.push(assemblyScaffoldHandler);
        while (undoStack.size() > 50) { //keeps stack at size of 50 or less
            undoStack.remove(0);
        }
        assemblyScaffoldHandler.updateAssembly(refreshMap);
        regenerateLayers(refreshMap);
    }


    private void regenerateLayers(boolean refreshMap) {
        AssemblyScaffoldHandler assemblyScaffoldHandler = undoStack.peek();
        if (refreshMap) {
            superAdapter.getMainLayer().getFeatureHandler().setLoopList(assemblyScaffoldHandler.getScaffoldFeature2DHandler().getAllVisibleLoops());
        }
        superAdapter.getGroupLayer().getFeatureHandler().setLoopList(assemblyScaffoldHandler.getSuperscaffoldFeature2DHandler().getAllVisibleLoops());
    }

    public boolean checkUndo() {
        return undoStack.size() > 1;
    }

    public void undo() {
        if (checkUndo()) {
            redoStack.push(undoStack.pop());

            undoStack.peek().updateAssembly(true);
            regenerateLayers(true);
            superAdapter.safeClearAllMZDCache();
        }
    }

    public boolean checkRedo() {
        return !redoStack.empty();
    }

    public void redo() {
        if (checkRedo()) {
            undoStack.push(redoStack.pop());

            undoStack.peek().updateAssembly(true);
            regenerateLayers(true);
            superAdapter.safeClearAllMZDCache();
        }
    }
}
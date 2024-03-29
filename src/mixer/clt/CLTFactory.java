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

package mixer.clt;

import mixer.MixerGlobals;
import mixer.algos.*;


/**
 * Factory for command line tools to call different functions
 *
 * @author Muhammad Shamim
 * @since 1/30/2015
 */
public class CLTFactory {

    public static void generalUsage() {
        System.out.println("Mixer Tools Version " + MixerGlobals.versionNum);
        System.out.println("Usage:");
        System.out.println("\t" + "-h, --help print help");
        System.out.println("\t" + "-v, --verbose verbose mode");
        System.out.println("\t" + "-V, --version print version");
        System.out.println("Tool(s): slice");
        System.out.println("Type mixer_tools <commandName> for more detailed usage instructions");
    }

    public static MixerCLT getCLTCommand(String cmd) {

        cmd = cmd.toLowerCase();
        if (cmd.startsWith("slice") || cmd.startsWith("dice")) {
            return new Slice(cmd);
        } else if (cmd.equals("walk")) {
            return new ChromosomeWalker();
        } else if (cmd.contains("shuffle") || cmd.contains("umap")) {
            return new ShuffleAndUMAP(cmd);
        } else if (cmd.contains("rename")) {
            return new Rename(cmd);
        } else if (cmd.equals("network")) {
            return new Network();
        } else if (cmd.equals("aba")) {
            return new ABA();
        } else if (cmd.equals("finetune")) {
            return new FineTune();
        }

        return null;
    }
}

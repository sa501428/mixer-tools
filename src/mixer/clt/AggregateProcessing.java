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

package mixer.clt;

import mixer.MixerGlobals;
import mixer.MixerTools;


/**
 * Created for testing multiple CLTs at once
 * Basically scratch space
 */
@SuppressWarnings({"UnusedAssignment", "unused"})
public class AggregateProcessing {

    public static void main(String[] argv) throws Exception {

        String[] strings;

        String refs = "/Users/mshamim/Desktop/SLICE.Reboot/existing/GSE63525_GM12878_subcompartments.bed," +
                "/Users/mshamim/Desktop/SLICE.Reboot/existing/GM12878_SCI_sub_compartments.bed," +
                "/Users/mshamim/Desktop/SLICE.Reboot/existing/GM12878_track_hg19.bed";
        String file14 = "/Users/mshamim/Desktop/hicfiles/gm12878_rh14_30.hic";

        for (int res : new int[]{100000}) {
            String folder = "new_test_res_" + res + "v" + MixerGlobals.versionNum;
            strings = new String[]{"slice", "-r", res + "", "-k", "KR", //"--verbose",
                    file14, "2,11,3",
                    "/Users/mshamim/Desktop/SLICE.Reboot/" + folder, folder + "_"
            };
            System.out.println("-----------------------------------------------------");
            MixerTools.main(strings);
            System.gc();

            strings = new String[]{"shuffle", "-r", res + "", "-k", "KR", "-w", "" + (16 * (100000 / res)), file14,
                    //refs + "," +
                    "/Users/mshamim/Desktop/SLICE.Reboot/" + folder +
                            "/" + folder + "_5_clusters_gm12878_rh14_30.subcompartment.bed",
                    "/Users/mshamim/Desktop/SLICE.Reboot/shuffle_GM_res_" + res + "v" + MixerGlobals.versionNum,
                    //"RH2014,SCI,SNIPER," +
                    "SLICE3"
            };
            System.out.println("-----------------------------------------------------");
            MixerTools.main(strings);
            System.gc();
        }

        //ExtractionB4.extract();

        /*
        String[] strings;

        String refs = "/Users/mss/Desktop/SLICE.work/subcompartment_analysis/slice/existing/GSE63525_GM12878_subcompartments.bed," +
                "/Users/mss/Desktop/SLICE.work/subcompartment_analysis/slice/existing/GM12878_SCI_sub_compartments.bed," +
                "/Users/mss/Desktop/SLICE.work/subcompartment_analysis/slice/existing/GM12878_track_hg19.bed";
        String file14 = "/Users/mss/Desktop/hic_files/gm12878_rh14_30.hic";

        String folder = "new_change_test";
        strings = new String[]{"slice", "-r", "100000", "-k", "KR",
                "--type", "zscore-cosine", "--subsample", "50",
                file14, "2,11,10",
                "/Users/mss/Desktop/SLICE.work/tempslice/" + folder, folder + "_"
        };
        System.out.println("-----------------------------------------------------");
        //MixerTools.main(strings);
        System.gc();


        strings = new String[]{"shuffle", "-r", "100000", "-k", "KR", "-w", "16",
                "/Volumes/AidenLabWD7/Backup/AidenLab/LocalFiles/collins/map1_30.hic",
                //"/Users/mss/Desktop/SLICE.work/tempslice/collins_dice_test/collins_dice_test_4_clusters_map1_30.subcompartment.bed," +
                //       "/Users/mss/Desktop/SLICE.work/tempslice/collins_dice_test/collins_dice_test_5_clusters_map1_30.subcompartment.bed," +
                "/Users/mss/Desktop/SLICE.work/tempslice/collins_dice_test/collins_dice_test_10_clusters_map1_30.subcompartment.bed",
                "/Users/mss/Desktop/SLICE.work/tempslice/shuffle10",
                //"collins1_c4,collins1_c5," +
                "collins1_c10"
        };
        System.out.println("-----------------------------------------------------");
        //MixerTools.main(strings);
        System.gc();

        strings = new String[]{"shuffle", "-r", "100000", "-k", "KR", "-w", "16",
                "/Volumes/AidenLabWD7/Backup/AidenLab/LocalFiles/collins/map2_30.hic",
                "/Users/mss/Desktop/SLICE.work/tempslice/collins_dice_test/collins_dice_test_4_clusters_map1_30.subcompartment.bed," +
                        "/Users/mss/Desktop/SLICE.work/tempslice/collins_dice_test/collins_dice_test_5_clusters_map1_30.subcompartment.bed," +
                        "/Users/mss/Desktop/SLICE.work/tempslice/collins_dice_test/collins_dice_test_10_clusters_map1_30.subcompartment.bed",
                "/Users/mss/Desktop/SLICE.work/tempslice/shuffle10",
                "collins2_c4,collins2_c5,collins2_c10"
        };
        System.out.println("-----------------------------------------------------");
        //MixerTools.main(strings);
        System.gc();

        strings = new String[]{"shuffle", "-r", "100000", "-k", "KR", "-w", "16",
                file14,
                "/Users/mss/Desktop/SLICE.work/subcompartment_analysis/slice/existing/GSE63525_GM12878_subcompartments.bed," +
                        "/Users/mss/Desktop/SLICE.work/subcompartment_analysis/slice/existing/GM12878_track_hg19.bed," +
                        "/Users/mss/Desktop/SLICE.work/subcompartment_analysis/slice/existing/GM12878_SCI_sub_compartments.bed," +
                        "/Users/mss/Desktop/new_baseline_slice_gm2014/sliceNB_gm12878_rh14_100000_zscore-cosine2_clusters_gm12878_rh14_30.subcompartment.bed," +
                        "/Users/mss/Desktop/new_baseline_slice_gm2014/sliceNB_gm12878_rh14_100000_zscore-cosine5_clusters_gm12878_rh14_30.subcompartment.bed," +
                        "/Users/mss/Desktop/new_baseline_slice_gm2014/sliceNB_gm12878_rh14_100000_zscore-cosine6_clusters_gm12878_rh14_30.subcompartment.bed," +
                        "/Users/mss/Desktop/new_baseline_slice_gm2014/sliceNB_gm12878_rh14_100000_zscore-cosine4_clusters_gm12878_rh14_30.subcompartment.bed," +
                        "/Users/mss/Desktop/new_baseline_slice_gm2014/sliceNB_gm12878_rh14_100000_zscore-cosine10_clusters_gm12878_rh14_30.subcompartment.bed",
                "/Users/mss/Desktop/SLICE.work/tempslice/shuffle8",
                "rh2014,sniper,sci,sliceNB2,sliceNB5,sliceNB6,sliceNB4,sliceNB10"
        };
        System.out.println("-----------------------------------------------------");
        //MixerTools.main(strings);
        System.gc();

        strings = new String[]{"shuffle", "-r", "100000", "-k", "KR", "-w", "16",
                file14,
                "/Users/mss/Desktop/new_baseline_slice_gm2014/sliceNB_5_clusters_gm12878_rh14_reordered.bed," +
                        "/Users/mss/Desktop/SLICE.work/subcompartment_analysis/slice/existing/GM12878_track_hg19.bed," +
                        "/Users/mss/Desktop/SLICE.work/subcompartment_analysis/slice/existing/GM12878_SCI_sub_compartments.bed," +
                        "/Users/mss/Desktop/SLICE.work/subcompartment_analysis/slice/existing/GSE63525_GM12878_subcompartments.bed",
                "/Users/mss/Desktop/SLICE.work/tempslice/shuffle20",
                //"sliceNB_5clusters_reorder_1_skipsmall"
                "SLICE5,SNIPER,SCI,RH2014"
        };
        System.out.println("-----------------------------------------------------");
        //MixerTools.main(strings);
        System.gc();

        ///////////////////////////////////////////
        ///////////////////////////////////////////
        ///////////////////////////////////////////

        String stem = "/Users/mss/Desktop/prosen/";

        String dice50K = stem + "dice_pro_sen/50000/diceNB50000_pro_sen_zscore-cosine10_clusters_mss_pro_30.subcompartment.bed," +
                stem + "dice_pro_sen/50000/diceNB50000_pro_sen_zscore-cosine2_clusters_mss_pro_30.subcompartment.bed," +
                stem + "dice_pro_sen/50000/diceNB50000_pro_sen_zscore-cosine3_clusters_mss_pro_30.subcompartment.bed," +
                stem + "dice_pro_sen/50000/diceNB50000_pro_sen_zscore-cosine4_clusters_mss_pro_30.subcompartment.bed," +
                stem + "dice_pro_sen/50000/diceNB50000_pro_sen_zscore-cosine5_clusters_mss_pro_30.subcompartment.bed," +
                stem + "dice_pro_sen/50000/diceNB50000_pro_sen_zscore-cosine6_clusters_mss_pro_30.subcompartment.bed," +
                stem + "dice_pro_sen/50000/diceNB50000_pro_sen_zscore-cosine7_clusters_mss_pro_30.subcompartment.bed," +
                stem + "dice_pro_sen/50000/diceNB50000_pro_sen_zscore-cosine8_clusters_mss_pro_30.subcompartment.bed," +
                stem + "dice_pro_sen/50000/diceNB50000_pro_sen_zscore-cosine9_clusters_mss_pro_30.subcompartment.bed";

        String slicepro = stem + "slice_pro/50000/slc50000_pro_zscore-cosine10_clusters_mss_pro_30.subcompartment.bed," +
                stem + "slice_pro/50000/slc50000_pro_zscore-cosine2_clusters_mss_pro_30.subcompartment.bed," +
                stem + "slice_pro/50000/slc50000_pro_zscore-cosine3_clusters_mss_pro_30.subcompartment.bed," +
                stem + "slice_pro/50000/slc50000_pro_zscore-cosine4_clusters_mss_pro_30.subcompartment.bed," +
                stem + "slice_pro/50000/slc50000_pro_zscore-cosine5_clusters_mss_pro_30.subcompartment.bed," +
                stem + "slice_pro/50000/slc50000_pro_zscore-cosine6_clusters_mss_pro_30.subcompartment.bed," +
                stem + "slice_pro/50000/slc50000_pro_zscore-cosine7_clusters_mss_pro_30.subcompartment.bed," +
                stem + "slice_pro/50000/slc50000_pro_zscore-cosine8_clusters_mss_pro_30.subcompartment.bed," +
                stem + "slice_pro/50000/slc50000_pro_zscore-cosine9_clusters_mss_pro_30.subcompartment.bed";

        String slicesen = stem + "slice_sen/50000/slc50000_sen_zscore-cosine10_clusters_mss_sen_30.subcompartment.bed," +
                stem + "slice_sen/50000/slc50000_sen_zscore-cosine2_clusters_mss_sen_30.subcompartment.bed," +
                stem + "slice_sen/50000/slc50000_sen_zscore-cosine3_clusters_mss_sen_30.subcompartment.bed," +
                stem + "slice_sen/50000/slc50000_sen_zscore-cosine4_clusters_mss_sen_30.subcompartment.bed," +
                stem + "slice_sen/50000/slc50000_sen_zscore-cosine5_clusters_mss_sen_30.subcompartment.bed," +
                stem + "slice_sen/50000/slc50000_sen_zscore-cosine6_clusters_mss_sen_30.subcompartment.bed," +
                stem + "slice_sen/50000/slc50000_sen_zscore-cosine7_clusters_mss_sen_30.subcompartment.bed," +
                stem + "slice_sen/50000/slc50000_sen_zscore-cosine8_clusters_mss_sen_30.subcompartment.bed," +
                stem + "slice_sen/50000/slc50000_sen_zscore-cosine9_clusters_mss_sen_30.subcompartment.bed";

        String filePro = "/Users/mss/Desktop/hic_files/mss_pro_30.hic";
        String fileSen = "/Users/mss/Desktop/hic_files/mss_sen_30.hic";

        strings = new String[]{"shuffle", "-r", "50000", "-k", "KR", "-w", "30",
                filePro,
                dice50K,
                "/Users/mss/Desktop/prosen/shuffle/dice_Pro",
                "C10,C2,C3,C4,C5,C6,C7,C8,C9"
        };
        System.out.println("-----------------------------------------------------");
        //MixerTools.main(strings);
        System.gc();

        strings = new String[]{"shuffle", "-r", "50000", "-k", "KR", "-w", "30",
                fileSen,
                dice50K,
                "/Users/mss/Desktop/prosen/shuffle/dice_Sen",
                "C10,C2,C3,C4,C5,C6,C7,C8,C9"
        };
        System.out.println("-----------------------------------------------------");
        //MixerTools.main(strings);
        System.gc();

        strings = new String[]{"shuffle", "-r", "50000", "-k", "KR", "-w", "30",
                filePro,
                slicepro,
                "/Users/mss/Desktop/prosen/shuffle/slice_Pro",
                "C10,C2,C3,C4,C5,C6,C7,C8,C9"
        };
        System.out.println("-----------------------------------------------------");
        //MixerTools.main(strings);
        System.gc();


        strings = new String[]{"shuffle", "-r", "50000", "-k", "KR", "-w", "30",
                fileSen,
                slicesen,
                "/Users/mss/Desktop/prosen/shuffle/slice_Sen",
                "C10,C2,C3,C4,C5,C6,C7,C8,C9"
        };
        System.out.println("-----------------------------------------------------");
        //MixerTools.main(strings);
        System.gc();


        ///////////////////////////////////////////
        ///////////////////////////////////////////
        ///////////////////////////////////////////

        folder = "covid_dice";
        strings = new String[]{"dice", "-r", "100000", "-k", "KR",
                "/Users/mss/Desktop/hic_files/marianna/HIC68_30.hic," +
                        "/Users/mss/Desktop/hic_files/marianna/HIC80_30.hic",
                "2,7,10", // 2,11,10
                "/Users/mss/Desktop/hic_files/marianna/new_dice_V6_100K", folder + "_"
        };
        System.out.println("-----------------------------------------------------");
        MixerTools.main(strings);
        System.gc();


        stem = "/Users/mss/Desktop/hic_files/marianna/new_dice_V6_100K/";
        String dice100K = stem + "covid_dice_10_clusters_HIC68_30.subcompartment.bed," +
                stem + "covid_dice_2_clusters_HIC68_30.subcompartment.bed," +
                stem + "covid_dice_3_clusters_HIC68_30.subcompartment.bed," +
                stem + "covid_dice_4_clusters_HIC68_30.subcompartment.bed," +
                stem + "covid_dice_5_clusters_HIC68_30.subcompartment.bed," +
                stem + "covid_dice_6_clusters_HIC68_30.subcompartment.bed," +
                stem + "covid_dice_7_clusters_HIC68_30.subcompartment.bed," +
                stem + "covid_dice_8_clusters_HIC68_30.subcompartment.bed," +
                stem + "covid_dice_9_clusters_HIC68_30.subcompartment.bed";

        strings = new String[]{"shuffle", "-r", "100000", "-k", "KR", "-w", "30",
                "/Users/mss/Desktop/hic_files/marianna/HIC68_30.hic",
                dice100K,
                "/Users/mss/Desktop/hic_files/marianna/new_dice_V6_100K/shuffle_dice_100K_68.hic",
                "C10,C2,C3,C4,C5,C6,C7,C8,C9"
        };
        System.out.println("-----------------------------------------------------");
        MixerTools.main(strings);
        System.gc();

        strings = new String[]{"shuffle", "-r", "100000", "-k", "KR", "-w", "30",
                "/Users/mss/Desktop/hic_files/marianna/HIC80_30.hic",
                dice100K,
                "/Users/mss/Desktop/hic_files/marianna/new_dice_V6_100K/shuffle_dice_100K_80.hic",
                "C10,C2,C3,C4,C5,C6,C7,C8,C9"
        };
        System.out.println("-----------------------------------------------------");
        MixerTools.main(strings);
        System.gc();

         */
    }
}

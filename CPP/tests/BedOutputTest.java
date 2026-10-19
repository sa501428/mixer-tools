import mixer.utils.slice.structures.*;
import mixer.utils.slice.cleaning.CisCoverage;
import java.util.*;

public class BedOutputTest {
    private static void check(boolean ok) {
        if (!ok) throw new AssertionError("BED output contract failed");
    }
    public static void main(String[] args) {
        double[] counts=new double[100];Arrays.fill(counts,100);counts[0]=1;counts[1]=0;
        CisCoverage coverage=new CisCoverage(counts);
        check(coverage.bad[0]&&coverage.bad[1]&&coverage.z[0]<-2.5);
        check(Math.abs(coverage.mean-99)<1e-10&&Math.abs(coverage.sd-Math.sqrt(98))<1e-10);
        for(int i=2;i<counts.length;i++)check(!coverage.bad[i]);
        CisCoverage flat=new CisCoverage(new double[]{10,0,10});
        check(flat.sd==0&&!flat.bad[0]&&flat.bad[1]&&!flat.bad[2]);
        Set<String> colors = new HashSet<>();
        for (int id=0; id<=128; id++) check(colors.add(SubcompartmentColors.getColorString(id)));
        check(SubcompartmentColors.getColorString(1).equals("0,114,178"));
        List<SubcompartmentInterval> input = Arrays.asList(
            new SubcompartmentInterval(1,"chr1",10,20,1),
            new SubcompartmentInterval(1,"chr1",0,10,1),
            new SubcompartmentInterval(1,"chr1",30,40,1),
            new SubcompartmentInterval(1,"chr1",40,53,2),
            new SubcompartmentInterval(2,"chr2",0,10,2));
        List<SubcompartmentInterval> merged = SliceUtils.collapseSubcompartmentIntervals(input);
        check(merged.size()==4);
        check(merged.get(0).getX1()==0 && merged.get(0).getX2()==20);
        for (SubcompartmentInterval interval:merged) {
            String[] row=interval.toString().split("\t");
            check(row.length==9 && row[4].equals("0") && row[5].equals("."));
            check(!row[0].startsWith("chrchr"));
            check(row[1].equals(row[6]) && row[2].equals(row[7]));
        }
        List<SubcompartmentInterval> split=merged.get(2).splitByWidth(10);
        check(split.get(split.size()-1).getX2()==53);
        System.out.println("Java BED9 merging, gap preservation, clipping and palette passed");
    }
}

package mixer.utils.slice.cleaning;

/** Raw off-diagonal cis row nonzeros, evaluated independently per chromosome. */
public final class CisCoverage {
    public static final double CUTOFF = -2.5;
    public final double[] nonzero, z;
    public final boolean[] bad;
    public final double mean, sd;
    public CisCoverage(double[] counts) {
        nonzero=counts.clone(); z=new double[counts.length];bad=new boolean[counts.length];
        int n=0;double mu=0,m2=0;
        for(double value:counts)if(value>0){n++;double delta=value-mu;mu+=delta/n;m2+=delta*(value-mu);}
        mean=mu;sd=n>0?Math.sqrt(Math.max(0,m2/n)):0;
        for(int i=0;i<counts.length;i++){
            z[i]=sd>0?(counts[i]-mean)/sd:0;
            bad[i]=counts[i]<=0||z[i]<CUTOFF;
        }
    }
}

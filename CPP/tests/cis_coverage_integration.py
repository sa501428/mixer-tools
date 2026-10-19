"""Sparse cis rows must stay excluded despite valid norms and dense trans."""
import csv
from pathlib import Path
import subprocess
import sys
import tempfile
from hic_fixture import fixture


def run(binary,*args):
    result=subprocess.run([str(binary),*map(str,args)],capture_output=True,text=True,timeout=120)
    assert result.returncode==0,(args,result.stdout,result.stderr)


def check(out,patterns):
    report=next(out.glob('cis_coverage*.tsv'))
    rows=list(csv.DictReader(report.open(),delimiter='\t'))
    sparse=[r for r in rows if int(r['start'])//100000 in (55,66)]
    assert len(sparse)>=6
    for r in sparse:
        assert r['retained']=='0',r
        if int(r['start'])//100000==55:
            assert float(r['cis_nonzero'])==1 and float(r['zscore']) < -2.5,r
        else:assert float(r['cis_nonzero'])==0,r
    for pattern in patterns:
        files=list(out.glob(pattern));assert files,(out,pattern)
        for path in files:
            positive=0
            for line in path.read_text().splitlines():
                if not line or line.startswith('#'):continue
                r=line.split('\t')
                if r[3]=='0':continue
                positive+=1
                for bin in (55,66):
                    assert int(r[2])<=bin*100000 or int(r[1])>=(bin+1)*100000,(path,r,bin)
            assert positive,path
    return rows


def main():
    binary=Path(sys.argv[1]).resolve()
    with tempfile.TemporaryDirectory() as tmp:
        root=Path(tmp);hic=root/'cis_sparse.hic';fixture(hic,cis_sparse=True)
        for adaptive in (False,True):
            out=root/str(adaptive)
            extra=['--cluster-policy','adaptive','--normalization-views','dual','--k-range','2,3'] if adaptive else []
            run(binary,'--no-umap','--random-seeds',42,'-z',2,*extra,hic,'2,3,1',out,'test')
            check(out,['inter.bed','cis.bed','combined.bed'] if adaptive else ['test_2_kmedians_clusters.bed'])
    print('Sparse and empty raw cis rows excluded from final calls and pooling')

if __name__=='__main__':main()

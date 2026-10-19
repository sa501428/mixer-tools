#!/usr/bin/env python3
"""End-to-end .hic tests; no downloads, NumPy, Java, or preexisting data needed."""
import csv
import ast
import math
from pathlib import Path
import struct
import subprocess
import sys
from bed_contract import validate_bed
import tempfile
from hic_fixture import fixture


def npy(path):
    data = path.read_bytes()
    assert data[:8] == b'\x93NUMPY\x01\x00'
    size = struct.unpack_from('<H', data, 8)[0]
    header = ast.literal_eval(data[10:10+size].decode())
    fmt = {'<f8': 'd', '<f4': 'f', '<i4': 'i'}[header['descr']]
    return header, struct.unpack('<'+str(math.prod(header['shape']))+fmt, data[10+size:])


def run(binary, *args, ok=True):
    p = subprocess.run([str(binary), *map(str, args)], capture_output=True, text=True, timeout=90)
    assert (p.returncode == 0) == ok, (args, p.stdout, p.stderr)
    return p


def calls(path):
    rows = [line.split('\t') for line in path.read_text().splitlines() if not line.startswith('#')]
    assert rows and all(len(row) >= 9 for row in rows)
    for row in rows:
        assert 0 <= int(row[1]) < int(row[2]) <= 12000000
    return rows


def main():
    binary = Path(sys.argv[1]).resolve()
    with tempfile.TemporaryDirectory() as temp:
        root = Path(temp)
        hic = root/'data.hic'
        fixture(hic)
        # Short invocation selects K automatically and exports every screened K.
        for criterion in ('bic','aic'):
            out=root/('auto-'+criterion)
            run(binary,'--k-selection',criterion,'--k-range','1,4','--normalization-views','dual',hic,out,'auto')
            for view in ('inter','cis'):
                scores=list(csv.DictReader((out/(view+'_model_selection.tsv')).open(),delimiter='\t'))
                assert {int(r['k']) for r in scores}=={1,2,3,4}
                fit=list(csv.DictReader((out/'adaptive_fits.tsv').open(),delimiter='\t'))[view=='cis']
                assert fit['criterion']==criterion
                assert abs(float(fit['selection_score'])-min(float(r['score']) for r in scores))<.05
            validate_bed(out/'combined.bed')
        run(binary,'--k-selection','invalid',hic,root/'bad-auto',ok=False)
        for threads in (1, 3):
            out = root/str(threads)
            run(binary, 'slice', '--no-umap', '--both', '--encode-mode', '--random-seeds', '17,42',
                '-z', threads, hic, '2,6,2', out, 'test')
            for method in ('kmeans', 'kmedians'):
                for k in range(2, 6):
                    rows = validate_bed(out/f'test_{k}_{method}_clusters.bed')
                    assert {row[0] for row in rows} == {'chr1', 'chr2', 'chr3'}
                    assert {int(row[3]) for row in rows} == set(range(1, k+1))
                    assert all(int(rows[i][1]) >= int(rows[i-1][2]) for i in range(1, len(rows)) if rows[i][0] == rows[i-1][0])
                header, values = npy(out/f'{method}_cluster_post_stats.npy')
                assert header['shape'] == (6, 4) and values[:4] == (2, 3, 4, 5)
                assert all(math.isfinite(v) for v in values)
            header, values = npy(out/'slice_matrix.npy')
            assert header['descr'] == '<f4' and header['shape'][0] == 360
            assert any(math.isnan(v) for v in values) and any(math.isfinite(v) for v in values)
            assert npy(out/'genome_indices.npy')[0]['descr'] == '<i4'
            validate_bed(out/'slice_subcompartment_clusters.bed')
            assert all(len(row) == 12 for row in calls(out/'slice_subcompartment_clusterings.tsv'))
        for path in (root/'1').iterdir():
            if path.name=='analysis_parameters.tsv':
                left=[r for r in path.read_text().splitlines() if not r.startswith('threads\t')]
                right=[r for r in (root/'3'/path.name).read_text().splitlines() if not r.startswith('threads\t')]
                assert left==right
                continue
            assert path.read_bytes() == (root/'3'/path.name).read_bytes(), path.name
        # Projection is diagnostic: calls and saved matrix must be unchanged.
        run(binary, '--umap', '--random-seeds', 42, hic, '2,3,2', root/'projection', 'p')
        run(binary, '--no-umap', '--random-seeds', 42, hic, '2,3,2', root/'plain', 'p')
        for name in ('slice_matrix.npy', 'p_2_kmedians_clusters.bed', 'kmedians_cluster_post_stats.npy'):
            assert (root/'projection'/name).read_bytes() == (root/'plain'/name).read_bytes()
        assert npy(root/'projection/umap_projection.npy')[0]['shape'] == (360, 2)
        assert '<svg' in (root/'projection/umap_chrom.svg').read_text()
        # Adaptive paired fits use one compression, preserve both outputs and are deterministic.
        for threads in (1, 3):
            out=root/f'dual-{threads}'
            run(binary,'--cluster-policy','adaptive','--normalization-views','dual',
                '--model','gmm','--k-range','2,3','-z',threads,'--random-seeds',42,
                hic,'2,3,2',out,'dual')
            for name in ('inter.bed','cis.bed','combined.bed','state_lineage.tsv','view_footprint.tsv'):
                assert (out/name).stat().st_size>0,name
        for name in ('inter.bed','cis.bed','combined.bed','adaptive_summary.tsv','state_lineage.tsv'):
            assert (root/'dual-1'/name).read_bytes()==(root/'dual-3'/name).read_bytes(),name
        assert not (root/'dual-1/umap_projection.npy').exists()
        for name in ('inter.bed','cis.bed','combined.bed'):
            validate_bed(root/'dual-1'/name)
        def footprint(path):
            return {(r[0],x) for r in calls(path) for x in range(int(r[1]),int(r[2]),100000)}
        assert footprint(root/'dual-1/combined.bed') == footprint(root/'dual-1/inter.bed') | footprint(root/'dual-1/cis.bed')

        # A missing companion normalization preserves the inter result.
        run(binary,'--cis-normalization','MISSING','--cluster-policy','adaptive',
            '--normalization-views','dual','--k-range','2,3',hic,'2,3,1',root/'missing-cis','m')
        assert (root/'missing-cis/combined.bed').stat().st_size>0
        run(binary,'--cis-normalization','SCALE','--trans-normalization','SCALE','--post-balance',
            '--cluster-policy','adaptive','--k-range','2,3',hic,'2,3,1',root/'post-balanced','p')
        assert (root/'post-balanced/combined.bed').stat().st_size>0
        kr = root/'kr.hic'
        fixture(kr, norms=('KR', 'INTER_KR', 'GW_KR'), shift=3)
        run(binary, '--no-umap', f'{hic},{kr}', '2,3,1', root/'multi', 'm')
        assert npy(root/'multi/slice_matrix.npy')[0]['shape'][1] > npy(root/'plain/slice_matrix.npy')[0]['shape'][1]
        assert (root/'multi/initial_split_1.bed').is_file()
        run(binary, '--no-umap', '-c', 'chr1,2', hic, '2,3,1', root/'subset', 's')
        assert {row[0] for row in calls(root/'subset/s_2_kmedians_clusters.bed')} == {'chr1','chr2'}
        run(binary, '--no-umap', '--no-filter-outliers', hic, '2,3,1', root/'unfiltered', 'u')
        exact = root/'exact.hic'
        fixture(exact, exact_length=True)
        run(binary, '--no-umap', exact, '2,3,1', root/'exact', 'e')
        calls(root/'exact/e_2_kmedians_clusters.bed')
        scale = root/'scale.hic'
        fixture(scale, norms=('SCALE',))
        run(binary, '--no-umap', scale, '2,3,1', root/'missingnorm', 's')
        params=dict(line.split('\t') for line in (root/'missingnorm/analysis_parameters.tsv').read_text().splitlines()[1:])
        assert params['trans_source_normalization_0']=='NONE' and float(params['balance_error_0'])<=1e-4
        run(binary, '--java-compatible', '--no-umap', scale, '2,3,1', root/'missingnorm-java', 's', ok=False)
        run(binary, '--legacy-scale', '--no-umap', scale, '2,3,1', root/'legacy-scale', 's')
        run(binary, '--cis-normalization', 'SCALE', '--row-standardization', 'weighted-zscore', '--no-umap', scale, '2,3,1', root/'rowz', 's')
        run(binary, '--max-memory-gb', '.000001', '--no-umap', hic, '2,3,1', root/'memory', 's', ok=False)
        run(binary, '--no-umap', '-k', 'SCALE', scale, '2,3,1', root/'override', 's')
        assert (root/'legacy-scale/slice_matrix.npy').read_bytes()==(root/'override/slice_matrix.npy').read_bytes()
        vc=root/'inter-vc.hic';fixture(vc,norms=('SCALE','INTER_VC'))
        run(binary,'--no-umap','--trans-normalization','AUTO',vc,'2,3,1',root/'vc','s')
        params=dict(line.split('\t') for line in (root/'vc/analysis_parameters.tsv').read_text().splitlines()[1:])
        assert params['trans_source_normalization_0']=='INTER_VC'
        run(binary,'--no-umap','--trans-normalization','AUTO',scale,'2,3,1',root/'auto-raw','s')
        params=dict(line.split('\t') for line in (root/'auto-raw/analysis_parameters.tsv').read_text().splitlines()[1:])
        assert params['trans_source_normalization_0']=='NONE' and 'balance_error_0' in params
        rearranged = root/'rearranged.hic'
        fixture(rearranged, trans_value=lambda a,b,x,y: 1e25 if a==1 and b==2 and 5<=x<15 and 5<=y<15 else 10+(x+y)%8)
        run(binary, '--no-umap', '--has-translocation', rearranged, '2,3,1', root/'translocations', 't')
        detections = list((root/'translocations').glob('suspected_translocations_*.bedpe'))
        assert detections and any(p.stat().st_size for p in detections)
        for args in ([], ['--help'], ['--version']):
            run(binary, *args)
        for args in (['--bogus'], [hic, '2,2,1', root, 'x'], [hic, '2,3,0', root, 'x'],
                     ['-r', 30000, hic, '2,3,1', root, 'x'], ['-r', 50000, hic, '2,3,1', root, 'x'],
                     ['-c', 'missing', hic, '2,3,1', root, 'x'], ['-c', '1', hic, '2,3,1', root, 'x'],
                     [root/'absent.hic', '2,3,1', root, 'x'], [hic,'2,3,1',root,'../bad']):
            run(binary, *args, ok=False)
        missing = root/'missing.hic'
        fixture(missing, missing_pair=True)
        run(binary, '--no-umap', missing, '2,3,1', root/'missing', 'x', ok=False)
    print('Native .hic pipeline, both methods, ENCODE, projection, multi-input, normalization fallback, translocations, deterministic threads, CLI validation passed')


if __name__ == '__main__':
    main()

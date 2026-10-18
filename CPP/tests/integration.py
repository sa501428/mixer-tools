#!/usr/bin/env python3
"""End-to-end .hic tests; no downloads, NumPy, Java, or preexisting data needed."""
import ast
import math
from pathlib import Path
import struct
import subprocess
import sys
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
        for threads in (1, 3):
            out = root/str(threads)
            run(binary, 'slice', '--no-umap', '--both', '--encode-mode', '--random-seeds', '17,42',
                '-z', threads, hic, '2,6,2', out, 'test')
            for method in ('kmeans', 'kmedians'):
                for k in range(2, 6):
                    rows = calls(out/f'test_{k}_{method}_clusters.bed')
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
            assert all(len(row) == 12 for row in calls(out/'slice_subcompartment_clusters.bed'))
        for path in (root/'1').iterdir():
            assert path.read_bytes() == (root/'3'/path.name).read_bytes(), path.name
        # Projection is diagnostic: calls and saved matrix must be unchanged.
        run(binary, '--random-seeds', 42, hic, '2,3,2', root/'projection', 'p')
        run(binary, '--no-umap', '--random-seeds', 42, hic, '2,3,2', root/'plain', 'p')
        for name in ('slice_matrix.npy', 'p_2_kmedians_clusters.bed', 'kmedians_cluster_post_stats.npy'):
            assert (root/'projection'/name).read_bytes() == (root/'plain'/name).read_bytes()
        assert npy(root/'projection/umap_projection.npy')[0]['shape'] == (360, 2)
        assert '<svg' in (root/'projection/umap_chrom.svg').read_text()
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
        run(binary, '--no-umap', scale, '2,3,1', root/'missingnorm', 's', ok=False)
        run(binary, '--no-umap', '-k', 'SCALE', scale, '2,3,1', root/'override', 's')
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

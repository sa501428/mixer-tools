#!/usr/bin/env python3
"""Differential tests against the adjacent, checked-in Java source and JARs.

The reference source is copied into a temporary directory. Only diagnostic UMAP
is disabled; the second build exposes the existing USE_KMEANS switch. No Java
production source is modified, and Java is never used by the native executable.
"""
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
from hic_fixture import fixture
from integration import npy, run


def execute(args):
    p = subprocess.run(list(map(str, args)), capture_output=True, text=True, timeout=120)
    assert p.returncode == 0, (args, p.stdout, p.stderr)
    return p


def compare_array(a, b):
    ha, va = npy(a)
    hb, vb = npy(b)
    assert ha == hb, (a.name, ha, hb)
    for i, (x, y) in enumerate(zip(va, vb)):
        if math.isnan(x):
            assert math.isnan(y), (a.name, i, x, y)
        else:
            assert math.isclose(x, y, rel_tol=2e-6, abs_tol=2e-6), (a.name, i, x, y)


def assignments(path):
    """Ignore Java collapseGWList's duplicate first interval; compare bin calls."""
    values = {}
    for line in path.read_text().splitlines():
        if line.startswith('#'):
            continue
        row = line.split('\t')
        for start in range(int(row[1]), int(row[2]), 100000):
            key = (row[0], start)
            value = tuple(row[3:6]+row[8:])
            if key in values:
                assert values[key] == value
            values[key] = value
    return values


def compile_reference(reference, work, both=False, project=False):
    source = work/'src'
    shutil.copytree(reference/'src', source)
    entry = source/'mixer/algos/Slice.java'
    text = entry.read_text()
    if not project:
        text = text.replace('PROJECT_TO_UMAP = true', 'PROJECT_TO_UMAP = false')
    if both:
        text = text.replace('USE_KMEANS = false', 'USE_KMEANS = true')
    entry.write_text(text)
    classes = work/'classes'
    classpath = os.pathsep.join(str(reference/'lib'/part/'*') for part in ('general', 'umap'))
    execute(['javac', '-proc:none', '-cp', classpath, '-sourcepath', source, '-d', classes, source/'mixer/MixerTools.java'])
    return ['java', '-XX:ActiveProcessorCount=2', '-Djava.awt.headless=true', '-cp', str(classes)+os.pathsep+classpath, 'mixer.MixerTools', 'slice']


def main():
    binary = Path(sys.argv[1]).resolve()
    reference = Path(__file__).resolve().parents[2]
    with tempfile.TemporaryDirectory() as temp:
        root = Path(temp)
        hic = root/'a.hic'
        second = root/'b.hic'
        fixture(hic)
        fixture(second, norms=('KR', 'INTER_KR', 'GW_KR'), shift=3)
        java = compile_reference(reference, root/'default')
        for name, files, cluster_range in [('single', str(hic), '2,4,2'), ('multiple', f'{hic},{second}', '2,4,2')]:
            jout, cout = root/(name+'_java'), root/(name+'_cpp')
            args = ['--random-seeds', '17,42', '-z', '2', files, cluster_range]
            execute(java+args+[jout, 'test'])
            run(binary, '--java-compatible', '--no-umap', *args, cout, 'test')
            for file in ('slice_matrix.npy', 'genome_indices.npy', 'kmedians_cluster_post_stats.npy'):
                compare_array(jout/file, cout/file)
            for k in (2, 3):
                file = f'test_{k}_kmedians_clusters.bed'
                assert assignments(jout/file) == assignments(cout/file), file
            # The Java exporter overwrites initial_split.bed for each input;
            # the native port preserves each input's diagnostics separately.
            cpp_split = 'initial_split_1.bed' if name == 'multiple' else 'initial_split.bed'
            assert assignments(jout/'initial_split.bed') == assignments(cout/cpp_split)
        java_both = compile_reference(reference, root/'both', both=True)
        jout, cout = root/'both_java', root/'both_cpp'
        args = ['--random-seeds', '42', '--encode-mode', str(hic), '2,6,2']
        execute(java_both+args+[jout, 'test'])
        run(binary, '--java-compatible', '--both', '--no-umap', *args, cout, 'test')
        for method in ('kmeans', 'kmedians'):
            compare_array(jout/f'{method}_cluster_post_stats.npy', cout/f'{method}_cluster_post_stats.npy')
            for k in range(2, 6):
                file = f'test_{k}_{method}_clusters.bed'
                assert assignments(jout/file) == assignments(cout/file), file
        file = 'slice_subcompartment_clusters.bed'
        # ENCODE's alternative-clustering fields include interval widths, so
        # compare complete row sets after removing Java's duplicate first row.
        rows = lambda path: {tuple(line.split('\t')) for line in path.read_text().splitlines() if not line.startswith('#')}
        assert rows(jout/file) == rows(cout/file)
        java_projection = compile_reference(reference, root/'projection_reference', project=True)
        jout, cout = root/'projection_java', root/'projection_cpp'
        args = ['--random-seeds', '42', str(hic), '2,3,1']
        execute(java_projection+args+[jout, 'test'])
        run(binary, '--java-compatible', *args, cout, 'test')
        compare_array(jout/'corr_slice_matrix.npy', cout/'corr_slice_matrix.npy')
    print('Java differential parity passed: compression, cleaned matrix, coordinates, multi-input, both clustering methods, WCSS/AIC/BIC/correlation, ENCODE')


if __name__ == '__main__':
    main()

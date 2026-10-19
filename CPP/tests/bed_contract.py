"""Validate browser BED9 and expand merged records for per-bin assertions."""
def validate_bed(path, step=None):
    rows = [line.split('\t') for line in path.read_text().splitlines() if line and not line.startswith('#')]
    previous = None
    colors = {}
    for row in rows:
        assert len(row) == 9, (path, row)
        chrom, start, end, label, score, strand, thick_start, thick_end, rgb = row
        assert 0 <= int(start) < int(end), (path, row)
        assert score == '0' and strand == '.' and [thick_start, thick_end] == [start, end], (path, row)
        channels = list(map(int, rgb.split(',')))
        assert len(channels) == 3 and all(0 <= x <= 255 for x in channels), (path, row)
        assert colors.setdefault(label, rgb) == rgb, (path, label)
        if previous and previous[0] == chrom:
            assert int(previous[2]) <= int(start), (path, previous, row)
            assert previous[2] != start or previous[3] != label, (path, previous, row)
        previous = row
    assert len(set(colors.values())) == len(colors), (path, colors)
    if step is None:
        return rows
    return [[r[0], str(x), str(min(x+step, int(r[2]))), r[3]]
            for r in rows for x in range(int(r[1]), int(r[2]), step)]

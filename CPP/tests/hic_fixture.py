#!/usr/bin/env python3
"""Independent .hic v9 fixture: normalized cis O/E and trans SLICE end to end."""
import math
from pathlib import Path
import struct
import subprocess
import sys
import tempfile
import zlib


def pack(fmt, *values):
    return struct.pack('<' + fmt, *values)


def string(value):
    return value.encode() + b'\0'


def fixture(path, norms=('SCALE', 'INTER_SCALE', 'GW_SCALE'), resolution=100000, missing_pair=False, shift=0, exact_length=False, chromosome_count=3, trans_value=None):
    # Non-integral chromosome lengths exercise the historical floor(length/r)+1
    # bin convention without relying on an external writer or downloaded data.
    n = 120
    length = n * resolution - (0 if exact_length else 1)
    chroms = [('All', 36000)] + [(str(c), length) for c in range(1, chromosome_count + 1)]
    data = bytearray(b'HIC\0' + pack('iQ', 9, 0) + string('synthetic'))
    nvi_patch = len(data)
    data += bytes(16) + pack('ii', 0, len(chroms))
    for name, length in chroms:
        data += string(name) + pack('Q', length)
    data += pack('iii', 1, resolution, 0)
    # Readers can safely prefetch header bytes.
    data += bytes(131072 - len(data))
    matrices = []
    for a in range(1, chromosome_count + 1):
        for b in range(a, chromosome_count + 1):
            if missing_pair and (a, b) == (1, 2):
                continue
            rows = []
            total = 0
            count = 0
            for y in range(n):
                points = []
                for x in range(y + 1) if a == b else range(n):
                    # Six contiguous cis domains; inter contacts follow two
                    # broader states with variation within each state.
                    if a == b:
                        value = 70 if x // 20 == y // 20 else 3
                        value += (x * 7 + y * 3 + a + shift) % 9
                    else:
                        value = 55 if (x // 20) % 2 == (y // 20) % 2 else 8
                        value += (x * 3 + y * 5 + a + b + shift) % 13
                        if trans_value is not None:
                            value = trans_value(a, b, x, y)
                    points.append(pack('hf', x, value))
                    total += value
                    count += 1
                rows.append(pack('hh', y, len(points)) + b''.join(points))
            # Sparse row representation, float counts, short X/Y coordinates.
            raw = pack('iii4Bh', count, 0, 0, 1, 0, 0, 1, n) + b''.join(rows)
            block = zlib.compress(raw, level=0)
            block_at = len(data)
            data += block
            matrix_at = len(data)
            matrix = pack('iii', a, b, 1) + string('BP') + pack('i4f4i', 0, total, count, 0, 0, resolution, n, 1, 1)
            matrix += pack('iQi', 0, block_at, len(block))
            data += matrix
            matrices.append((f'{a}_{b}', matrix_at, len(matrix)))
    master = len(data)
    struct.pack_into('<Q', data, 8, master)
    directory = pack('i', len(matrices))
    for name, at, length in matrices:
        directory += string(name) + pack('Qi', at, length)

    def expected(norm=None):
        return ((string(norm) if norm else b'') + string('BP') + pack('iQ', resolution, n)
                + pack(f'{n}f', *([10.] * n)) + pack('i', chromosome_count)
                + b''.join(pack('if', c, 1.) for c in range(1, chromosome_count + 1)))

    footer = directory + pack('i', 1) + expected()
    data += pack('Q', len(footer)) + footer + pack('i', len(norms))
    for norm in norms:
        data += expected(norm)
    nvi = len(data)
    index_size = 4 + sum(len(string(norm)) + 4 + 3 + 4 + 8 + 8 for norm in norms for _ in range(chromosome_count))
    index, payload = bytearray(pack('i', len(norms) * chromosome_count)), bytearray()
    for norm in norms:
        for c in range(1, chromosome_count + 1):
            # Different vectors make incorrect normalization choices observable.
            value = 2. if norm in ('SCALE', 'KR') else 3.
            vector = pack('Q', n) + pack(f'{n}f', *([value] * n))
            index += string(norm) + pack('i', c) + string('BP') + pack('iQQ', resolution, nvi + index_size + len(payload), len(vector))
            payload += vector
    assert len(index) == index_size
    data += index + payload + bytes(131072)
    struct.pack_into('<QQ', data, nvi_patch, nvi, len(index))
    path.write_bytes(data)


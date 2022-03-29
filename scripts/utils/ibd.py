import numpy
import pandas
import os
import gzip


class Segment:
    def __init__(self, id1, id2, chrom, cm_start=None, cm_end=None, bp_start=None, bp_end=None):
        self.id1 = id1
        self.id2 = id2
        self.chrom = chrom
        self.cm_start = cm_start
        self.cm_end = cm_end
        if self.cm_start is not None and self.cm_end is not None:
            self.cm_len = cm_end - cm_start

        self.bp_start = bp_start
        self.bp_end = bp_end
        if self.bp_start is not None and self.bp_end is not None:
            self.bp_len = self.bp_end - self.bp_start

    def overlap(self, other):
        if self.cm_start is None or self.cm_end is None:
            raise ValueError('Segment is not interpolated, please use "interpolate(map)"')

        if self.id1 != other.id1 or self.id2 != other.id2:
            return 0

        if self.chrom != other.chrom:
            return 0

        if self.cm_start >= other.cm_end or self.cm_end <= other.cm_start:
            return 0

        left = max(self.cm_start, other.cm_start)
        right = min(self.cm_end, other.cm_end)

        return right - left

    def get_key(self):
        return tuple(sorted((self.id1, self.id2)))


def interpolate_all(segments: dict, maps_dir: str) -> dict:
    chromosomes = list(range(1, 23))
    cm_maps = {}
    for chrom in chromosomes:
        cm_map_path = os.path.join(maps_dir, f'chr{chrom}.cm.map')
        cm_map = read_cm_map(cm_map_path)
        cm_maps[chrom] = cm_map

    if len(segments) == 0:
        return {}
    # here we group ALL segments by chromosome for faster interpolation
    chr_segments = {chrom: [] for chrom in chromosomes}
    for segs in segments.values():
        for s in segs:
            chr_segments[s.chrom].append(s)

    interpolated = []
    for chrom, segs in chr_segments.items():
        cm_map = cm_maps[chrom]
        interpolated.extend(interpolate(segs, cm_map))

    int_segments = {key: [] for key in segments.keys()}
    for seg in interpolated:
        int_segments[seg.get_key()].append(seg)

    return int_segments


def interpolate(segments, cm_map):

    bp_starts = numpy.array([seg.bp_start for seg in segments])
    bp_ends = numpy.array([seg.bp_end for seg in segments])
    if None in bp_starts:
        raise ValueError('Cannot interpolate segments with bp_start=None')
    if None in bp_ends:
        raise ValueError('Cannot interpolate segments with bp_end=None')

    cm_starts = numpy.interp(bp_starts, cm_map.bp_pos.values, cm_map.cm_pos.values)
    cm_ends = numpy.interp(bp_ends, cm_map.bp_pos.values, cm_map.cm_pos.values)
    # print('cm_map max is ', cm_map.cm_pos.values.max())
    for seg, cm_start, cm_end in zip(segments, cm_starts, cm_ends):
        seg.cm_start = cm_start
        seg.cm_end = cm_end
        seg.cm_len = cm_end - cm_start

    return segments


def read_cm_map(path):
    # 1       1:249233056_C_T 293.397 249233056
    cm_map = pandas.read_table(path, header=None, names=['chrom', 'snp', 'cm_pos', 'bp_pos'])
    return cm_map


def read_king_segments(path):
    # FID1    ID1     FID2    ID2     IBDType Chr StartMB StopMB StartSNP StopSNP N_SNP Length
    data = pandas.read_table(path, compression='gzip', dtype={'FID1': str, 'ID1': str, 'FID2': str, 'ID2': str})
    segments = {}
    for i, row in data.iterrows():
        id1 = row['FID1'] + '_' + row['ID1']
        id2 = row['FID2'] + '_' + row['ID2']
        seg = Segment(id1, id2, row['Chr'],
                      bp_start=row['StartMB']*1e+6, bp_end=row['StopMB']*1e+6)
        key = tuple(sorted((seg.id1, seg.id2)))
        if key not in segments:
            segments[key] = [seg]
        else:
            segments[key].append(seg)

        if i < 5:
            print(f'segment number {i}:', seg.id1, seg.id2, seg.chrom, seg.bp_start, seg.bp_end)

    print(f'len of segments is {len(segments)}')
    return segments


def read_pedsim_segments(path: str) -> dict:
    # first1_g1-b1-s1	first1_g2-b1-i1	1	752721	249170711	IBD1	0.000000	261.713366	261.713366
    data = pandas.read_table(path, header=None, names=['id1', 'id2', 'chrom', 'gen_start', 'gen_end', 'ibd_type', 'cm_start', 'cm_end', 'cm_len'])
    segments = {}
    for i, row in data.iterrows():
        seg = Segment(row['id1'], row['id2'], row['chrom'],
                      cm_start=row['cm_start'], cm_end=row['cm_end'],
                      bp_start=row['gen_start'], bp_end=row['gen_end'])
        key = tuple(sorted((seg.id1, seg.id2)))
        if key not in segments:
            segments[key] = [seg]
        else:
            segments[key].append(seg)

    return segments


def read_germline_segments(path):
    # 0 first1_g2-b1-i1       0 first1_g2-b2-i1       1       102586656 120506032     1:102586656:T:C 1:120506032:G:A 4416019.386  cM      0       0       0
    germline_names = [
        'fid_iid1',
        'fid_iid2',
        'chrom',
        'start_end_bp',
        'start_end_snp',
        'snp_len',
        'genetic_len',
        'len_units',
        'mismatches',
        'is_homozygous1',
        'is_homozygous2'
    ]

    data = pandas.read_table(path, header=None, names=germline_names)
    segments = {}
    for i, row in data[data.len_units == 'cM'].iterrows():

        id1 = row['fid_iid1'].split()[1]
        id2 = row['fid_iid2'].split()[1]

        bp_start, bp_end = [int(b) for b in row['start_end_bp'].split()]
        seg = Segment(id1, id2, row['chrom'], bp_start=bp_start, bp_end=bp_end)
        key = tuple(sorted((seg.id1, seg.id2)))
        if key not in segments:
            segments[key] = [seg]
        else:
            segments[key].append(seg)

    return segments


def read_refined_ibd_segments(path):
    # first3_g5-b2-s1	1	first3_g7-b3-i1	2	20	12227899	13105474	19.45	1.788
    refined_ibd_names = [
        'fid_iid1',
        'haplotype1',
        'fid_iid2',
        'haplotype2',
        'chrom',
        'bp_start',
        'bp_end',
        'lod_score',
        'genetic_len'
    ]

    data = pandas.read_table(path, header=None, names=refined_ibd_names)
    segments = {}
    for i, row in data.iterrows():

        id1 = row['fid_iid1']
        id2 = row['fid_iid2']

        seg = Segment(id1, id2, row['chrom'], bp_start=row['bp_start'], bp_end=row['bp_end'])
        key = tuple(sorted((seg.id1, seg.id2)))
        if key not in segments:
            segments[key] = [seg]
        else:
            segments[key].append(seg)

    return segments


def total_overlap(true_segments, found_segments):

    overlaps = {}
    for key, true_segs in true_segments.items():
        overlap = 0.0
        true_len = 0.0
        found_len = 0.0

        if key not in found_segments:
            true_len += sum([seg.cm_len for seg in true_segs])
            continue

        found_segs = found_segments[key]
        for ts in true_segs:
            for fs in found_segs:
                overlap += ts.overlap(fs)

        true_len += sum([seg.cm_len for seg in true_segs])
        found_len += sum([seg.cm_len for seg in found_segs])
        overlaps[key] = {'overlap': overlap, 'true_len': true_len, 'found_len': found_len}

    return overlaps


def segments_to_germline(segments, dest_path):

    records = []
    for key, segs in segments.items():
        for seg in segs:
            record = {
                'fid_iid1': f'{key[0]} {key[0]}',
                'fid_iid2': f'{key[1]} {key[1]}',
                'chrom': seg.chrom,
                'start_end_bp': f'{seg.bp_start} {seg.bp_end}',
                'start_end_snp': f'{seg.bp_start} {seg.bp_end}',  # ersa should not use it
                'snp_len': 1,  # ersa should not use it
                'genetic_len': seg.cm_len,
                'len_units': 'cM',
                'mismatches': 0,  # ersa should not use it
                'is_homozygous1': 0,  # ersa should not use it
                'is_homozygous2': 0  # ersa should not use it
            }
            records.append(record)

    frame = pandas.DataFrame.from_records(records)
    frame.to_csv(dest_path, sep='\t', index=False, header=None)
    return len(records)


def merge_germline_segments(segments, gap=0.6):

    new_segments = {key: [] for key in segments.keys()}

    for key, segs in segments.items():

        merged = None
        sorted_segs = sorted(segs, key=lambda s: (s.chrom, s.cm_start))
        for i, seg in enumerate(sorted_segs):
            if merged is None:
                merged = seg
                continue

            if merged.chrom == seg.chrom and seg.cm_start - merged.cm_end <= gap:

                merged = Segment(seg.id1, seg.id2, seg.chrom,
                                 cm_start=merged.cm_start, cm_end=seg.cm_end, bp_start=merged.bp_start, bp_end=seg.bp_end)

            else:
                new_segments[key].append(merged)
                merged = seg

    return new_segments


def line_generator(matchfiles):
    for i in matchfiles:
        iterator = gzip.open(i, 'rt')
        for line in iterator:
            yield line

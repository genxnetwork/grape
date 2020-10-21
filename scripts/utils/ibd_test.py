import ibd
import tempfile
import pandas


def test_interpolate():
    seg1 = ibd.Segment('id1', 'id2', '1', bp_start=10000, bp_end=100000)
    seg2 = ibd.Segment('id1', 'id2', '1', bp_start=10000, bp_end=705882)
    cm_map = ibd.read_cm_map('data\\chr1.cm.map')
    interpolated = ibd.interpolate([seg1, seg2], cm_map)
    assert interpolated[0].cm_len == 0.0
    assert interpolated[1].cm_len >= 0.35
    assert interpolated[1].cm_len <= 0.40


def test_read_germline_segments():
    path = 'data\\all_err4.tsv'
    segments = ibd.read_germline_segments(path)
    assert len(segments) == 634
    total = []
    for segs in segments.values():
        total.extend(segs)
    assert len(total) == 39901


def test_segment_overlap():
    seg1 = ibd.Segment('id1', 'id2', '1', cm_start=10.0, cm_end=15.0)
    seg2 = ibd.Segment('id1', 'id2', '1', cm_start=5.0, cm_end=10.0)
    assert seg1.overlap(seg2) == 0.0

    seg2 = ibd.Segment('id1', 'id2', '1', cm_start=5.0, cm_end=11.0)
    assert seg1.overlap(seg2) == 1.0

    seg2 = ibd.Segment('id1', 'id2', '1', cm_start=10.0, cm_end=16.0)
    assert seg1.overlap(seg2) == 5.0

    seg2 = ibd.Segment('id1', 'id2', '1', cm_start=10.0, cm_end=15.0)
    assert seg1.overlap(seg2) == 5.0

    seg2 = ibd.Segment('id1', 'id2', '1', cm_start=13.0, cm_end=20.0)
    assert seg1.overlap(seg2) == 2.0


def test_total_overlap():
    seg11 = ibd.Segment('id1', 'id2', '1', cm_start=10.0, cm_end=15.0)
    seg12 = ibd.Segment('id1', 'id2', '2', cm_start=20.0, cm_end=30.0)

    seg21 = ibd.Segment('id2', 'id3', '1', cm_start=10.0, cm_end=15.0)
    seg22 = ibd.Segment('id2', 'id3', '2', cm_start=20.0, cm_end=30.0)

    true_seg11 = ibd.Segment('id1', 'id2', '1', cm_start=10.5, cm_end=14.5)
    true_seg12 = ibd.Segment('id1', 'id2', '2', cm_start=25.0, cm_end=31.0)

    true_seg21 = ibd.Segment('id2', 'id3', '1', cm_start=10.5, cm_end=14.5)
    true_seg22 = ibd.Segment('id2', 'id3', '2', cm_start=25.0, cm_end=31.0)

    true_segments = {('id1', 'id2'): [true_seg11, true_seg12], ('id2', 'id3'): [true_seg21, true_seg22]}
    found_segments = {('id1', 'id2'): [seg11, seg12], ('id2', 'id3'): [seg21, seg22]}

    overlaps = ibd.total_overlap(true_segments, found_segments)

    assert overlaps[('id1', 'id2')]['overlap'] == 9.0
    assert overlaps[('id2', 'id3')]['overlap'] == 9.0


def test_segments_to_germline():
    # empty
    records = []
    frame = pandas.DataFrame.from_records(records)
    with tempfile.NamedTemporaryFile('w') as tf:
        frame.to_csv(tf, sep='\t', index=False, header=None)

import tempfile
from pathlib import Path

import transform_ibis_segments


def test_output_produced():
    with tempfile.NamedTemporaryFile('w') as tf:
        transform_ibis_segments.transform('test_data/ibis_merged.seq', tf.name)
        with open(tf.name, 'r') as f:
            assert 'first1_g1-b1-i1' in f.read()


def test_empty_input_ignored():
    with tempfile.NamedTemporaryFile('w') as tf:
        tf.close()
        assert not Path(tf.name).exists()
        transform_ibis_segments.transform('test_data/empty.seq', tf.name)
        assert Path(tf.name).exists()

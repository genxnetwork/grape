import tempfile
from pathlib import Path

import postprocess_ersa


def test_output_produced():
    with tempfile.NamedTemporaryFile('w') as tf:
        postprocess_ersa.postprocess(
            ibd_path='test_data/ibis_merged.seq',
            ersa_path='test_data/ersa_relatives.tsv',
            output_path=tf.name
        )
        with open(tf.name, 'r') as f:
            assert 'first1_g1-b1-i1' in f.read()


def test_empty_input_ignored():
    with tempfile.NamedTemporaryFile('w') as tf:
        tf.close()
        assert not Path(tf.name).exists()
        postprocess_ersa.postprocess(
            ibd_path="test_data/empty.seq",
            ersa_path="test_data/ersa_relatives.tsv",
            output_path=tf.name
        )
        assert Path(tf.name).exists()

import relatives
import os


def test_read_ersa():
    ersa = relatives.read_ersa(os.path.join('test_data', 'test_ersa_output.tsv'))
    assert len(ersa.edges) == 5
    edge1 = ersa.get_edge_data('first1_g5-b3-i1', 'first1_g6-b6-i1')
    assert edge1['Rel_est1'] == 'Sibling'
    assert edge1['Rel_est2'] == 'Sibling'
    assert edge1['d_est'] == '1'

    edge_na = ersa.get_edge_data('first1_g1-b1-i1', 'first1_g1-b4-i1', default={'NA': True})
    assert edge_na['NA']


def test_read_across_kin():
    king = relatives.read_across_kin(os.path.join('test_data', 'test_king_output.seg'))
    assert len(king.edges) == 5
    edge1 = king.get_edge_data('first1_g1-b1-s1', 'first1_g2-b3-i1')
    assert edge1['degree'] == 'PO'

    edge_na = king.get_edge_data('first1_g5-b3-i1', 'first1_g6-b3-i1', default={'NA': True})
    assert edge_na['NA']
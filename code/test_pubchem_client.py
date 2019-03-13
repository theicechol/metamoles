
import pandas as pd

from pandas.util.testing import assert_frame_equal

import pubchem_client


def test_sid_to_smiles():
    """Unit test for pubchem_client.py sid_to_smiles."""

    sids = ['3489', '3990']
    expected = ['C(CO)N', 'C1CSSC1CCCCC(=O)O']
    actual = []

    for sid in sids:
        result_smile = pubchem_client.sid_to_smiles(sid)

        assert len(
            result_smile) >= 1, 'SMILES string is very short. Check SMILES.'
        isinstance(result_smile, str), 'SMILES not returned as string.'

        actual.append(result_smile[0])

    assert expected == actual, 'Actual SMILES are not the expected SMILES.'

    return


def test_kegg_df_to_smiles():
    """Unit test for pubchem_client.py kegg_df_to_smiles."""

    test_frame = pd.DataFrame([['space fill', 'ethanolamine', '3489'], [
                              'space fill', 'pyruvate', '3324']], columns=['Filler', 'Compound Name', 'SID'])

    expected_frame = pd.DataFrame([['space fill',
                                    'ethanolamine',
                                    '3489',
                                    int(700),
                                    'C(CO)N'],
                                   ['space fill',
                                    'pyruvate',
                                    '3324',
                                    int(1060),
                                    'CC(=O)C(=O)O']],
                                  columns=['Filler',
                                           'Compound Name',
                                           'SID',
                                           'CID',
                                           'SMILES'])

    result_frame = pubchem_client.kegg_df_to_smiles(test_frame)

    assert_frame_equal(
        result_frame[0], expected_frame), 'Did not generate expected df.'

    return
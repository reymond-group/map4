
import os
import pickle

from rdkit import Chem
from .map4 import MAP4Calculator
from .PDGA_class_natural import PDGA
from mhfp.encoder import MHFPEncoder
import tmap as tm


def calc_mhfp(ENC, mol):
    """calculates the minhashed fingerprint

    Arguments:
        mol

    Returns:
        tmap VectorUint -- minhashed fingerprint
    """
    
    smiles = Chem.MolToSmiles(mol)
    fp = tm.VectorUint(ENC.encode(smiles))
    return fp


def calc_map(mol):
    """calculates the atom pair minhashed fingerprint

    Arguments:
        mol

    Returns:
        tmap VectorUint -- minhashed fingerprint
    """

    MAP4 = MAP4Calculator(dimensions=512)
    fp = MAP4.calculate(mol)
    
    return fp


def similarity_search(chos_fp, db, mol, number_of_hits, lf1, lf2, lf3, findID1, findID2, findID3):
    """returns n hits of the given query ligand  

    Arguments:
        smile {string} -- smile of the query ligand
        number_of_hits {integer} -- number of required hits

    Returns:
        list -- n NN of the query molecule according to MXfp 
    """

    results = []

    if db == 'ChEMBL':
        lf = lf1
        findID = findID1

    elif db == 'SwissProt':
        lf = lf2
        findID = findID2

    else:
        lf = lf3
        findID = findID3

    mhfp_encoder = MHFPEncoder(512)

    if chos_fp == 'MAP4':
        fp = calc_map(mol)
    else:
        fp = calc_mhfp(mhfp_encoder, mol)

    NNs = lf.query_linear_scan(fp, int(number_of_hits))


    for i, NN in enumerate(NNs):
        results.append([findID[NN[1]][1].split(';'), findID[NN[1]]
                        [0], round(NN[0], 3), findID[NN[1]][2]])

    return results


def seq_to_smiles(seq):
    ga = PDGA()
    smiles = ga.smiles_from_seq(seq)
    return smiles

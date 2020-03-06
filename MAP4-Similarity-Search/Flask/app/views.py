import os
import pickle
import time
from multiprocessing import Pool

import tmap as tm
from app import app
from flask import (Flask, escape, flash, redirect, render_template, request,
                                      session, url_for)
from flask_restful import reqparse
from flask_wtf import FlaskForm
from mhfp.encoder import MHFPEncoder
from rdkit import Chem
from rdkit.Chem import AllChem
from wtforms import (Form, StringField, SubmitField, TextAreaField, TextField,
                                          validators)
from wtforms.validators import DataRequired

from .functs import seq_to_smiles, similarity_search


data_path_prod = '/MAP-SearchData/'
data_path_dev = '/data/MHAPSearchData/'
if os.path.exists(data_path_prod):
    data_path = data_path_prod
else:
    data_path = data_path_dev

print(data_path)
lf1 = tm.LSHForest(512, 32)
lf2 = tm.LSHForest(512, 32)
lf3 = tm.LSHForest(512, 32)
lf4 = tm.LSHForest(512, 32)
lf5 = tm.LSHForest(512, 32)
lf6 = tm.LSHForest(512, 32)

lf1.restore(data_path + 'chembl_25_wAct_comp.mhap4-512_LSHforest')
findID1 = pickle.load(
    open(data_path + 'chembl_25_wAct_comp.mhap4-512_dictionary', 'rb'))

lf2.restore(data_path + 'uniprot_sprot_smiles_compact_mhap4-512_LSHforest')
findID2 = pickle.load(
    open(data_path + 'uniprot_sprot_smiles_compact_mhap4-512_dictionary', 'rb'))

lf3.restore(data_path + 'metabolome_compact.mhap4-512_LSHforest')
findID3 = pickle.load(
    open(data_path + 'metabolome_compact.mhap4-512_dictionary', 'rb'))

lf4.restore(data_path + 'chembl_25_wAct_comp.mhfp6_LSHforest')
findID4 = pickle.load(
    open(data_path + 'chembl_25_wAct_comp.mhfp6_dictionary', 'rb'))

lf5.restore(data_path + 'uniprot_sprot_smiles_compact_mhfp6_LSHforest')
findID5 = pickle.load(
    open(data_path + 'uniprot_sprot_smiles_compact_mhfp6_dictionary', 'rb'))

lf6.restore(data_path + 'metabolome4.smi_compact_mhfp6_LSHforest')
findID6 = pickle.load(
    open(data_path + 'metabolome4.smi_compact_mhfp6_dictionary', 'rb'))


@app.route('/', methods=['GET'])
def index():
    """passes the data to the index HTML page

    Returns:
        render_template -- template for index.html
    """

    return render_template("index.html", DB=['ChEMBL', 'Metabolome', 'SwissProt'], FP=['MAP4', 'MHFP6'])


@app.route('/wait', methods=['POST', 'GET'])
def wait():
    """gets the data from the user choices in index.html and 
        passes them to results.html: database, query SMILES and number of required NNs

    Returns:
        render_template -- wait.html
    """

    body = render_template('wait.html')

    db = request.form['db']
    fp = request.form['fp']
    smiles = request.form['smiles']
    seq = request.form['seq']
    hits = request.form['hits']

    if smiles == '' and seq != '':
        smiles = seq_to_smiles(seq)
    elif smiles == '' and seq == '':
        smiles = 'c1ccccc1'

    session['db'] = db
    session['fp'] = fp
    session['smiles'] = smiles
    session['hits'] = hits

    predurl = "1; url=results"

    return (body, 200, {'Refresh': predurl})


@app.route('/results', methods=['POST', 'GET'])
def results():
    """gets the data from wait.html, calculates the NNs, 
    passes the results and the number of NN to results.html.

    Returns:
        render_template -- results.html
    """
    db = session.get('db', None)
    fp = session.get('fp', None)
    smiles = session.get('smiles', None)
    hits = session.get('hits', None)

    smiles_list = smiles.split('.')
    if len(smiles) > 1:
        smiles = max(smiles_list, key=len)

    mol = Chem.MolFromSmiles(smiles)
    if mol == None:
        return render_template("rdkiterror.html")
    smiles_tmp = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=False)
    mol = Chem.MolFromSmiles(smiles_tmp)

    if fp == 'MAP4':
        results = similarity_search(
            fp, db, mol, hits, lf1, lf2, lf3, findID1, findID2, findID3)
    else:
        results = similarity_search(
            fp, db, mol, hits, lf4, lf5, lf6, findID4, findID5, findID6)

    if results == []:
        return render_template('results.html', results=[],
                               hits=0, query=smiles, db=db, fp=fp, ID=[], url=[],
                               title='There are not Near Neighbours in the ' + db + ' Space')

    results_ = []

    if db == 'SwissProt':
        for result in results:
            result_ = result
            ID = [result[3].split('|')[1]]

            add_info = 'Sequence: ' + result[0][0]
            result_.append(ID)
            result_.append(add_info)
            results_.append(result_)
        url = 'https://www.uniprot.org/uniprot/'

    elif db == 'Metabolome':
        for result in results:
            result_ = result
            ID = result[0]
            add_info = ''
            for source in result[3].split(';'):
                if source != 'NA':
                    add_info += ' ' + source
            if add_info != '':
                info = 'Source annotation:' + add_info
            else:
                info = ''
            result_.append(ID)
            result_.append(info)
            results_.append(result_)
        url = 'http://www.hmdb.ca/metabolites/'

    else:
        for result in results:
            result_ = result
            ID = result[0]
            add_info = ''

            if result[3] != 'NA':
                for target in result[3].split(';'):
                    target = target.split(',')
                    if target[0] not in add_info:
                        add_info += ' ' + target[0]

            if add_info != '':
                info = 'Target annotation:' + add_info
            else:
                info = ''
            result_.append(ID)
            result_.append(info)
            results_.append(result_)
        url = 'https://www.ebi.ac.uk/chembl/compound/inspect/'

    return render_template('results.html', results=results_,
                                                    hits=str(len(results)), query=smiles, db=db, fp=fp, ID=ID, url=url,
                           title='the ' + str(len(results)) + ' Nearest Neighbours in the ' + fp + ' ' + db + ' Space')

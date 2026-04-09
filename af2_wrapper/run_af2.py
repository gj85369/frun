'''
Author : Alphafold / Usman Ghani
Date : 11/3/2021

This script is based on alphafold notebook for multimer prediction:
    https://github.com/deepmind/alphafold/blob/main/notebooks/AlphaFold.ipynb
'''
# The code below is based on AlphaFold Multimer Colab implementation.

# Copyright 2021 DeepMind Technologies Limited
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import os
import sys
import csv
import json
import pickle
import warnings
import argparse
import numpy as np
import random
import pathlib
from absl import logging
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
import tensorflow as tf
warnings.filterwarnings('ignore')
logging.set_verbosity("error")
tf.get_logger().setLevel('ERROR')


AF2_PATH = os.getenv('ALPHAFOLD_MULTIMER_PATH')
AF2_PARAMS_PATH = os.getenv('ALPHAFOLD_PARAMS_PATH')
#Add alphafold to PATH
sys.path.append(pathlib.Path(__file__).parent.resolve())
PATHs = sys.path
if AF2_PATH not in PATHs:
    sys.path.append(AF2_PATH)

from alphafold.data import feature_processing

#from alphafold.model import config
from alphafold.model import model
#from alphafold.model import data

from alphafold.common import protein
#from alphafold.common import confidence
#from alphafold.common.confidence import predicted_tm_score

from af2tmplt import global_var
from af2tmplt.inputs import process_inputs
from af2tmplt.msa_mmseq import get_MSAs_mmseq
from af2tmplt.msa_std import get_MSAs_std
from af2tmplt.msa_std import get_MSAs_googleapi
from af2tmplt.msa_std import get_MSAs_single_seq
from af2tmplt.featurize import build_feature_dict_default

from af2tmplt.scoring import get_confidence_metrics
from af2tmplt.predict import predict_structure_tmplt

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Predict structure.')
    parser.add_argument('input_json', type=str, help='json file containing job parameters.')
    parser.add_argument('-mp', '--msa_path', default="./")
    parser.add_argument('-fo', '--features_only', action='store_true')
    parser.add_argument('-ff', '--from_features', action='store_true')
    args = parser.parse_args()

    with open(args.input_json, 'r') as f:
        inputs = json.load(f)

    query_seq_dict, chain_names = process_inputs(inputs)
    templates = inputs.get('templates', [None])
    #Set the global variable for possible chain names. The hyphen at the start is to avoid using the
    #0 index for the 'asym_id' feature going forward. It seems that alphafold used 0 for padding for a lot of its features and
    #so they themselves have avoided using 0s
    protein.PDB_CHAIN_IDS = ('-') + ''.join(chain_names)
    #feature_processing.MAX_TEMPLATES = 1
    model.get_confidence_metrics = get_confidence_metrics


    if not args.from_features:
        if inputs['msa'] == 'standard':
            query_seq_to_msa = get_MSAs_std(query_seq_dict,
                                            AF2_PARAMS_PATH,
                                            hhblits_binary_path=inputs['hhblits_path'],
                                            hmmer_binary_path=inputs['hmmer_path'],
                                            output_path=args.msa_path)
        elif inputs['msa'] == 'google_api':
            query_seq_to_msa = get_MSAs_googleapi(query_seq_dict,
                                                  output_path=args.msa_path)
        elif inputs['msa'] == 'mmseqs':
            query_seq_to_msa = get_MSAs_mmseq(query_seq_dict,
                                              output_path=args.msa_path)
        elif inputs['msa'] == 'single_seq':
            query_seq_to_msa = get_MSAs_single_seq(query_seq_dict,
                                                   output_path=args.msa_path)
        else:
            print("unknown MSA parameter value")
            exit()

        for i, template in enumerate(templates):
            #print(i, template['pdb_file'])
            feature_dict =  build_feature_dict_default(
                                inputs['sequences'],
                                query_seq_to_msa,
                                allow_msa_pairing = False,
                                msa_crop_size=inputs['msa_crop_size'],
                                custom_template = template)
            file_suffix = '' if template == None else '.tmplt_{}'.format(i)
            features_fname = 'features{}.pkl'.format(file_suffix)
            with open(features_fname, 'wb') as f:
                pickle.dump(feature_dict, f)
        if args.features_only:
            exit()

    for i, template in enumerate(templates):
        template_name = None if template == None else template['pdb_file']
        file_suffix = '' if template == None else '.tmplt_{}'.format(i)
        features_fname = 'features{}.pkl'.format(file_suffix)
        with open(features_fname, 'rb') as f:
            feature_dict = pickle.load(f)
        predict_structure_tmplt(feature_dict,
                                inputs['job_name'] + file_suffix,
                                AF2_PARAMS_PATH,
                                number_of_models=inputs['num_models'],
                                number_of_recycles=inputs['num_recycle'],
                                number_of_predictions=inputs['num_predictions'],
                                number_of_msa_clusters=inputs['num_msa_clusters'],
                                seed=inputs['seed'],
                                custom_template_name=template_name,
#                                merge_chains=False)
                                merge_chains=True)
                                

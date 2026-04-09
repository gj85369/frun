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
import csv
import json
import pickle
import warnings
import argparse
import numpy as np
import scipy.special
import random
import copy

from typing import Optional

from alphafold.data import feature_processing

from alphafold.model import config
from alphafold.model import model
from alphafold.model import data

from alphafold.common import protein
from alphafold.common import confidence
#from alphafold.common.confidence import predicted_tm_score
from alphafold.common.confidence import _calculate_bin_centers

from af2tmplt import global_var
from af2tmplt.inputs import process_inputs
from af2tmplt.featurize import build_feature_dict_default
class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, (np.floating, np.complexfloating)):
            return float(obj)        
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return super().default(obj)


def _renum_residues(residue_index):
    # For multi-chain targets ran through monomer parameters
    # and multi-chain targets ran throuh multimer parameters
    # (particularly important for targets with templates),
    # reindex residue ids as a single chain with breaks.
    residue_index_modif = residue_index.copy()
    #'''
    asym_id_values_set = set(global_var.ASYM_IDs)
    tot_chain_len = 0
    for i in sorted(list(asym_id_values_set)):
        chain_len = np.sum(global_var.ASYM_IDs[:] == i)
        tot_chain_len += chain_len

        if 'monomer' in global_var.MODEL_SET:
            residue_index_modif[:,tot_chain_len:] += chain_len + 200
        else:
            residue_index_modif[tot_chain_len:] += chain_len + 200
    #'''
    return residue_index_modif

def modify_processed_features(processed_feature_dict,
                              np_example,
                              merge_chains=True):
    # Save the 'asym_id' features in a global variable.
    # This way they can accessed from the get_confidence_metrics() function.
    if 'monomer' in global_var.MODEL_SET:
        global_var.ASYM_IDs = np_example['asym_id'][0]
    else:
        global_var.ASYM_IDs = np_example['asym_id']

    # For multimer parameters,
    # create a modified version of the processed features,
    # where the 'entity_id', 'sym_id' and 'asym_id' features are the same for
    # all chains (The whole system is one big pseudochain with breaks)
    processed_feature_dict_modif = {}
    if 'asym_id' in np_example and merge_chains:
        chain_id_features = ['entity_id', 'sym_id', 'asym_id']
    else:
        chain_id_features = []
    for k, v in processed_feature_dict.items():
        if k in chain_id_features:
            #print(k, v.shape, np_example.get(k, np.array([None])).shape)
            processed_feature_dict_modif[k] = np.ones(v.shape)
        elif k in ['residue_index']:
            processed_feature_dict_modif[k] = _renum_residues(v)
        else:
            #print(k, v.shape, np_example.get(k, np.array([None])).shape)
            processed_feature_dict_modif[k] = v

    return processed_feature_dict_modif


def predict_structure_tmplt(np_example,
                            prefix,
                            params_path,
                            number_of_models = 5,
                            number_of_recycles = 3,
                            number_of_predictions = 1,
                            number_of_msa_clusters = 'default', #512 for monomer, 252 for multimer
                            custom_template_name = None,
                            merge_chains = True,
                            seed = 0):
    if merge_chains == False and custom_template_name != None:
        print("[Warning] Using custom template with merge_chains == False.")
        print("[Warning]    If you are handling a multimeric template,")
        print("[Warning]    make sure you are using a custom version of AF2")
        print("[Warning]    which doesn't mask interchain template contacts.")
        print("[Warning]    Otherwise, interchain template orientation")
        print("[Warning]    will be ignored.")
    model_names = global_var.MODEL_PRESETS[global_var.MODEL_SET]

    for model_name in model_names[:number_of_models]:
        print(f'Prediciting {model_name}')
        cfg = config.model_config(model_name)
        if 'monomer' in global_var.MODEL_SET:
            cfg.data.eval.num_ensemble = 1
            cfg.data.common.num_recycle = number_of_recycles
            cfg.model.num_recycle = number_of_recycles
            if number_of_msa_clusters != 'default':
                cfg.data.eval.max_msa_clusters = number_of_msa_clusters
        elif 'multimer' in global_var.MODEL_SET:
            cfg.model.num_ensemble_eval = 1
            cfg.model.num_recycle = number_of_recycles
            if number_of_msa_clusters != 'default':
                cfg.model.num_msa = number_of_msa_clusters
            cfg.model.embeddings_and_evoformer.template.intrachain_mask = False
        params = data.get_model_haiku_params(model_name, params_path)
        model_runner = model.RunModel(cfg, params)

        for predi in range(number_of_predictions):
            model_metrics_output_path = f'{prefix}.{model_name}.{predi}.metrics.pkl'
            jsonmodel_metrics_output_path = f'{prefix}.{model_name}.{predi}.metrics.json'
            
            if os.path.exists(model_metrics_output_path):
                continue

            predi_seed = seed + predi
            processed_feature_dict = model_runner.process_features(
                                            np_example,
                                            random_seed=predi_seed)
            processed_feature_dict_modif = modify_processed_features(
                                            processed_feature_dict,
                                            np_example,
                                            merge_chains=merge_chains)

            prediction = model_runner.predict(processed_feature_dict_modif,
                                          random_seed=predi_seed)
            #prediction = prediction1[0]
            print(type(prediction))
            print(prediction)
            # Get the quality metrix for the prediction
            metric_keys = ['plddt',
                           'predicted_aligned_error',
                           'max_predicted_aligned_error',
                           'ptm',
                           'iptm',
                           'multimer_confidence',
                           'interfaces',
                           'ranking_confidence']
            model_metrics = {}
            for key in metric_keys:
                model_metrics[key] = prediction.get(key, None)

            model_metrics['mean_plddt'] = model_metrics['plddt'].mean()
            mm1 = copy.deepcopy(model_metrics)
            for inst in list(model_metrics['interfaces'].keys()):
                newkey = f'{inst[0]}_{inst[1]}'
                model_metrics['interfaces'][newkey] = model_metrics['interfaces'][inst]
                del model_metrics['interfaces'][inst]
            
            print(f'model metrics keys {list(model_metrics["interfaces"].keys())}')
            #for now, also dump pae logits, asym ids and model coord:
            #model_metrics['raw_pae'] = prediction['raw_pae']
            #model_metrics['raw_atom_positions'] = prediction['structure_module']['final_atom_positions']
            #model_metrics['raw_atom_mask']      = prediction['structure_module']['final_atom_mask']
            #model_metrics['raw_asym_id'] = global_var.ASYM_IDs

            with open(jsonmodel_metrics_output_path, 'w') as f:
                json.dump(model_metrics, f, indent=4, cls=NumpyEncoder)
            f.close()
            with open(model_metrics_output_path, 'wb') as f:
                pickle.dump(mm1, f)
            f.close()
            # Set the b-factors to the per-residue plddt.
            final_atom_mask = prediction['structure_module']['final_atom_mask']
            b_factors = prediction['plddt'][:, None] * final_atom_mask

            # Extract the protein atomic coordinates and save as pdb.
            if 'monomer' in global_var.MODEL_SET:
                # Processing of monomer features does not carry over the 'asym_id'
                # features which are used for correct chain assignment.
                # We add those manually.
                processed_feature_dict['asym_id'] = np.expand_dims(global_var.ASYM_IDs, axis=0)
            unrelaxed_protein = protein.from_prediction(
                                        processed_feature_dict,
                                        prediction,
                                        b_factors=b_factors,
                                        remove_leading_feature_dimension=('monomer' in global_var.MODEL_SET))
            pdb_output_path = f'{prefix}.{model_name}.{predi}.pdb'
            pdb = protein.to_pdb(unrelaxed_protein)
            with open(pdb_output_path, 'w') as f:
                f.write(pdb)

            del prediction

        # Delete unused outputs to save memory.
        del model_runner
        del params

    return(True)

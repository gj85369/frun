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

import pickle
import numpy as np

from typing import MutableMapping

from alphafold.notebooks import notebook_utils
from alphafold.data import pipeline
from alphafold.data import pipeline_multimer
from alphafold.data import msa_pairing
from alphafold.data import feature_processing

from af2tmplt import global_var
from af2tmplt.featurize_template import featurize_custom_template


def cust_pair_and_merge(
    all_chain_features,
    enable_pairing = False
    ):
    """Runs processing on features to augment, pair and merge.
    Args:
        all_chain_features: A MutableMap of dictionaries of features for each chain.
    Returns:
        A dictionary of features.
    """

    feature_processing.process_unmerged_features(all_chain_features)

    np_chains_list = list(all_chain_features.values())

    #pair_msa_sequences = not _is_homomer_or_monomer(np_chains_list)
    pair_msa_sequences = enable_pairing

    if pair_msa_sequences:
        np_chains_list = msa_pairing.create_paired_features(
            chains=np_chains_list)
        np_chains_list = msa_pairing.deduplicate_unpaired_sequences(np_chains_list)
    np_chains_list = feature_processing.crop_chains(
        np_chains_list,
        msa_crop_size=feature_processing.MSA_CROP_SIZE,
        pair_msa_sequences=pair_msa_sequences,
        max_templates=feature_processing.MAX_TEMPLATES)
    np_example = msa_pairing.merge_chain_features(
        np_chains_list=np_chains_list, pair_msa_sequences=pair_msa_sequences,
        max_templates=feature_processing.MAX_TEMPLATES)
    np_example = feature_processing.process_final(np_example)
    return np_example

def cust_pair_and_merge_pseudomonomer(
    features_for_chain,
    enable_pairing = False
    ):

    merged_features = {}
    first_chain = list(features_for_chain.keys())[0]
    for feature_name in features_for_chain[first_chain]:
        #for chain in features_for_chain:
            #print(feature_name, features_for_chain[chain][feature_name].shape)
        if feature_name in ['domain_name', 'num_alignments', 'seq_length']:
            feats = [np.asarray(x[feature_name][0], dtype=x[feature_name].dtype) for x in features_for_chain.values()]
        else:
            feats = [x[feature_name] for x in features_for_chain.values()]
        if feature_name in msa_pairing.MSA_FEATURES:
            if enable_pairing:
                print("Can't pair yet")
                exit()
            else:
                merged_features[feature_name] = msa_pairing.block_diag(
                                                        *feats,
                                                        pad_value = msa_pairing.MSA_PAD_VALUES[feature_name])
        elif feature_name in msa_pairing.SEQ_FEATURES:
            merged_features[feature_name] = np.concatenate(feats, axis=0)
        elif feature_name in msa_pairing.TEMPLATE_FEATURES + ('template_all_atom_masks', 'template_confidence_scores',):
            merged_features[feature_name] = np.concatenate(feats, axis=1)
        elif feature_name in msa_pairing.CHAIN_FEATURES:
            merged_features[feature_name] = np.sum(x for x in feats).astype(np.int32)
        elif feature_name in ('sequence',):
            merged_features[feature_name] = np.array([b''.join([x[0] for x in feats])])
        else:
            merged_features[feature_name] = feats[0]
        #print(feature_name, merged_features[feature_name].shape)

    tot_seq_len = merged_features['seq_length']
    for feature_name in msa_pairing.CHAIN_FEATURES:
        merged_features[feature_name] = np.array([merged_features[feature_name]] * tot_seq_len)

    asym_ids = []
    for i, chain in enumerate(features_for_chain):
        asym_ids.append(np.array([i+1] * features_for_chain[chain]['seq_length'][0]))
    merged_features['asym_id'] = np.expand_dims(np.concatenate(asym_ids), axis=0)
    #print(merged_features['asym_id'])
    return merged_features


def build_feature_dict_default(
                        seq_data_list,
                        seq_to_msa_dict,
                        custom_template=None,
                        allow_msa_pairing = False,
                        msa_crop_size=2048,
                        min_num_seq=512):
    features_for_chain = {}

    if custom_template != None:
        #print("Using template")
        template_features_by_chain_dict = featurize_custom_template(
                                                seq_data_list,
                                                custom_template)
    else:
        #Make a blank for all the chains
        template_features_by_chain_dict = {}
        for seq_data in seq_data_list:
            for chain in seq_data['chains']:
                template_features_by_chain_dict[chain] = notebook_utils.empty_placeholder_template_features(
                    num_templates=0, num_res=len(seq_data['sequence']))

    features_for_chain = {}

    L_prev = 0
    for seq_data in seq_data_list:
        #print(seq_data)
        seq_name = seq_data['name']
        sequence = seq_data['sequence']
        for chain in seq_data['chains']:
            #parsed_MSA = seq_name_to_parsed_MSAs[seq_name].deepcopy()
            unpairable_MSA_list = seq_to_msa_dict[seq_name].get('unpairable', None)[:msa_crop_size]
            pairable_MSA_list   = seq_to_msa_dict[seq_name].get('pairable', None)

            # Turn the raw data into model features.
            feature_dict = {}
            # make_sequence_features creates the following features:
            # 'aatype' - one-hot
            # 'between_segment_residues' 0 x num_res
            # 'domain_name' ???
            # 'residue_index' i for i in range(num_res)
            # 'seq_length' num_res fo i in range(num_res)
            # 'sequence" - sequence in utf-8
            feature_dict.update(pipeline.make_sequence_features(sequence=sequence,
                                                                description='query',
                                                                num_res=len(sequence)))
            # Following the idea by Minkyung Baek:
            # add breaks to residue index to indicate chain breaks
            # (Needed when running with monomer parameters
            #  or with multimer parameters in pseudo-monomer regime)
            ##feature_dict['residue_index'][:] += L_prev
            ##L_prev += feature_dict['residue_index'].shape[0] + int(200)
            # make_msa_features creates the following features:
            # 'deletion_matrix_int'
            # 'msa'
            # 'num_alignments' [msa_depth] x num_res
            # 'msa_species_identifiers'
            feature_dict.update(pipeline.make_msa_features(unpairable_MSA_list))

            #Add template features
            # 'template_aatype' [num_templates, num_res, num_restypes=(20+2)]
            # 'template_all_atom_masks' [num_templates, num_res, atom_tpe_num]
            # 'template_all_atom_positions' [num_templates, num_res, atom_type_num, 3]
            # 'template_domain_names' 0 x num_templates
            # 'template_sequence' 0 x num_templates
            # 'template_sum_probs' 0 x num_templates
            feature_dict.update(template_features_by_chain_dict[chain])

            #if global_var.MULTIMER and not global_var.HOMOMER:
            if allow_msa_pairing and pairable_MSA_list != None:
                # MSA_FEATURES =
                # ['msa', 'msa_mask', 'deletion_matrix', 'deletion_matrix_int']
                valid_feats = msa_pairing.MSA_FEATURES + ('msa_species_identifiers',)
                all_seq_features = {f'{k}_all_seq': v for k, v in pipeline.make_msa_features(pairable_MSA_list).items()
                    if k in valid_feats}
                feature_dict.update(all_seq_features)
            #print(chain, feature_dict['msa'].shape)
            features_for_chain[chain] = feature_dict

    if 'monomer' in global_var.MODEL_SET:
        np_example = cust_pair_and_merge_pseudomonomer(features_for_chain)
    else:
        feature_processing.MSA_CROP_SIZE = msa_crop_size
        # Remove empty leading dimension from:
        #   'sequence', 'domain_name', 'num_alignments', 'seq_length'
        # Convert from one-hot to id for 'aatype'
        # Reparametrize one-hot representation for 'template_aatype'
        # Rename 'template_all_atom_masks' to 'template_all_atom_mask'
        all_chain_features = {}
        for chain_id, chain_features in features_for_chain.items():
            all_chain_features[chain_id] = pipeline_multimer.convert_monomer_features(
                chain_features, chain_id)
            #print(chain_id, all_chain_features[chain_id]['msa'].shape)
        # add_assembly_features does the following:
        # 1. change dictionary keys from {chain_id} to {seq_id}_{sym_id}
        # 2. adds per-residue features 'chain_id', 'sym_id', 'entity_id'
        #all_chain_features = custom_add_assembly_features(all_chain_features)
        all_chain_features = pipeline_multimer.add_assembly_features(all_chain_features)
        #np_example = feature_processing.pair_and_merge(all_chain_features=all_chain_features)
        np_example = cust_pair_and_merge(all_chain_features=all_chain_features,
                                         enable_pairing=allow_msa_pairing)
        #np_example['bert_mask'] = np.ones(np_example['bert_mask'].shape)
        # Pad MSA to avoid zero-sized extra_msa.
        np_example = pipeline_multimer.pad_msa(np_example,
                                               min_num_seq=min_num_seq)

        #with open('test_features.pkl', 'wb') as f:
        #    pickle.dump(np_example, f)
        #exit()

    return(np_example)

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

AF2_PATH = os.getenv('ALPHAFOLD_MULTIMER_PATH')
AF2_PARAMS_PATH = os.getenv('ALPHAFOLD_PARAMS_PATH')
#Add alphafold to PATH
PATHs = sys.path
if AF2_PATH not in PATHs:
    sys.path.append(AF2_PATH)

import csv
import json
import pickle
import warnings
import argparse
import numpy as np
import scipy.special
import random

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
#import global_var

def predicted_tm_score(
    logits: np.ndarray,
    breaks: np.ndarray,
    pair_weights: Optional[np.ndarray] = None) -> np.ndarray:
    """Computes predicted TM alignment or predicted interface TM alignment score.

    The main difference from the standard af2 function is that
    here we fixed the bug in the initialization of pair_mask (wrong dimensions)

 	Args:
      logits: [num_res, num_res, num_bins] the logits output from
        PredictedAlignedErrorHead.
      breaks: [num_bins] the error bins.
      pair_weights: [num_res] the per residue weights to use for the
        expectation.

    Returns:
      ptm_score: The predicted TM alignment or the predicted iTM score.
    """

    bin_centers = _calculate_bin_centers(breaks)

    if pair_weights is None:
        pair_weights = np.ones(shape=(logits.shape[0], logits.shape[0]), dtype=float)
    # here we implicitly assume that residue weights are binary 0/1 values
    num_res = int(np.sum(pair_weights.max(axis=-1)))
    # Clip num_res to avoid negative/undefined d0.
    clipped_num_res = max(num_res, 19)

    # Compute d_0(num_res) as defined by TM-score, eqn. (5) in Yang & Skolnick
    # "Scoring function for automated assessment of protein structure template
    # quality", 2004: http://zhanglab.ccmb.med.umich.edu/papers/2004_3.pdf
    d0 = 1.24 * (clipped_num_res - 15) ** (1./3) - 1.8
    if d0 < 0.5:
        d0 = 0.02*num_res
    # Convert logits to probs.
    probs = scipy.special.softmax(logits, axis=-1)

    # TM-Score term for every bin.
    tm_per_bin = 1. / (1 + np.square(bin_centers) / np.square(d0))
    # E_distances tm(distance).
    predicted_tm_term = np.sum(probs * tm_per_bin, axis=-1)

    predicted_tm_term *= pair_weights

    normed_residue_mask = pair_weights / (1e-8 + np.sum(
        pair_weights, axis=-1, keepdims=True))
    per_alignment = np.sum(predicted_tm_term * normed_residue_mask, axis=-1)
    return np.asarray(per_alignment[(per_alignment * pair_weights.max(axis=-1)).argmax()])

def __get_min_atom_dist(resi_i_pos, resi_j_pos):
    return ((resi_i_pos[None, :] - resi_j_pos[:, None])**2).sum(axis=-1).min()
    #return scipy.spatial.distance.cdist(resi_i_pos, resi_j_pos, metric='sqeuclidean').min()

def __get_contact_map(atom_positions, atom_masks, asym_ids, cutoff = 4.5):

    ''' Slower version for validation
    contact_map2 = np.zeros((atom_positions.shape[0], atom_positions.shape[0]))
    for i in range(contact_map2.shape[0]):
        asym_i = asym_ids[i]
        pos_i = atom_positions[i]
        mask_i = atom_masks[i]
        pos_i = pos_i[mask_i > 0.5]
        for j in range(i+1, contact_map2.shape[1]):
            asym_j = asym_ids[j]
            pos_j = atom_positions[j]
            mask_j = atom_masks[j]
            pos_j = pos_j[mask_j > 0.5]
            if asym_i == asym_j:
                continue

            dist_ij = scipy.spatial.distance.cdist(pos_i, pos_j, metric='sqeuclidean').min()
            if dist_ij < (cutoff*cutoff):
                contact_map2[i,j] = 1.0
                contact_map2[j,i] = 1.0
    '''
    # first, get a rough contact map based on CA atoms
    ca_pos = atom_positions[:,1,:]
    ca_dist_map = ((ca_pos[None, :] - ca_pos[:, None])**2).sum(axis=-1)
    filt_asym_ids = asym_ids[None,:] != asym_ids[:, None]
    prefilt_pairs = np.where((ca_dist_map < 20*20) & filt_asym_ids)
    contact_map = np.zeros((atom_positions.shape[0], atom_positions.shape[0]))

    # go through the identified pairs and refine distance using all atoms
    dist_sq = ((atom_positions[prefilt_pairs[0]][:, None, :, :]
                - atom_positions[prefilt_pairs[1]][:, :, None, :])**2).sum(axis=-1)
    filt_atom_mask = (atom_masks[prefilt_pairs[0]][:, None, :] < 0.5) | (atom_masks[prefilt_pairs[1]][:, :, None] < 0.5)
    dist_sq += 100500 * filt_atom_mask

    if dist_sq.shape[0] != 0:
        dist_sq = dist_sq.min(axis=(-1,-2))
        final_pairs = (prefilt_pairs[0][np.where(dist_sq < cutoff*cutoff)],
                       prefilt_pairs[1][np.where(dist_sq < cutoff*cutoff)])
        contact_map[final_pairs] = 1.0

    #print('diff:', (contact_map - contact_map2).sum())

    return contact_map

def get_confidence_metrics(prediction_result, multimer_mode=True):
    """
    Post processes prediction_result to get confidence metrics.
    Populates 'plddt', 'ptm', 'iptm', 'ranking_confidence' and
    'multimer_confidence' fields.
    The main difference from the standard af2 function is the
    passing of gloabal ASYM_IDs variable to predict_tm_score()
    """
    confidence_metrics = {}
    if 'predicted_lddt' in prediction_result.keys():
        confidence_metrics['plddt'] = confidence.compute_plddt(
            prediction_result['predicted_lddt']['logits'])
    num_chains = len(set(global_var.ASYM_IDs))
    if 'predicted_aligned_error' in prediction_result:
        confidence_metrics.update(confidence.compute_predicted_aligned_error(
            logits=prediction_result['predicted_aligned_error']['logits'],
            breaks=prediction_result['predicted_aligned_error']['breaks']))
        confidence_metrics['ptm'] = predicted_tm_score(
            logits=prediction_result['predicted_aligned_error']['logits'],
            breaks=prediction_result['predicted_aligned_error']['breaks'])
        # Save raw pae to dump it to file later
        confidence_metrics['raw_pae'] = prediction_result['predicted_aligned_error']

        if num_chains > 1: #multimer_mode:
            mask_interchain = global_var.ASYM_IDs[:, None] != global_var.ASYM_IDs[None, :]
            atom_positions = prediction_result['structure_module']['final_atom_positions']
            atom_masks = prediction_result['structure_module']['final_atom_mask']
            ca_dist_mask = __get_contact_map(atom_positions,
                                             atom_masks,
                                             global_var.ASYM_IDs)
            mask_interface =  ca_dist_mask * mask_interchain

            # Compute the ipTM only for the multimer model.
            confidence_metrics['iptm'] = predicted_tm_score(
                logits=prediction_result['predicted_aligned_error']['logits'],
                breaks=prediction_result['predicted_aligned_error']['breaks'],
                pair_weights = mask_interface)
            confidence_metrics['multimer_confidence'] = (
                0.8 * confidence_metrics['iptm'] + 0.2 * confidence_metrics['ptm'])
            confidence_metrics['ranking_confidence'] = (
                0.8 * confidence_metrics['iptm'] + 0.2 * confidence_metrics['ptm'])

            asym_id_set = sorted(list(set(global_var.ASYM_IDs)))
            interfaces = {}
            for i, chain1 in enumerate(asym_id_set):
                resi_i = (global_var.ASYM_IDs==chain1)
                mask_ii = resi_i[:, None] * resi_i[None, :]
                for j, chain2 in enumerate(asym_id_set[i+1:]):
                    resi_j = (global_var.ASYM_IDs==chain2)
                    mask_jj = resi_j[:, None] * resi_j[None, :]

                    resi_ij = resi_i | resi_j
                    mask_ij_full = resi_ij[:, None] * resi_ij[None, :]

                    mask_ij_interchain = mask_ij_full ^ (mask_ii | mask_jj)
                    mask_ij_interface  = ca_dist_mask * mask_ij_interchain

                    ptm_ij = predicted_tm_score(
                        logits=prediction_result['predicted_aligned_error']['logits'],
                        breaks=prediction_result['predicted_aligned_error']['breaks'],
                        pair_weights = mask_ij_full)
                    iptm_ij = predicted_tm_score(
                        logits=prediction_result['predicted_aligned_error']['logits'],
                        breaks=prediction_result['predicted_aligned_error']['breaks'],
                        pair_weights = mask_ij_interface)
                    multimer_conf_ij = 0.8*iptm_ij + 0.2*ptm_ij

                    interface_ij = {}
                    interface_ij['ptm'] = ptm_ij
                    interface_ij['iptm'] = iptm_ij
                    interface_ij['multimer_confidence'] = multimer_conf_ij
                    print(chain1, chain2, ptm_ij, iptm_ij, multimer_conf_ij)
                    interfaces[(chain1, chain2)] = interface_ij
            confidence_metrics['interfaces'] = interfaces

    # Monomer models use mean pLDDT for model ranking.
    if (num_chains == 1) and ('plddt' in confidence_metrics.keys()):
        confidence_metrics['ranking_confidence'] = np.mean(confidence_metrics['plddt'])

    return confidence_metrics


if __name__ == '__main__':

    model_metrics_path = 'O60832.Q9NY12.mmseqs.model_1_multimer_v2.metrics.backup.pkl'
    with open(model_metrics_path, 'rb') as f:
        model_metrics = pickle.load(f)

    global_var.ASYM_IDs = model_metrics['raw_asym_id']
    prediction_result = {}
    prediction_result['predicted_aligned_error'] = model_metrics['raw_pae']
    prediction_result['structure_module'] = {'final_atom_positions': model_metrics['raw_atom_positions'],
                                             'final_atom_mask': model_metrics['raw_atom_mask']}
    new_metrics = get_confidence_metrics(prediction_result, multimer_mode=True)

    score_file = [['Model_name', 'Confidence_scores', "Interfaces", 'Template', 'Refined_pdb']]
    for model_name in ['test']:
        interface_str = "|".join(["{}-{}:{:.3f}".format(
                            protein.PDB_CHAIN_IDS[int(k[0])],
                            protein.PDB_CHAIN_IDS[int(k[1])],
                            v['multimer_confidence'])
                                for k,v in new_metrics["interfaces"].items()])
        row = [model_name,
               new_metrics['multimer_confidence'],
               interface_str,
               'dummy_template',
               'dummy_refined']
        score_file.append(row)

    print(score_file)
    #with open('new_confidence_scores.csv', 'w') as f:
    #    csv_writer = csv.writer(f)
    #    csv_writer.writerows(score_file)


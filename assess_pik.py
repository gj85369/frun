#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  2 15:54:43 2026

@author: g4code
"""

import pickle
from pathlib import Path
import numpy as np
import json
import argparse
import multiprocessing as mp

parser = argparse.ArgumentParser()
parser.add_argument("--target-dir",
                    help="Directory which includes pkl structures which should be analyzed",
                    default=None,
                    type=Path)
parser.add_argument("--input-json",
                    help="INput json file to get sequences info",
                    default=None,
                    type=Path)
parser.add_argument("--np",
                    help="np",
                    default=1,
                    type=int)
cli_args = parser.parse_args()

def get_dict_from_pkl(pkl_file: Path):
    with open(pkl_file, 'rb') as f:
        data = pickle.load(f)
    return data

def get_map_chain_to_idx(input_json: Path):
    map_chain_to_idx = {}
    # load data
    with open(input_json) as f:
        data = json.load(f)
    sequences = data["sequences"]
    # get mapping
    for i, seq_dict in enumerate(sequences, 1):
        chain = seq_dict["chains"][0]
        map_chain_to_idx[chain] = i
    return map_chain_to_idx


def get_dict_of_af2_data(pkl_path: Path, map_chain_to_idx: dict):
    # get all the data
    partial_data = {}
    data = get_dict_from_pkl(pkl_path)

    # get all combinations
    all_combinations = []
    for i in range(1, len(map_chain_to_idx) + 1):
        for j in range(1, len(map_chain_to_idx) + 1):
            if i < j:
                all_combinations.append((i, j))

    # Get Antibody-antigen-data
    rec_chains = list('ABCDEFG')
    ant_chains = list('HIJKLMNOPQRSTUVWXYZ')
    aindx = []
    ridx = []
    rev_map = {}
    ptm_outs = {}
    for inst in list(map_chain_to_idx.keys()):
        if inst in rec_chains:
            ridx.append(map_chain_to_idx[inst])
            rev_map[map_chain_to_idx[inst]] = inst
        elif inst in ant_chains:
            aindx.append(map_chain_to_idx[inst])     
            rev_map[map_chain_to_idx[inst]] = inst   
            ptm_outs[inst] = {}
        else:
            print(f'some thing has gone very wrong for chain {inst}')
            quit()
        
        
    for inst in list(ptm_outs.keys()):
        #tidx = map_chain_to_idx[inst]
        ptm_outs[f'{inst}_ptms'] = []
        ptm_outs[f'{inst}_iptms'] = []
        ptm_outs[f'{inst}_multimer_confs'] = []
    
    for comb in all_combinations:
        if comb[0] in aindx:
            if comb[1] not in aindx:
                ptm_outs[f'{rev_map[comb[0]]}_ptms'].append(data["interfaces"][comb]["ptm"].item())
                ptm_outs[f'{rev_map[comb[0]]}_iptms'].append(data["interfaces"][comb]["iptm"].item())
                ptm_outs[f'{rev_map[comb[0]]}_multimer_confs'].append(data["interfaces"][comb]["multimer_confidence"].item())
    
    
    print(ptm_outs)
    # h_idx = map_chain_to_idx["H"]
    # l_idx = map_chain_to_idx["L"]
    # h_ptms = []
    # h_iptms = []
    # h_multimer_confs = []
    # l_ptms = []
    # l_iptms = []
    # l_multimer_confs = []
    # for comb in all_combinations:
    #     if (h_idx in comb) and (l_idx not in comb): # Hchain-antigen case
    #         h_ptms.append(data["interfaces"][comb]["ptm"].item())
    #         h_iptms.append(data["interfaces"][comb]["iptm"].item())
    #         h_multimer_confs.append(data["interfaces"][comb]["multimer_confidence"].item())
    #     if (l_idx in comb) and (h_idx not in comb): # Lchain-antigen case
    #         l_ptms.append(data["interfaces"][comb]["ptm"].item())
    #         l_iptms.append(data["interfaces"][comb]["iptm"].item())
    #         l_multimer_confs.append(data["interfaces"][comb]["multimer_confidence"].item())

    # # write all the data
    # partial_data["filename"] = pkl_path.stem.split(".metrics")[0]
    # partial_data["all_ptm"] = data["ptm"].item()
    # partial_data["all_iptm"] = data["iptm"].item()
    # partial_data["all_multimer_confidence"] = data["multimer_confidence"].item()
    # partial_data["H_ptm"] = np.array(h_ptms).mean()
    # partial_data["H_iptm"] = np.array(h_iptms).mean()
    # partial_data["H_multimer_confidence"] = np.array(h_multimer_confs).mean()
    # partial_data["L_ptm"] = np.array(l_ptms).mean()
    # partial_data["L_iptm"] = np.array(l_iptms).mean()
    # partial_data["L_multimer_confidence"] = np.array(l_multimer_confs).mean()
    # partial_data["H-L_ptm"] = (partial_data["H_ptm"] + partial_data["L_ptm"])/2
    # partial_data["H-L_iptm"] = (partial_data["H_iptm"] + partial_data["L_iptm"])/2
    # partial_data["H-L_multimer_confidence"] = (partial_data["H_multimer_confidence"] + partial_data["L_multimer_confidence"])/2
    # return partial_data


def mp_wrapper_get_dict_of_af2_data(args):
    return get_dict_of_af2_data(*args)


if __name__ == "__main__":
    # get mapping
    map_chain_to_idx = get_map_chain_to_idx(cli_args.input_json)
    # mp
    arglist_for_mp = [(pkl_file, map_chain_to_idx) for pkl_file in cli_args.target_dir.glob("*.pkl")] # assume that A and B are receptor antibody
    with mp.Pool(cli_args.np) as p:
        all_data_array = p.map(mp_wrapper_get_dict_of_af2_data, arglist_for_mp)
    # dump
    with open(cli_args.target_dir/"af2_info.json", "w") as f:
        json.dump(all_data_array, f, indent=4)

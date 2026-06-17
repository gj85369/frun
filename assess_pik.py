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

def get_chains(json_path: Path):
    with open(f'{json_path}/ligdic.json', 'r') as f:
        ltd = json.load(f)
    f.close()
    with open(f'{json_path}/recdic.json', 'r') as f:
        rtd = json.load(f)
    f.close()    
    print(ltd)
    print(rtd)
    
    

def get_dict_of_af2_data(pkl_path: Path, map_chain_to_idx: dict):
    # get all the data
    partial_data = {}
    data = get_dict_from_pkl(pkl_path)
    get_chains(pkl_path.parent.parent)
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
            ptm_outs[inst.lower()] = {}
        else:
            print(f'some thing has gone very wrong for chain {inst}')
            quit()
        
        
    for inst in list(ptm_outs.keys()):
        #tidx = map_chain_to_idx[inst]
        ptm_outs[inst][f'{inst}_ptms'] = []
        ptm_outs[inst][f'{inst}_iptms'] = []
        ptm_outs[inst][f'{inst}_multimer_confs'] = []
    ct = 0
    for comb in all_combinations:
        if comb[1] in aindx:
            if comb[0] not in aindx:
                ct += 1
                c0 = rev_map[comb[1]].lower()
                ptm_outs[c0][f'{c0}_ptms'].append(data["interfaces"][comb]["ptm"].item())
                ptm_outs[c0][f'{c0}_iptms'].append(data["interfaces"][comb]["iptm"].item())
                ptm_outs[c0][f'{c0}_multimer_confs'].append(data["interfaces"][comb]["multimer_confidence"].item())
    if ct == 0:
        for comb in all_combinations:
            if comb[0] in aindx:
                if comb[1] not in aindx:
                    
                    c0 = rev_map[comb[0]].lower()
                    ptm_outs[c0][f'{c0}_ptms'].append(data["interfaces"][comb]["ptm"].item())
                    ptm_outs[c0][f'{c0}_iptms'].append(data["interfaces"][comb]["iptm"].item())
                    ptm_outs[c0][f'{c0}_multimer_confs'].append(data["interfaces"][comb]["multimer_confidence"].item())
        
    
    
    partial_data["filename"] = pkl_path.stem.split(".metrics")[0]
    partial_data["all_ptm"] = data["ptm"].item()
    partial_data["all_iptm"] = data["iptm"].item()
    partial_data["all_multimer_confidence"] = data["multimer_confidence"].item()    
    
    for inst in list(ptm_outs.keys()):
        td = ptm_outs[inst]
        partial_data[f"{inst}_ptm"] = np.array(td[f'{inst}_ptms']).mean()
        partial_data[f"{inst}_iptm"] = np.array(td[f'{inst}_iptms']).mean()
        partial_data[f"{inst}_multimer_confidence"] = np.array(td[f'{inst}_multimer_confs']).mean()
    pch = list(ptm_outs.keys())
    if len(pch) > 1:
        
        for i in range(1,len(list(pch))):
            for j in range(0,i):
                partial_data[f"{pch[i].upper()}-{pch[j].upper()}_ptm"] = (partial_data[f"{pch[i]}_ptm"] + partial_data[f"{pch[j]}_ptm"])/2
                partial_data[f"{pch[i].upper()}-{pch[j].upper()}_iptm"] = (partial_data[f"{pch[i]}_iptm"] + partial_data[f"{pch[j]}_iptm"])/2
                partial_data[f"{pch[i].upper()}-{pch[j].upper()}_multimer_confidence"] = (partial_data[f"{pch[i]}_multimer_confidence"] + partial_data[f"{pch[j]}_multimer_confidence"])/2
    return partial_data
                


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

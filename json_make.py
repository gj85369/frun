#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1 09:34:58 2026

@author: g4code
"""


import copy, json

templatejson = {             
    "sequences": [ 

    ],
    "templates": [

    ],
    "parameters": "multimer_v3",
    "job_name": "",
    "num_models": 1,
    "num_recycle": 3,
    "num_predictions": 1,
    "msa_crop_size": 2048,
    "num_msa_clusters": "default",
    "seed": 0,
    "msa": "mmseqs"
}

base_template_json =         {
            "pdb_file": "",
            "mapping_method": "alignment",
            "atom_types": "all",
            "query_to_template": {
            }
        }




base_seq_entry = {
            "name": "",
            "sequence": "",
            "chains": [
                ""
            ]
        }

def making_json(inseqdic, opm, opdb, ojson, nme=None):
    bj = copy.deepcopy(templatejson)
    if nme:
        bj['job_name'] = f'{nme}'
    else:
        bj['job_name'] = 'modded'
    ttj = copy.deepcopy(base_template_json)
    ttj['pdb_file'] = opdb
    chdic = {}
    for chn, seq in inseqdic.items():
        tsj = copy.deepcopy(base_seq_entry)
        tsj['sequence'] = ''.join(seq)
        chdic[chn] = chn
        if nme:
            
            tsj['name'] = f'{opm}/{nme}_{chn}'
        else:
            tsj['name'] = f'{opm}/{chn}'
        tsj['chains'] = [chn]
        bj['sequences'].append(tsj)
    ttj['query_to_template'] = chdic
    bj['templates'].append(ttj)
    with open(ojson, 'w') as f:
        json.dump(bj, f, indent=4)
    f.close()
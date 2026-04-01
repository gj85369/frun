#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 31 22:36:44 2026

@author: g4code
"""

from protein_parser import protein, parse_line

pp = '/root/darya_test/output/v22/58B1/b1_short2_ag_ab_ACABabonly_v1_58B1_iptm_08/complex_58B1_unrelaxed_rank_001_alphafold2_multimer_v1_model_4_seed_52320497.pdb'
from pathlib import Path
ppp = Path(pp)
ok = protein(ppp)
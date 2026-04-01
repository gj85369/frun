#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 31 19:20:41 2026

@author: g4code
"""

anarci = '/workspace/colabfold/localcolabfold/v1.5.5_old_installers/localcolabfold/colabfold-conda/bin/ANARCI'
import subprocess

def check_antibody(fasta_in):
    cmd = [anarci, '-i', fasta_in]
    out = subprocess.check_output(cmd, universal_newlines=True)
    if len(out.split("\n")) > 20:
        chn = []
        for inst in out.split("\n"):
            if len(inst) > 0:
                
                chn.append(inst[0])
        return [True, max(chn)]
    else:
        return [False]
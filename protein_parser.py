#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 31 22:25:07 2026

@author: g4code
"""

d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
from itertools import product
from collections import deque


def needleman_wunsch(x, y):
    N, M = len(x), len(y)
    s = lambda a, b: int(a == b)

    DIAG = -1, -1
    LEFT = -1, 0 
    UP = 0, -1
    F = {}
    Ptr = {}

    F[-1, -1] = 0
    for i in range(N):
        F[i, -1] = -i
    for j in range(M):
        F[-1, j] = -j

    option_Ptr = DIAG, LEFT, UP
    for i, j in product(range(N), range(M)):
        option_F = (
            F[i - 1, j - 1] + s(x[i], y[j]),
            F[i - 1, j] - 1,
            F[i, j - 1] - 1,
        )
        F[i, j], Ptr[i, j] = max(zip(option_F, option_Ptr))
    alignment = deque()
    i, j = N - 1, M - 1
    while i >= 0 and j >= 0:
        direction = Ptr[i, j]
        if direction == DIAG:
            element = i, j
        elif direction == LEFT:
            element = i, None
        elif direction == UP:
            element = None, j
        alignment.appendleft(element)
        di, dj = direction
        i, j = i + di, j + dj
    while i >= 0:
        alignment.appendleft((i, None))
        i -= 1
    while j >= 0:
        alignment.appendleft((None, j))
        j -= 1

    return list(alignment)


def parse_line(inline):
    indic = {}
    if 'ATOM' in inline:
        indic['tp'] = inline[:6].strip()
        if len(inline) > 6:
            indic['atomnum'] = int(inline[6:11].strip())
            indic['atomtype'] = inline[12:16].strip()
            indic['altind'] = inline[16:17].strip()
            indic['resname'] = inline[17:20].strip()
            indic['chainname'] = inline[21:22].strip()
            indic['resnum'] = int(inline[22:26].strip())
            indic['insres'] = inline[26:27].strip()
            indic['x'] = float(inline[30:38].strip())
            indic['y'] = float(inline[38:46].strip())
            indic['z'] = float(inline[46:54].strip())
            indic['occup'] = inline[54:60].strip()
            indic['tfac'] = inline[60:66].strip()
            if len(inline) > 76:
                
                indic['elesym'] = inline[76:78].strip()
            if len(inline)> 78:
                indic['charge'] = inline[78:80].strip()
    else:
        indic['whole'] = inline
    return indic


class protein:
    def __init__(self, filename):
        self.filename = filename
        self.residues = []
        self.seq = {}
        self.outseq = {}
        self.totseq = {}
        self.file_lines = []
        self.base_lines = []        
        self.resdic = {}
        self.read_file()
        self.split_res()
        self.fas_seq = {}
        
    
    def read_file(self):
        with open(self.filename, 'r') as f:
            for line in f:
                self.base_lines.append(line)
                self.file_lines.append(parse_line(line))
        f.close()
    
    def split_res(self):
        for inst in self.file_lines:
            if 'chainname' in list(inst.keys()):
                if f"{inst['resnum']}_{inst['chainname']}_{inst['insres']}" not in list(self.resdic.keys()):
                    if inst['chainname'] not in list(self.seq.keys()):
                        self.seq[inst['chainname']] = []
                        self.outseq[inst['chainname']] = []
                        self.totseq[inst['chainname']] = []
                        self.fas_seq[inst['chainname']] = []
                        
                    self.resdic[f"{inst['resnum']}_{inst['chainname']}_{inst['insres']}"] = []
                    self.fas_seq[inst['chainname']].append(d[inst['insres']])
                self.resdic[f"{inst['resnum']}_{inst['chainname']}_{inst['insres']}"].append(inst)
        

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 31 22:59:23 2026

@author: g4code
"""

from protein_parser import protein

from ab_det import check_antibody
from protein_parser import parse_line
class complex_cut:
    def __init__(self, argsin, pdb_out_dir):
        self.argsin = argsin
        self.pdb_out_dir = pdb_out_dir
        self.fas_dic = {}
        self.receptor_chains = []
        self.ligand_chains = []
        self.nfas_dic = {}
        self.runner()
        self.lch = 'LMNO'
        self.hch = 'HIJK'
        self.l_int = 0
        self.h_int = 0

    def sel_name(self, tres):
        if res[1] == 'H':
            tnm = self.hch[self.h_int]
            self.h_int += 1
            return tnm
        if res[1] == 'L':
            tnm = self.lch[self.l_int]
            self.l_int += 1
            return tnm
    
    def parse_complex(self):
        self.prot = protein(self.argsin.complex)
        rint = 0
        rch = "ABCD"
        self.chdic = {}
        for inst in list(self.prot.fas_seq.keys()):
            self.fas_dic[inst] = ''.join(self.prot.fas_seq[inst])
            res = check_antibody(self.fas_dic[inst])
            if res[0] == True:
                self.ligand_chains.append(inst)
                self.chdic[inst] = self.sel_name(res)
                self.nfas_dic[res[1]] = ''.join(self.prot.fas_seq[inst])
            else:
                self.receptor_chains.append(inst)
                self.chdic[inst] = rch[rint]
                self.nfas_dic[rch[rint]] = ''.join(self.prot.fas_seq[inst])
                rint +=1


    def rename_chain(self, inline):
        tp = parse_line(inline)
        retline = inline[:21] + self.chdic[tp['chainname']] + inline[22:]
        return retline
            
    def make_pdbs(self):
        self.ligpdb = f'{self.pdb_out_dir}/lig.pdb'
        self.recpdb = f'{self.pdb_out_dir}/rec.pdb'
        lig_file = open(f'{self.pdb_out_dir}/lig.pdb', 'w')
        rec_file = open(f'{self.pdb_out_dir}/rec.pdb', 'w')
        for i in range(0,len(self.prot.base_lines)):
            if 'chainname' in list(self.prot.file_lines[i].keys()):
                if self.prot.file_lines[i]['chainname'] in self.receptor_chains:
                    rec_file.write(self.rename_chain(self.prot.base_lines[i]))
                else:
                    lig_file.write(self.rename_chain(self.prot.base_lines[i]))
        lig_file.close()
        rec_file.close()
    
    def runner(self):
        self.parse_complex()
        self.make_pdbs()
    
    
        
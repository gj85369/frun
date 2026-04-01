#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 31 22:59:23 2026

@author: g4code
"""

from protein_parser import protein

from ab_det import check_antibody

class complex_cut:
    def __init__(self, argsin, pdb_out_dir):
        self.argsin = argsin
        self.pdb_out_dir = pdb_out_dir
        self.fas_dic = {}
        self.receptor_chains = []
        self.ligand_chains = []
        self.runner()
    
    def parse_complex(self):
        self.prot = protein(self.argsin.complex)
        for inst in list(self.prot.fas_seq.keys()):
            self.fas_dic[inst] = ''.join(self.prot.fas_seq[inst])
            res = check_antibody(self.fas_dic[inst])
            if res[0] == True:
                self.ligand_chains.append(inst)
            else:
                self.receptor_chains.append(inst)
        
            
    def make_pdbs(self):
        lig_file = open(f'{self.pdb_out_dir}/lig.pdb', 'w')
        rec_file = open(f'{self.pdb_out_dir}/rec.pdb', 'w')
        for i in range(0,len(self.prot.base_lines)):
            if 'chainname' in list(self.prot.file_lines[i].keys()):
                if self.prot.file_lines[i]['chainname'] in self.receptor_chains:
                    rec_file.write(self.prot.base_lines[i])
                else:
                    lig_file.write(self.prot.base_lines[i])
        lig_file.close()
        rec_file.close()
    
    def runner(self):
        self.parse_complex()
        self.make_pdbs()
    
    
        
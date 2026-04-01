#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1 09:55:53 2026

@author: g4code
"""


from sequence_similarity import needleman_wunsch
from utils import read_fasta

class compare:
    def __init__(self, tincomplex, tlig_dic, trec_dic):
        self.tincomplex = tincomplex
        self.tlig_dic = tlig_dic
        self.trec_dic = trec_dic
        self.compare_fastas()
    def parse_ned_out(self, ndout):
        ct = 0
        for inst in ndout:
            if None in inst:
                ct +=1
        return ct


        
    def compare_fastas(self):
        compo = {}
        ccomp = {}
        for inst in list(self.tincomplex.nfas_dic.keys()):
            ccomp[inst] = []
            for sinst in list(self.tlig_dic.keys()):
                compo[f'{inst}_l{sinst}'] = self.parse_ned_out(needleman_wunsch(self.tincomplex.nfas_dic[inst], read_fasta(self.tlig_dic[sinst]['nfasta'])))
                ccomp[inst].append(compo[f'{inst}_l{sinst}'])
            for sinst in list(self.trec_dic.keys()):
                compo[f'{inst}_r{sinst}'] = self.parse_ned_out(needleman_wunsch(self.tincomplex.nfas_dic[inst], read_fasta(self.trec_dic[sinst]['nfasta'])))            
                ccomp[inst].append(compo[f'{inst}_r{sinst}'])
        dones = []
        for inst in list(compo.keys()):
            cn, fn = inst.split('_')
            if inst not in dones:
                if compo[inst] == min(ccomp[cn]):
                    if fn[0] == 'r':
                        self.trec_dic[int(fn[1])]['chain_name'] = cn
                        dones.append(inst)
                    if fn[0] == 'l':
                        self.tlig_dic[int(fn[1])]['chain_name'] = cn 
                        dones.append(inst)
        print(compo)    
        print(ccomp)
        print(self.tlig_dic)
        print(self.trec_dic)
        
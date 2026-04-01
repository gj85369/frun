#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1 09:30:35 2026

@author: g4code
"""

from protein_parser import protein
from utils import write_rules, write_msa, mod_mismatch
import os


class interact:
    def __init__(self, rec_file, lig_file, outfile, msa_indir, msa_outdir, name=None, outname=None, intrange=16):
        self.lig_obj = protein(lig_file)
        self.rec_obj = protein(rec_file)
        self.intrange = intrange
        self.name = name
        self.outname = outname
        self.lig_dic = self.lig_obj.resdic
        self.rec_dic = self.rec_obj.resdic
        self.lig_lines = self.lig_obj.base_lines
        self.rec_lines = self.rec_obj.base_lines
        self.lig_inseq = self.lig_obj.seq
        self.lig_outseq = self.lig_obj.outseq
        self.rec_inseq = self.rec_obj.seq
        self.rec_outseq = self.rec_obj.outseq    
        self.lig_tot_seq = self.lig_obj.totseq
        self.rec_tot_seq = self.rec_obj.totseq
        self.msa_indir = msa_indir
        self.msa_outdir = msa_outdir
        self.outfile = outfile        
        self.lig_res_int = []
        self.rec_res_int = []
        self.calc_int()
        self.write_out()  #uncomment this you shit
        self.all_seqs = {}
        self.all_out_seqs = {}
        self.make_seq_dic()
        self.msa_dic = {}
        self.get_msas()
        
    def calc_int(self):
        for rinst in list(self.rec_dic.keys()):
            for linst in list(self.lig_dic.keys()):
                if self.calc_each(self.lig_dic[linst], self.rec_dic[rinst]):
                    self.rec_res_int.append(rinst)
                    self.lig_res_int.append(linst)
        self.rec_res_int = list(set(self.rec_res_int))
        self.lig_res_int = list(set(self.lig_res_int))
    
    def calc_each(self, lig_list, rec_list):
        for lig_inst in lig_list:
            for rec_inst in rec_list:
                tx = (lig_inst['x'] - rec_inst['x'])**2
                ty = (lig_inst['y'] - rec_inst['y'])**2            
                tz = (lig_inst['z'] - rec_inst['z'])**2
                if tx + ty + tz <= self.intrange:
                    return True
        return False
    def write_out(self):
        with open(self.outfile, 'w') as f:
            write_rules(self.lig_lines, self.lig_res_int, f, list(self.lig_dic.keys()), self.lig_inseq, self.lig_outseq, self.lig_tot_seq)
            f.write("TER\n")            
            write_rules(self.rec_lines, self.rec_res_int, f, list(self.rec_dic.keys()), self.rec_inseq, self.rec_outseq, self.rec_tot_seq)
            
        f.close()
     
    def make_seq_dic(self):
        for inst in list(self.lig_tot_seq.keys()):
            self.all_seqs[inst] = self.lig_tot_seq[inst]
            self.all_out_seqs[inst] = self.lig_outseq[inst]
        for inst in list(self.rec_tot_seq.keys()):
            self.all_seqs[inst] = self.rec_tot_seq[inst] 
            self.all_out_seqs[inst] = self.rec_outseq[inst]
     
    
    def get_msas(self):
        for inst in list(self.all_seqs.keys()):
            tlist = []
            if self.name == None:
                afile = f'{self.msa_indir}/{inst}/mmseqs/aggregated.a3m'
            else:
                afile = f'{self.msa_indir}/{self.name}_{inst}/mmseqs/aggregated.a3m'
            with open(afile, 'r') as f:
                for line in f:
                    tlist.append(line)
            f.close()
            self.msa_dic[inst] = {}
            self.msa_dic[inst]['raw'] = tlist
            msa_seq = tlist[1].strip()
            self.msa_dic[inst]['msa_seq'] = msa_seq

                
            self.msa_dic[inst]['mod'] = mod_mismatch(tlist, ''.join(self.all_seqs[inst]), self.all_out_seqs[inst])
            if self.outname == None:
                os.makedirs(f'{self.msa_outdir}/{inst}/mmseqs', exist_ok=True)
                ofile = f'{self.msa_outdir}/{inst}/mmseqs/aggregated.a3m'
            else:
                os.makedirs(f'{self.msa_outdir}/{self.outname}_{inst}/mmseqs', exist_ok=True)                
                ofile = f'{self.msa_outdir}/{self.outname}_{inst}/mmseqs/aggregated.a3m'
            write_msa(ofile, self.msa_dic[inst]['mod'])                

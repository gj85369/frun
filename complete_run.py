#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 31 10:32:55 2026

@author: g4code
"""
import argparse
import os
import tempfile
from pathlib import Path
from complex_cut import complex_cut
from sequence_similarity import needleman_wunsch
from ab_det import check_antibody
from comparing_fastas import compare

def main(argsin):
    # with tempfile.TemporaryDirectory() as tdname:
        
    #     try:
    #         singrun = prepare_runner(tdname, argsin)
    #         singrun.runner()
    #     except Exception as e:
    #         print(f'Error {e}')
    tdname = f'{os.getcwd()}/funny'
    os.makedirs(tdname, exist_ok=True)
    try:
        singrun = prepare_runner(tdname, argsin)
        singrun.runner()
    except Exception as e:
        print(f'Error {e}')       
    
import subprocess


    
class cleanup:
    def __init__(self, *args, **kwargs):
        pass
    

class prepare_runner:
    def __init__(self, workdir, argsp):
        self.workdir = workdir
        self.argsp = argsp
        
    def make_ind_fasta(self, a3m_path, oname, faspas):
        cmd = f'echo ">{oname}" > {faspas}'
        subprocess.run(cmd, shell=True)
        cmd = f'head -n 2 {a3m_path} | tail -n 1 >> {faspas}'
        subprocess.run(cmd, shell=True)
        cmd = f'cat {faspas}'
        subprocess.run(cmd, shell=True)
        cmd = f'ls -ltr {faspas}'
        subprocess.run(cmd, shell=True)        
    def make_fastas(self):
        
        ###   THIS NEEDS TO CHANGE - run parse complex first then assign chain names after, then move msas
        #self.make_ind_fasta(self.argsp.receptor, 'rec', f'{self.workdir}/fasta/rec.fasta')
        #self.rec_fasta = f'{self.workdir}/fasta/rec.fasta'
        self.rec_dic = {}
        if len(self.argsp.receptor) > 4:
            print('receptor limit hard coded to 4 search this text to change')
            quit()
        rch = 'ABCD'

        for i in range(0, len(self.argsp.receptor)):
            
            self.make_ind_fasta(self.argsp.receptor[i], f'rec_{i}', f'{self.workdir}/fasta/rec_{i}.fasta')
            self.rec_dic[i] = {}
            self.rec_dic[i]['nfasta'] = f'{self.workdir}/fasta/rec_{i}.fasta'
            self.rec_dic[i]['ofasta'] = self.argsp.receptor[i]    
            self.rec_dic[i]['chain_name'] = None
            
        self.lig_dic = {}
        for i in range(0, len(self.argsp.ligand)):
            
            self.make_ind_fasta(self.argsp.ligand[i], f'lig_{i}', f'{self.workdir}/fasta/lig_{i}.fasta')
            self.lig_dic[i] = {}
            self.lig_dic[i]['nfasta'] = f'{self.workdir}/fasta/lig_{i}.fasta'
            self.lig_dic[i]['ofasta'] = self.argsp.ligand[i]
            opt = check_antibody(self.lig_dic[i]['nfasta'])
            if opt[0] == False:
                print(f'ligand is supposed to be antibody {self.argsp.ligand[i]} does not match antibody')
                quit()
            else:
                print(f'{self.argsp.ligand[i]} is an antibody chain {opt[1]}')
                self.lig_dic[i]['chain_name'] = None
            
            print(self.lig_dic)
            print(self.rec_dic)
            
        
    def parse_complex(self):
        self.incomplex = complex_cut(self.argsp, f'{self.workdir}/pdbs')
        

    def compare_fastas(self):
        tcomp = compare(self.incomplex, self.lig_dic, self.rec_dic)

        
    def runner(self):
        os.makedirs(self.argsp.output_dir, exist_ok=True)
        os.makedirs(f'{self.workdir}/msa', exist_ok=True)
        os.makedirs(f'{self.workdir}/ligand', exist_ok=True)
        os.makedirs(f'{self.workdir}/fasta', exist_ok=True)
        os.makedirs(f'{self.workdir}/pdbs', exist_ok=True)
        
        self.make_fastas()
        print(os.getcwd())
        self.parse_complex()
        self.compare_fastas()
        
    



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--nanobody',help='is the system a nanobody',  action='store_true')    
    parser.add_argument("-r",
                        "--receptor",
                        help="Receptor msa",
                        required=True,  
                        action='append',                        
                        type=Path)
    parser.add_argument("-l",
                        "--ligand",
                        help="ligand msa",
                        action='append',
                        type=Path)    
    parser.add_argument("-c",
                        "--complex",
                        help="The input complex",
                        required=True,
                        type=Path)
    parser.add_argument("-o",
                        "--output-dir",
                        help="dir for work and output.",
                        required=True,
                        type=Path)
    args = parser.parse_args()
    if not args.ligand:
        print('ligand msas needed')
        quit()
    if args.nanobody:
        if len(args.ligand) > 1:
            print('only 1 input msa for nanobodies')
            quit()

    else:
        if len(args.ligand) != 2:
            print('2 msas needed for antibodies')
            quit()        
    
            
        
    main(args)
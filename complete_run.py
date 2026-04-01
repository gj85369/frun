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
from protein_parser import protein, parse_line
anarci = '/workspace/colabfold/localcolabfold/v1.5.5_old_installers/localcolabfold/colabfold-conda/bin/ANARCI'

def main(argsin):
    with tempfile.TemporaryDirectory() as tdname:
        
        try:
            singrun = prepare_runner(tdname, argsin)
            singrun.runner()
        except Exception as e:
            print(f'Error {e}')
        
    
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
            self.rec_dic[i]['chain_name'] = rch[i]
            
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
                self.lig_dic[i]['chain_name'] = opt[1]
            
            print(self.lig_dic)
            print(self.rec_dic)
            
        
    def parse_complex(self):
        incomplex = protein(self.argsp.complex)
        print(incomplex.fas_seq)
            
        
    def runner(self):
        os.makedirs(f'{self.workdir}/msa')
        os.makedirs(f'{self.workdir}/ligand')
        os.makedirs(f'{self.workdir}/fasta')
        self.make_fastas()
        print(os.getcwd())
        
    



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
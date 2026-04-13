#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 31 10:32:55 2026

@author: g4code
"""
import argparse
import os
import tempfile, subprocess
from pathlib import Path
from complex_cut import complex_cut
#from sequence_similarity import needleman_wunsch
from ab_det import check_antibody
from comparing_fastas import compare
import shutil
from interaction import interact
from json_make import making_json
from cleaning import cleaning
from glob import glob
from multiprocessing import Pool
import copy
bpath = Path(__file__).parent.resolve()



def run_inst(instance):
    instance.runner()
    return instance

def load_stuff():
    pass


def main(argsin):
    cfiles = glob(f'{argsin.complex_dir}/*.pdb')
    instances = []
    for i in range(0,len(cfiles)):
        finst = cfiles[i]
        
        argsin.complex = finst
        argsin.output_dir = f'{argsin.output_dir}/{i}'
        a1 = copy.deepcopy(argsin)
        instances.append(prepare_runner(a1))
    with Pool(processes=6) as p:
        instances = p.map(run_inst, instances)
    
    


    
class cleanup:
    def __init__(self, *args, **kwargs):
        pass
    

class prepare_runner:
    def __init__(self,  argsp):
        self.argsp = argsp
        print(f'running {self.argsp.complex}')
        
    def make_ind_fasta(self, a3m_path, oname, faspas):
        cmd = f'echo ">{oname}" > {faspas}'
        subprocess.run(cmd, shell=True)
        cmd = f'head -n 2 {a3m_path} | tail -n 1 >> {faspas}'
        subprocess.run(cmd, shell=True)
        #cmd = f'cat {faspas}'
        #subprocess.run(cmd, shell=True)
        #cmd = f'ls -ltr {faspas}'
        #subprocess.run(cmd, shell=True)        
    def make_fastas(self):
        
        ###   THIS NEEDS TO CHANGE - run parse complex first then assign chain names after, then move msas
        #self.make_ind_fasta(self.argsp.receptor, 'rec', f'{self.workdir}/fasta/rec.fasta')
        #self.rec_fasta = f'{self.workdir}/fasta/rec.fasta'
        self.rec_dic = {}
        if len(self.argsp.receptor) > 4:
            print('receptor limit hard coded to 4 search this text to change')
            quit()

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
                #print(f'{self.argsp.ligand[i]} is an antibody chain {opt[1]}')
                self.lig_dic[i]['chain_name'] = None
            

            
        
    def parse_complex(self):
        self.incomplex = complex_cut(self.argsp, f'{self.workdir}/pdbs')
        

    def compare_fastas(self):
        tcomp = compare(self.incomplex, self.lig_dic, self.rec_dic)
        self.lig_dic = tcomp.tlig_dic
        self.rec_dic = tcomp.trec_dic
    
    
    def make_new_msas(self):
        for inst in list(self.lig_dic.keys()):
            os.makedirs(f'{self.workdir}/msa/new_{self.lig_dic[inst]["chain_name"]}/mmseqs', exist_ok=True)
            shutil.copyfile(self.lig_dic[inst]['ofasta'], f'{self.workdir}/msa/new_{self.lig_dic[inst]["chain_name"]}/mmseqs/aggregated.a3m')
        for inst in list(self.rec_dic.keys()):
            os.makedirs(f'{self.workdir}/msa/new_{self.rec_dic[inst]["chain_name"]}/mmseqs', exist_ok=True)
            shutil.copyfile(self.rec_dic[inst]['ofasta'], f'{self.workdir}/msa/new_{self.rec_dic[inst]["chain_name"]}/mmseqs/aggregated.a3m')
                        

    def run_interaction(self):
        os.makedirs(f'{self.workdir}/newrun')
        os.makedirs(f'{self.workdir}/newrun/processed_templates')
        os.makedirs(f'{self.workdir}/newrun/processed_msas')
        os.makedirs(f'{self.workdir}/newrun/af_run')
        
        intout = interact(self.incomplex.recpdb, 
                          self.incomplex.ligpdb, 
                          f'{self.workdir}/newrun/processed_templates/temp.pdb', 
                          f'{self.workdir}/msa', 
                          f'{self.workdir}/newrun/processed_msas',
                          name='new', 
                          outname='modded')
        seq_dic = {}
        for inst in list(intout.msa_dic.keys()):
            seq_dic[inst] = intout.msa_dic[inst]['msa_seq']
        making_json(seq_dic, 
                    f'{self.workdir}/newrun/processed_msas', 
                    f'{self.workdir}/newrun/processed_templates/temp.pdb', 
                    f'{self.workdir}/newrun/af_run/job.json', 
                    nme='modded')
        
    def run_af(self):
        os.chdir(f'{self.workdir}/newrun')
        cmd = f'bash {bpath}/runit'
        subprocess.check_call(cmd.split())
        
    def run_cleaning(self):
        cln = cleaning(self.workdir, self.argsp)
        
    def runner(self):
        with tempfile.TemporaryDirectory() as tdname:
            self.workdir = tdname

            os.makedirs(self.argsp.output_dir, exist_ok=True)
            os.makedirs(f'{self.workdir}/msa', exist_ok=True)
            os.makedirs(f'{self.workdir}/ligand', exist_ok=True)
            os.makedirs(f'{self.workdir}/fasta', exist_ok=True)
            os.makedirs(f'{self.workdir}/pdbs', exist_ok=True)
            self.make_fastas()
            self.parse_complex()
            self.compare_fastas()
            self.make_new_msas()
            self.run_interaction()
            self.run_af()
            self.run_cleaning()
    



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
                        "--complex-dir",
                        help="The input complex",
                        required=True,
                        type=Path)
    parser.add_argument("-s",
                        "--complex",
                        help="The input complex",
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
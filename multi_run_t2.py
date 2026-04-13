#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 31 10:32:55 2026

@author: g4code
"""
import argparse
import sys
import os
import tempfile, subprocess
from pathlib import Path
#from sequence_similarity import needleman_wunsch

from glob import glob
from multiprocessing import Pool
import copy
bpath = Path(__file__).parent.resolve()



def run_inst(cmd):
    print(cmd)
    outs = subprocess.check_output(cmd, shell=True, universal_newlines=True)
    return outs



def main(argsin):
    cfiles = glob(f'{argsin.complex_dir}/*.pdb')
    pyinst = [sys.executable]
    vlist = copy.deepcopy(sys.argv)
    vlist[0] = f'{bpath}/complete_run.py'
    cind = vlist.index('-c')
    oind = vlist.index('-o')

    instances = []
    for i in range(0,len(cfiles)):
        finst = cfiles[i]
        vtmp = copy.deepcopy(vlist)
        vtmp[cind + 1] = finst
        vtmp[oind + 1] = f'{vlist[oind+1]}/{i}'
        tl = pyinst + vtmp
        instances.append(' '.join(tl))
        
        
    with Pool(processes=6) as p:
        instances = p.map(run_inst, instances)
    
    






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
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

def main(argsin):
    with tempfile.TemporaryDirectory() as tdname:
        
        try:
            singrun = prepare_runner(tdname)
            singrun.runner()
        except Exception as e:
            print(f'Error {e}')
        
    



class cleanup:
    def __init__(self, *args, **kwargs):
        pass
    

class prepare_runner:
    def __init__(self, workdir):
        self.workdir = workdir
        
    def runner(self):
        print(os.getcwd())
        
    



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--nanobody',help='is the system a nanobody',  action='store_true')    
    parser.add_argument("-r",
                        "--receptor",
                        help="Receptor msa",
                        required=True,                        
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
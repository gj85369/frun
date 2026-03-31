#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 31 10:32:55 2026

@author: g4code
"""
import argparse

from pathlib import Path

def main(argsin):
    pass



class cleanup:
    def __init__(self, *args, **kwargs):
        pass
    

class prepare_runner:
    def __init__(self, *args, **kwargs):
        pass
    



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--nanobody',help='is the system a nanobody',  action='store_true')    
    parser.add_argument("-r",
                        "--receptor",
                        help="Receptor msa",
                        required=True,                        
                        type=Path)
    parser.add_argument("-h",
                        "--heavy",
                        help="Heavy chain msa",
                        type=Path)  
    parser.add_argument("-l",
                        "--light",
                        help="Light chain msa",
                        type=Path)    
    parser.add_argument("-s",
                        "--single",
                        help="Nanobody msa",
                        type=Path)    
    parser.add_argument("-m",
                        "--msa-dir",
                        help="MSA dir that stores msas.",
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
    if args.nanobody:
        iss = 0
        if args.light:
            print('please only use the single option for msa input')
            iss += 1
        if args.heavy:
            print('please only use the single option for msa input')    
            iss += 1
        
        if not args.single:
            print('Single msa require for nanobody')    
            iss += 1
        
        if iss > 0:
            quit()
        
            
        
    main(args)
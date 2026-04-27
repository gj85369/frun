#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1 11:14:43 2026

@author: g4code
"""
from glob import glob
import shutil
import tarfile
from pathlib import Path

class cleaning:
    def __init__(self, workdir, argsin):
        self.argsin = argsin
        self.workdir = workdir
        self.copy_files()
        
    def copy_files(self):
        afrun_files = glob(f'{self.workdir}/newrun/af_run/*')
        files = [x for x in afrun_files if Path.is_file(Path(x))]
        
        for inst in files:
            fname = inst.split('/')[-1]
            if fname not in ['FUN']:
                shutil.copy(inst, f'{self.argsin.output_dir.absolute()}/{fname}')
        
        # a3ms = glob(f'{self.workdir}/newrun/processed_msas/*/*/*.a3m')
        # for inst in a3ms:
        #     alist = inst.split('/')
            
        #     shutil.copy(inst, f'{self.argsin.output_dir}/{alist[-3]}_{alist[-1]}')
            
        

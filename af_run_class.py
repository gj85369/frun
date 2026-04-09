#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr  5 13:58:29 2026

@author: g4code
"""


import subprocess
import json
import os


class af_runner:
    def __init__(self, in_json, rundir):
        self.in_json = in_json
        self.log = {}
        self.rundir = rundir
        
        
        
    def run_it(self):
        cmd = f'{self.in_json["python"]} {self.in_json["af2_script"]} job.json'
        ok = subprocess.check_output(cmd, universal_newlines=True, shell=True, cwd=self.rundir)
        self.log['output'] = ok

        cmd = f'{self.in_json["python"]} {self.in_json["af2_script"]} job.json'
        ok = subprocess.check_output(cmd, universal_newlines=True, shell=True, cwd=self.rundir)
        self.log['output'] = ok        
        
        cmd = f'{self.in_json["python"]} {self.in_json["af2_script"]} job.json'
        ok = subprocess.check_output(cmd, universal_newlines=True, shell=True, cwd=self.rundir)
        self.log['output'] = ok        
    
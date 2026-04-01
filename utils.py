#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1 09:31:58 2026

@author: g4code
"""


from protein_parser import parse_line, d
from sequence_similarity import needleman_wunsch

def write_rules(inlist, exlist, fobj, keylist, inseqdic, outseqdic, totseqdic):
    curchain = None
    startnum = {}
    chains = list(set([x.split('_')[1] for x in keylist]))
    dones = []
    for inst in chains:
        startnum[inst] = {}
        startnum[inst]['min'] = min([int(x.split('_')[0]) for x in keylist if x.split('_')[1] == inst])
        startnum[inst]['offset'] = 0
    for cur_line in inlist:
        if 'ATOM' in cur_line:
            rtmp = parse_line(cur_line)
            if curchain == None:
                curchain = rtmp['chainname']
                
            else:
                if curchain != rtmp['chainname']:
                    curchain = rtmp['chainname']
                    fobj.write("TER\n")
                    
            if f"{rtmp['resnum']}_{rtmp['chainname']}_{rtmp['insres']}" not in exlist:
                sn = startnum[rtmp['chainname']]['min']
                if f"{rtmp['resnum']}_{rtmp['chainname']}_{rtmp['insres']}"  not in dones:
                    
                    inseqdic[curchain].append(d[rtmp['resname']])
                    outseqdic[curchain].append('-')
                    totseqdic[curchain].append(d[rtmp['resname']])                    
                    dones.append(f"{rtmp['resnum']}_{rtmp['chainname']}_{rtmp['insres']}")
                
                if  sn == 1:
                    
                    fobj.write(cur_line)
                else:
                    cur_line_b = cur_line[:22]
                    cur_line_a = cur_line[26:]
                    numdiff = rtmp['resnum'] - sn + 1
                    modnum = '{:4s}'.format(str(numdiff))
                    tline = cur_line_b + modnum + cur_line_a
                    fobj.write(tline)
            else:
                if f"{rtmp['resnum']}_{rtmp['chainname']}_{rtmp['insres']}" not in dones:
                    
                    inseqdic[curchain].append('-')
                    outseqdic[curchain].append(d[rtmp['resname']])
                    totseqdic[curchain].append(d[rtmp['resname']])
                    dones.append(f"{rtmp['resnum']}_{rtmp['chainname']}_{rtmp['insres']}")
                
def getnewseq(msaseq, pdbseq, tempseq):
    fstmod = []
    nws = needleman_wunsch(msaseq.strip(), pdbseq)
    #msl = list(msaseq.strip())
    tsl = list(tempseq)
    for i in range(0,len(nws)):
        if nws[i][1] == None:
            fstmod.append('-')
        elif nws[i][0] == None:
            pass
        else:
            fstmod.append(tsl[nws[i][1]])
    return ''.join(fstmod)



def mod_mismatch(in_msa, pdbseq, exl_list):
    ret_list = []
    ret_list.append(in_msa[0])
    ret_list.append(in_msa[1])
    exlseq = ''.join(exl_list)
    nseq = list(getnewseq(in_msa[1], pdbseq, exlseq))
    for i in range(2,len(in_msa)):
        if '>' in in_msa[i]:
            ret_list.append(in_msa[i])
        else:
            nlist = []
            pre = list(in_msa[i].strip())
            msalinelist = [x for x in pre if not x.islower()]
            for j in range(0,len(nseq)):
                if nseq[j] == '-':
                    nlist.append('-')
                else:
                    nlist.append(msalinelist[j])
            nline = ''.join(nlist)
            ret_list.append(f"{nline}\n")
    return ret_list    
    


    

def write_msa(msa_name, msa_list):
    with open(msa_name, 'w') as f:
        for inst in msa_list:
            f.write(inst)
    f.close()
    
def parse_msains(ininst, names, chains):
    if len(ininst.split('_')) == 2:
        tnm, tch = ininst.split('_')
        names.append(tnm)
        chains.append(tch)
    elif len(ininst.split('_')) == 1:
        chains.append(ininst)
    else:
        raise ValueError(f'please change names so there isnt multiple underscores {ininst}')
        
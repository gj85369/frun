# The code below is based on AlphaFold Multimer Colab implementation.

# Copyright 2021 DeepMind Technologies Limited
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import os
import sys
import argparse
from glob import glob
from urllib import request
from concurrent import futures

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO

af2_path = os.getenv('ALPHAFOLD_MULTIMER_PATH')
#AF2_PARAMS_PATH = os.getenv('ALPHAFOLD_PARAMS_PATH')
#Add alphafold to PATH
PATHs = sys.path
if af2_path not in PATHs:
    sys.path.append(af2_path)

from alphafold.data import pipeline
from alphafold.data import parsers
from alphafold.data.tools import jackhmmer
from alphafold.data.tools import hhblits
from alphafold.data.parsers import Msa
from alphafold.notebooks import notebook_utils

def get_api_url():
    test_url_pattern = 'https://storage.googleapis.com/alphafold-colab{:s}/latest/uniref90_2021_03.fasta.1'
    ex = futures.ThreadPoolExecutor(3)
    def fetch(source):
        request.urlretrieve(test_url_pattern.format(source))
        return source
    fs = [ex.submit(fetch, source) for source in ['', '-europe', '-asia']]
    source = None
    for f in futures.as_completed(fs):
        source = f.result()
        ex.shutdown()
        break
    api_url =  f'https://storage.googleapis.com/alphafold-colab{source}/latest/'
    return api_url

#DB_ROOT_PATH =  f'https://storage.googleapis.com/alphafold-colab/latest/'

def init_db_configs_api(db_root_path, db_list):
    raw_db_dict = {
        'uniref90': {
            'db_name': 'uniref90',
            'runner': 'hmmer',
            'db_path': f'{db_root_path}uniref90_2021_03.fasta',
            'num_streamed_chunks': 59,
            #'num_streamed_chunks': 1,
            'z_value': 135_301_051,
            'max_hits': 10_000,
            'pairable': False},
        'small_bfd': {
            'db_name': 'small_bfd',
            'runner': 'hmmer',
            'db_path': f'{db_root_path}bfd-first_non_consensus_sequences.fasta',
            'num_streamed_chunks': 17,
            #'num_streamed_chunks': 1,
            'z_value': 65_984_053,
            'max_hits': 5_000,
            'pairable': False},
        'mgnify': {
            'db_name': 'mgnify',
            'runner': 'hmmer',
            'db_path': f'{db_root_path}mgy_clusters_2019_05.fasta',
            'num_streamed_chunks': 71,
            #'num_streamed_chunks': 1,
            'z_value': 304_820_129,
            'max_hits': 501,
            'pairable': False},
        'uniprot': {
            'db_name': 'uniprot',
            'runner': 'hmmer',
            'db_path': f'{db_root_path}uniprot_2021_03.fasta',
            'num_streamed_chunks': 98,
            #'num_streamed_chunks': 1,
            'z_value': 219_174_961 + 565_254,
            'max_hits': 50_000,
            'pairable': True},
    }
    db_configs = []
    for db in db_list:
        db_configs.append(raw_db_dict[db])
    return db_configs

def init_db_configs(db_root_path, db_list):
    raw_db_dict = {
        'uniref90': {
            'db_name': 'uniref90',
            'runner': 'hmmer',
            'db_path': os.path.join(db_root_path, 'uniref90', 'uniref90.fasta'),
            'num_streamed_chunks': None,
            'z_value': None,
            'max_hits': 10000,
            'pairable': False},
        'mgnify': {
            'db_name': 'mgnify',
            'runner': 'hmmer',
            'db_path': os.path.join(db_root_path, 'mgnify', 'mgy_clusters_2018_12.fa'),
            'num_streamed_chunks': None,
            'z_value': None,
            'max_hits': 501,
            'pairable': False},
        'small_bfd': {
            'db_name': 'small_bfd',
            'runner': 'hmmer',
            'db_path': os.path.join(db_root_path, 'small_bfd', 'bfd-first_non_consensus_sequences.fasta'),
            'num_streamed_chunks': None,
            'z_value': None,
            'max_hits': None,
            'pairable': False},
        'bfd': {
            'db_name': 'bfd',
            'runner': 'hhblits',
            'db_path': os.path.join(db_root_path, 'bfd', 'bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt'),
            'num_streamed_chunks': None,
            'z_value': None,
            'max_hits': None,
            'pairable': False},
        'uniclust30': {
            'db_name': 'uniclust30',
            'runner': 'hhblits',
            'db_path': os.path.join(db_root_path, 'uniclust30', 'uniclust30_2018_08', 'uniclust30_2018_08'),
            'num_streamed_chunks': None,
            'z_value': None,
            'max_hits': None,
            'pairable': False},
        'bfd_uniclust': {
            'db_name': 'bfd+uniclust30',
            'runner': 'hhblits',
            'db_path': [os.path.join(db_root_path, 'bfd', 'bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt'),
                        os.path.join(db_root_path, 'uniclust30', 'uniclust30_2018_08', 'uniclust30_2018_08')],
            'num_streamed_chunks': None,
            'z_value': None,
            'max_hits': None,
            'pairable': False},
        'uniprot': {
            'db_name': 'uniprot',
            'runner': 'hmmer',
            'db_path': os.path.join(db_root_path, 'uniprot', 'uniprot.fasta'),
            'num_streamed_chunks': None,
            'z_value': None,
            'max_hits': None,
            'pairable': True}
    }
    db_configs = []
    for db in db_list:
        db_configs.append(raw_db_dict[db])
    return db_configs


def get_MSAs_for_seq_api(seq_str,
                         output_dir,
                         db_configs,
                         hmmer_binary_path):
    fasta_path = os.path.join(output_dir, 'query.fasta')
    with open(fasta_path, 'wt') as f:
        f.write(f'>query\n{seq_str}')

    unpairable_msas = []
    pairable_msas = []
    for db_config in db_configs:
        db_name = db_config['db_name']
        pairable = db_config['pairable']
        jackhmmer_runner = jackhmmer.Jackhmmer(
            binary_path=hmmer_binary_path,
            database_path=db_config['db_path'],
            get_tblout=True,
            num_streamed_chunks=db_config['num_streamed_chunks'],
            z_value=db_config['z_value'])
        raw_msa_result = jackhmmer_runner.query(fasta_path)
        merged_msa = notebook_utils.merge_chunked_msa(
                            results=raw_msa_result,
                            max_hits=db_config['max_hits'])
        if merged_msa.sequences and not pairable:
            unpairable_msas.append(merged_msa)
        elif merged_msa.sequences and pairable:
            pairable_msas.append(merged_msa)
        msa_size = len(set(merged_msa.sequences))
        print(f'{msa_size} unique sequences found in {db_name}.')

    print('Msas: {} unpairable and {} pairable'.format(len(unpairable_msas), len(pairable_msas)))
    return {'unpairable':unpairable_msas, 'pairable':pairable_msas}


def get_MSAs_for_seq_std(seq_str,
                         output_dir,
                         db_configs,
                         hmmer_binary_path,
                         hhblits_binary_path,
                         small_bfd=False):
    fasta_path = os.path.join(output_dir, 'query.fasta')
    with open(fasta_path, 'wt') as f:
        f.write(f'>query\n{seq_str}')

    runner_configs = {
        'hmmer': {
            'runner': jackhmmer.Jackhmmer,
            'binary': hmmer_binary_path,
            'db_arg': 'database_path',
            'format': 'sto',
            'parser': parsers.parse_stockholm},
        'hhblits': {
            'runner': hhblits.HHBlits,
            'binary': hhblits_binary_path,
            'db_arg': 'databases',
            'format': 'a3m',
            'parser': parsers.parse_a3m},
    }

    unpairable_msas = []
    pairable_msas = []
    for db_config in db_configs:
        db_name = db_config['db_name']
        runner_type = db_config['runner']
        runner_config = runner_configs[runner_type]
        search_args = {'binary_path': runner_config['binary'],
                       runner_config['db_arg']: db_config['db_path']}
        search_runner = runner_config['runner'](**search_args)

        msa_path = os.path.join(output_dir, db_name+'.'+runner_config['format'])

        raw_msa_result = pipeline.run_msa_tool(
            msa_runner=search_runner,
            input_fasta_path=fasta_path,
            msa_out_path=msa_path,
            msa_format=runner_config['format'],#a3m
            max_sto_sequences=db_config['max_hits'], #only for jackhmmer dbs
            use_precomputed_msas=True)

        raw_msa = raw_msa_result[runner_config['format']]
        parsed_msa = runner_config['parser'](raw_msa)

        print('MSA depth: {} - {} sequences'.format(db_name, len(parsed_msa)))

        if db_config['pairable'] == False:
            unpairable_msas.append(parsed_msa)
        else:
            pairable_msas.append(parsed_msa)
    return {'unpairable': unpairable_msas, 'pairable': pairable_msas}


def get_MSAs_googleapi(seq_dict,
                       output_path = './',
                       hmmer_binary_path='/usr/bin/jackhmmer',
                       pairing=False,
                       debug=False):
    if not os.path.isdir('/tmp/ramdisk'):
        print("Expected to find a ramdisk at /tmp/ramdisk")
        print("Create one by executing:")
        print("sudo mkdir -m 777 --parents /tmp/ramdisk")
        print("sudo mount -t tmpfs -o size=9G ramdisk /tmp/ramdisk")
        exit()
    db_root_path = get_api_url()
    db_list = ['uniref90', 'mgnify', 'small_bfd']
    if pairing:
        db_list.append('uniprot')
    db_configs = init_db_configs_api(db_root_path, db_list)
    if debug:
        for i in range(len(db_configs)):
            db_configs[i]['num_streamed_chunks'] = 1
    msa_dict = {}
    for seq_name, seq_str in seq_dict.items():
        print(seq_name)
        output_dir = os.path.join(output_path, f'{seq_name}', 'api')
        os.makedirs(output_dir, exist_ok=True)
        seq_msas = get_MSAs_for_seq_api(seq_str,
                                        output_dir,
                                        db_configs,
                                        hmmer_binary_path)
        msa_dict[seq_name] = seq_msas
    return msa_dict


def get_MSAs_std(seq_dict,
                 db_root_path,
                 output_path='./',
                 hmmer_binary_path='/usr/bin/jackhmmer',
                 hhblits_binary_path='/storage/bin/hh-install/bin/hhblits',
                 small_bfd=False):
    if small_bfd:
        db_list = ['uniref90', 'mgnify', 'small_bfd', 'uniprot']
    else:
        db_list = ['uniref90', 'mgnify', 'bfd_uniclust', 'uniprot']
    db_configs = init_db_configs(db_root_path, db_list)
    msa_dict = {}
    for seq_name, seq_str in seq_dict.items():
        print(seq_name)
        output_dir = os.path.join(output_path, f'{seq_name}', 'std')
        print(output_dir)
        os.makedirs(output_dir, exist_ok=True)
        seq_msas = get_MSAs_for_seq_std(seq_str,
                                        output_dir,
                                        db_configs,
                                        hmmer_binary_path,
                                        hhblits_binary_path,
                                        small_bfd=small_bfd)
        msa_dict[seq_name] = seq_msas
    return msa_dict

def get_MSAs_single_seq(seq_dict, output_path='./'):
    msa_dict = {}
    for seq_name, seq_str in seq_dict.items():
        msa_dict[seq_name] = {'unpairable': [
                                Msa(sequences=[seq_str],
                                deletion_matrix=[[0]*len(seq_str)],
                                descriptions=[seq_name])]}
    return msa_dict

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Query sequence databases to obtain MSAs in a3m format.')
    parser.add_argument('-fas', '--fasta_dir',
                        help='Directory containing the fasta files to be queried. \
                              Only the first record from each fasta file will be processed',
                        default='./')
    parser.add_argument('-out', '--output_dir',
                        help='Directory where the results will be saved',
                        default='./')
    parser.add_argument('-db', '--db_root_dir',
                        help='The root directory of the MSA dbs',
                        default='./')
    parser.add_argument('-api', '--google_api', action='store_true')
    parser.add_argument('-small', '--small_bfd', action='store_true')
    args = parser.parse_args()

    seq_dict = {}
    fasta_files = glob(args.fasta_dir.strip('/') + '/*.fa')
    fasta_files += glob(args.fasta_dir.strip('/') + '/*.fas')
    fasta_files += glob(args.fasta_dir.strip('/') + '/*.fasta')

    for filename in fasta_files:
        print(filename)
        label = os.path.splitext(os.path.basename(filename))[0]
        with open(filename, 'r') as f:
            for record in SeqIO.parse(f, "fasta"):
                seq_dict[label] = record.seq
                break
    if args.google_api:
        msa_dict = get_MSAs_googleapi(seq_dict, output_path=args.output_dir)
    else:
        msa_dict = get_MSAs_std(seq_dict,
                                args.db_root_dir,
                                output_path=args.output_dir,
                                small_bfd=args.small_bfd)
        #print(a3m_dict)

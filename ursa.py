#!/bin/env python

#    Copyright 2013-2015 James S Blachly, MD and The Ohio State University
#
#    This file is part of Ursa miner.
#
#    Ursa miner is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    Ursa miner is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with Ursa miner.  If not, see <http://www.gnu.org/licenses/>.

from __future__ import print_function

import sys
import platform
import argparse

import numpy as np
import pandas as pd
import h5py

class bcolors(object):
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    @classmethod
    def disable(cls):
        cls.HEADER = ''
        cls.OKBLUE = ''
        cls.OKGREEN = ''
        cls.WARNING = ''
        cls.FAIL = ''
        cls.ENDC = ''
        cls.BOLD = ''
        cls.UNDERLINE = ''
        return(True)

class Kdata(object):
    def __init__(self,f):
        self.kversion = ''
        self.idxversion = 0
        self.num_bootstraps = 0
        self.start_time = ''
        self.call = ''
        
        num_targets = len(f['aux']['ids'])
        self.ids =      np.array(f['aux']['ids'])
        self.lengths =  np.array(f['aux']['lengths'])
        self.eff_lengths = np.array(f['aux']['eff_lengths'])
        
        return(None)

def main():
    # Disable ANSI colors in Windows terminal
    if platform.system == 'Windows': bcolors.disable()
    
    parser = argparse.ArgumentParser(description="Ursa miner: kallisto data miner")
    parser.add_argument('mode', choices=['describe', 'extract'])
    parser.add_argument('file')
    args = parser.parse_args()

    f = h5py.File(args.file, 'r')
    
    if args.mode == 'describe': describe(f)
    elif args.mode == 'extract': extract(f)
    else:
        parser.print_usage()
        sys.exit(1)

def get_aux_data(f):
    kdata = Kdata(f)
    kdata.kversion =    f['aux']['kallisto_version'][0]
    kdata.idxversion =  f['aux']['index_version'][0]
    kdata.num_targets = len(f['aux']['ids'])
    kdata.num_bootstrap=f['aux']['num_bootstrap'][0]
    kdata.start_time =  f['aux']['start_time'][0]
    kdata.call =        f['aux']['call'][0]
    return(kdata)

def describe(f):
    #print(bcolors.OKGREEN + "HDF5 groups:" + bcolors.ENDC)
    #for key in f.keys():
    #    print('\t' + f[key].name)
    #print("")

    # Auxillary data
    kdata = get_aux_data(f)
    
    print("")
    print(bcolors.OKGREEN + "Kallisto version     : " + bcolors.ENDC + str(kdata.kversion)      )
    print(bcolors.OKGREEN + "Index version        : " + bcolors.ENDC + str(kdata.idxversion)    )
    print(bcolors.OKGREEN + "Number of targets    : " + bcolors.ENDC + str(kdata.num_targets)   )
    print(bcolors.OKGREEN + "Number of bootstraps : " + bcolors.ENDC + str(kdata.num_bootstrap) )
    print(bcolors.OKGREEN + "Kallisto started     : " + bcolors.ENDC + str(kdata.start_time)    )
    print(bcolors.OKGREEN + "Command              : " + bcolors.ENDC + str(kdata.call)          )
    print("")
    
    

def extract(f):
    # Needs to be here even if not used (no bootstrap) for set subtraction at end
    bs_col_list = []
    
    sys.stderr.write(bcolors.OKGREEN + "Loading data..." + bcolors.ENDC + "\n")
    
    kdata = get_aux_data(f)
    d = { 'target_id': kdata.ids, 'length': kdata.lengths, 'eff_length': kdata.eff_lengths }
    df = pd.DataFrame(d)
    df.set_index('target_id', inplace=True)
    
    est_counts = pd.Series(f['est_counts'], index=df.index)
    df['est_counts'] = est_counts

    if kdata.num_bootstrap > 0:
        # Add all bootstrap count estimates as columns
        for i in range(kdata.num_bootstrap):
            bs_key = 'bs' + str(i)
            bs_col_list.append(bs_key)
            bs_counts = pd.Series(f['bootstrap'][bs_key], index=df.index)
            df[bs_key] = bs_counts
        # Now Compute Statistics
        sys.stderr.write(bcolors.OKGREEN + "Computing statistics..." + bcolors.ENDC + "\n")
        df['est_counts_mean'] = df[bs_col_list].mean(axis=1)
        df['est_counts_stddev'] = df[bs_col_list].std(axis=1)
        df['est_counts_sem'] = df[bs_col_list].sem(axis=1)
        df['est_counts_min'] = df[bs_col_list].min(axis=1)
        df['est_counts_med'] = df[bs_col_list].median(axis=1)
        df['est_counts_max'] = df[bs_col_list].max(axis=1)
    
    # Compute TPM
    if kdata.num_bootstrap == 0:
        divisor = (df.est_counts / df.eff_length ).sum()
        df['tpm'] = ( (df.est_counts / df.eff_length) / divisor ) * 1e6
    else:
        divisor = (df.est_counts_mean / df.eff_length).sum()
        df['tpm'] = ( (df.est_counts_mean / df.eff_length) / divisor) * 1e6
    
    sys.stderr.write(bcolors.OKGREEN + "Writing ursa.csv..." + bcolors.ENDC + "\n")
    #output_cols = list(set(df.columns) - set(bs_col_list))
    #output_cols = ['length','eff_length','est_counts','est_counts_mean','est_counts_stddev','est_counts_min','est_counts_med','est_counts_max']
    output_cols = ['length','eff_length']
    if kdata.num_bootstrap == 0:
        output_cols.append('est_counts')
    else:
        output_cols.extend(['est_counts_mean','est_counts_stddev','est_counts_sem','est_counts_min','est_counts_med','est_counts_max'])
    output_cols.append('tpm')

    df.to_csv('ursa.csv', columns=output_cols)
    return(True)

if __name__ == "__main__":
    main()





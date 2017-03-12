#!/usr/bin/env python
# encoding: utf-8

import lib.HawkEye as HE
import argparse

def main():
    args = parser.parse_args()
    if args.version:
        print "v0.1"
    else:
        print "\n\nWelcome to HawkEye! Let's align something!"
        print HE.grande_alignment(args.DNA)

if __name__ == '__main__':
    parser = argparse.ArgumentParser("HawkEye - sequence alignment for the brave", epilog="Written by Bodie Weedop and Will Pearse", usage="./HawkEye.py sequence_file.fasta")
    parser.add_argument("--version", action="store_true", help="Display version information")
    parser.add_argument("DNA", help="DNA (in FASTA format) to be aligned")
    main()


    

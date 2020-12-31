#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
download-proteome-fastas.py
Used to download whole proteomes from UniProt.
First you need to go to http://www.uniprot.org/proteomes/ and download lists of proteome accessions
    - Can be compressed or uncompressed, as many files as desired
    - Run this script in the same directory
    - Files should be named in format [unique-identifier].tsv(.gz - if compressed)
        - 'tsv' HAS to be present, separated from the identifier by a single period
The actual trick of the URL to use was found here: https://www.biostars.org/p/292993/ (top answer as of 2018-01-29)

__version__ = '0.1.0'
__author__ = 'Jamie Heather'

script downloaded 8/13/19 and modified by Kaden Southard
modified for use with python 3.7.3

example output_file_type: "fasta", "tab", "gff"

"""

import gzip
import os
import requests
import io

#if __name__ == '__main__':

def get_proteome(path,output_file_type, output_path):

    # Get tab files containing proteomes of interest
    all_files = [x for x in os.listdir(path) if '.tab' in x]
    all_files.sort()

    # Define the URL parameters
    url_prefix = 'https://www.uniprot.org/uniprot/?query=proteome:'
    url_suffix= '&format='+output_file_type
    url_tab= '&columns=id,protein_names,genes,length,comment(SUBCELLULAR LOCATION),feature(TOPOLOGICAL DOMAIN),feature(TRANSMEMBRANE)'

    print("Reading in accessions, downloading "+output_file_type+" data...")

    # Loop through all accession tsv files
    for f in all_files:
        print(f)

        # Determine opener
        if f.endswith(".gz"):
            opener = gzip.open
        else:
            opener = open

        #create empty list for accession numbers
        accessions = []

        #the name of the file we are exracting the accessions from
        base_name = f.split('.')[0]

        print('\t' + base_name)

        # Open file, read in accessions
        with opener(path+'/'+f, mode = 'r') as in_file:
            line_count = 0

            #generate list of accessions
            if opener == open:

                for line in in_file:

                    #first line is header
                    if line_count == 0:
                        headers = line.rstrip().split('\t')
                    else:
                        bits = line.rstrip().split('\t')
                        accessions.append(bits[0])

                    line_count += 1

            #compensation for a gzip.open not reading as str
            else:

                for line in in_file:

                    #first line is header
                    if line_count == 0:
                        headers = line.decode('utf8').rstrip().split('\t')
                    else:
                        bits = line.decode('utf8').rstrip().split('\t')
                        accessions.append(bits[0])

                    line_count += 1

        for accession in accessions:

            #create individual file name for proteome
            out_name = base_name + accession + '.'+output_file_type

            print("retrieving "+accession+" data...")

            # Determine full URL, pull the data and write to output file
            if output_file_type == 'tab':
                #need to specify columns for tab format download
                url = requests.get(url_prefix + accession + url_tab + url_suffix)
                f= open(output_path+'/'+out_name, 'w')
                f.write(url.text)
                f.close()

            else:
                url = requests.get(url_prefix + accession + url_suffix)
                f= open(output_path+'/'+out_name, 'w')
                f.write(url.text)
                f.close()

            #for url_line in url:
             #   gzip.open(path+'/'+out_name, 'w').write(url_line)

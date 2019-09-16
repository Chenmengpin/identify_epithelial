
# coding: utf-8

# # Analysis of 10X Molecule Information from count_BCRTCR_f5500
# 
# ## Setup
# 
# To run the notebook, download data from `/share/ScratchGeneral/sunwu/Chromium_10X/BCRTCR/count_BCRTCR_f5500/outs` and put in `outs` in this working directory. 

# In[1]:

import collections
import scipy.sparse as sp_spare
import tables
import pandas as pd
import sys

home_dir = '/share/ScratchGeneral/jamtor/'
project_dir = home_dir + 'projects/single_cell/identify_epithelial/results/cellranger/CID4386'
in_dir = project_dir, 'results/cellranger/' + sample_name + '/'


# In[113]:

MOLECULE_INFO_COLS = ['barcode', 
                      'barcode_corrected_reads', 
                      'conf_mapped_uniq_read_pos', 
                      'gem_group', 
                      'gene',
                      'genome',
                      'nonconf_mapped_reads',
                      'reads',
                      'umi',
                      'umi_corrected_reads',
                      'unmapped_reads']

MOLECULE_INFO_POS = dict((x, i) for (i, x) in enumerate(MOLECULE_INFO_COLS))

MOLECULE_REF_COLS = ['gene_ids', 'gene_names', 'genome_ids']


# In[127]:

# input and output files
molecule_info = in_dir + 'molecule_info.h5'
out_molecule_table = in_dir + 'UMI_dist.txt'


# In[105]:

CELL_BARCODE_LENGTH = 16
UMI_LENGTH = 10
DNA = ['A', 'C', 'G', 'T']

# the decoding scheme is based on knowing the length of the barcode ahead of time
# thus, for example, the number 0 is AAAAAAAAAAAAAAAA rather than A for the cell barcode

def decode(enc, barcode_length):
    ans = [] # append so O(1) instead of O(n)
    for i in range(barcode_length):
        base = DNA[enc % 4]
        ans.append(base)
        enc = enc // 4
    return ''.join(reversed(ans))


# In[128]:

with open(out_molecule_table, 'w') as outfile, tables.open_file(molecule_info, 'r') as infile:

    # total number of molecules
    num_molecules = len(getattr(infile.root, MOLECULE_INFO_COLS[0]).read()) 

    # subsample infile:
    to_keep = 4/5 * num_molecules
    infile = 
    
    # create a dict for looking up genes where { i -> (gene_id, gene_names) }
    gene_ids = getattr(infile.root, 'gene_ids').read()
    gene_names = getattr(infile.root, 'gene_names').read()
    gene_dict = dict((i, [x.decode()]) for (i, x) in enumerate(gene_ids))    
    for (i, x) in enumerate(gene_names):
        gene_dict[i].append(x.decode())
        gene_dict[i] = tuple(gene_dict[i])
    #print(gene_dict)

    
    # total number of genes
    num_genes = len(gene_ids)
    
    # create a dummy gene for max gene index + 1: this row describes reads that did not map confidently to any gene
    gene_dict[num_genes] = (None, None)
    
    # print header
    print('{},{},{},{},{}'.format('BARCODE', 'UMI','READS', 'GENE_ID', 'GENE_NAME'), 
          file=outfile)

    # iterate through each molecule
    for i in range(num_molecules):
        vector = [getattr(infile.root, x)[i] for x in MOLECULE_INFO_COLS]
        gene_id, gene_name = gene_dict[int(vector[MOLECULE_INFO_POS['gene']])]
        print('{},{},{},{},{}'
              .format(decode(int(vector[MOLECULE_INFO_POS['barcode']]), CELL_BARCODE_LENGTH), # decode cell barcode
                      decode(int(vector[MOLECULE_INFO_POS['umi']]), UMI_LENGTH), # decode umi
                      vector[MOLECULE_INFO_POS['reads']],
                      gene_id, 
                      gene_name),
              file=outfile)


# In[110]:

# validate the first entry in barcodes.tsv: AAACCTGAGAAACCAT-1
#decode(24674387, CELL_BARCODE_LENGTH)


# In[29]:

#AAACCTGAGAAACCAT

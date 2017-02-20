from Bio import AlignIO
import numpy as np
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator
import dendropy
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
import matplotlib.pyplot as plt
import scipy.cluster.hierarchy as hac

def aligner ():
    align = AlignIO.read(open('align.fasta'), 'fasta')
    return(align)

def raw_calc (aligned_seqs):
    raw_calculator = DistanceCalculator('identity')
    distance_matrix = raw_calculator.get_distance(aligned_seqs)
    return(distance_matrix)
    
def corr_calc (aligned_seqs):
    corr_calculator = DistanceCalculator('blosum62')
    correcting_matrix = corr_calculator.get_distance(aligned_seqs)
    return(correcting_matrix)

def saturation (d_matrix, c_matrix):
    d_mat = np.array(d_matrix)
    c_mat = np.array(c_matrix)
    diff_matrix = d_mat - c_mat
    rows = diff_matrix.shape[0]
    cols = diff_matrix.shape[1]
    for i in range(0, rows):
        for j in range(0, cols):
            med = np.median(abs(diff_matrix[i,j]-np.median(diff_matrix)))
    mad = 1.4826 * med
    return mad
saturation(raw_calc(aligner()), corr_calc(aligner()))

def sci_hier (d_matrix, c_matrix):
    d_mat = np.array(d_matrix)
    c_mat = np.array(c_matrix)
    diff_matrix = d_mat - c_mat
    hier = hac.linkage(diff_matrix)
    dend = hac.fcluster(hier,1.1505)
    for j in range(0, max(dend)+1):
        for i in [i for i,x in enumerate(dend) if x == j]:
            if j == 1:
                print(i)
sci_hier(raw_calc(aligner()), corr_calc(aligner()))
    

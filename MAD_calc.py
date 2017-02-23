from Bio import AlignIO
from Bio.Align.Applications import MafftCommandline
from Bio import SeqIO
import numpy as np
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator
import dendropy
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
import matplotlib.pyplot as plt
import scipy.cluster.hierarchy as hac
from Bio.Seq import Seq

def starting_pt (file):
    records = list(SeqIO.parse(file, "fasta"))
    return records
full_seqs = starting_pt("birds_shortened.fasta")

def align(seqs):
    
    SeqIO.write(seqs,"very_unlikely_to_be_called_this.fasta","fasta")
    file = "very_unlikely_to_be_called_this.fasta"
    in_file = "/home/god/Desktop/HawkEye/" + file
    mafft_cline = MafftCommandline(input=in_file)
    stdout, stderr = mafft_cline()
    handle = open(file, "w")
    handle.write(stdout)
    handle.close()
    upper_file = "VERY_UNLIKELY_TO_CALLED_THIS1.fasta"
    records = (rec.upper() for rec in SeqIO.parse(file, "fasta"))
    SeqIO.write(records, upper_file, "fasta")
    # parse the sequences and save in a list
    aligned_list = AlignIO.read(open(upper_file), 'fasta')
    return aligned_list
align_seqs = align(full_seqs)

def raw_calc (aligned_seqs):
    raw_calculator = DistanceCalculator('identity')
    distance_matrix = raw_calculator.get_distance(aligned_seqs)
    return(distance_matrix)
raw_mat = raw_calc(align_seqs)

def corr_calc (aligned_seqs):
    corr_calculator = DistanceCalculator('blosum62')
    correcting_matrix = corr_calculator.get_distance(aligned_seqs)
    return(correcting_matrix)
corr_mat = corr_calc(align_seqs)

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
saturation(raw_mat, corr_mat)

def sci_hier (d_matrix, c_matrix):
    d_mat = np.array(d_matrix)
    c_mat = np.array(c_matrix)
    diff_matrix = d_mat - c_mat
    hier = hac.linkage(diff_matrix)
    dend = hac.fcluster(hier,0.08, criterion = 'distance')
    num_clusters = max(dend)+1
    cluster = [[] for x in range(num_clusters)]
    k = 1
    while k < num_clusters:
        for j in range(num_clusters):
            for i in [i for i,x in enumerate(dend) if x == j]:
                if j == k:
                    cluster[j].append(i)
        k+=1
    cluster.remove(cluster[0])
    return cluster
output = sci_hier(raw_mat, corr_mat)




        




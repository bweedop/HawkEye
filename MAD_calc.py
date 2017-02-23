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

def align(seqs):
    SeqIO.write(seqs,"very_unlikely_to_be_called_this.fasta","fasta")
    
    file = "very_unlikely_to_be_called_this.fasta"
    in_file = "/home/god/Desktop/HawkEye/" + file
    mafft_cline = MafftCommandline(input=in_file)
    stdout, stderr = mafft_cline()
    handle = open(file, "w")
    handle.write(stdout)
    handle.close()
    #the mafft alignment put thesequences in lower case. For the rest of this
    #   program, upper case lettering is needed.
    upper_file = "VERY_UNLIKELY_TO_CALLED_THIS1.fasta"
    records = (rec.upper() for rec in SeqIO.parse(file, "fasta"))
    SeqIO.write(records, upper_file, "fasta")
    #parse the sequences and save in a list
    aligned_list = AlignIO.read(open(upper_file), 'fasta')
    return aligned_list

def raw_calc (aligned_seqs):
    raw_calculator = DistanceCalculator('identity')
    distance_matrix = raw_calculator.get_distance(aligned_seqs)
    return(distance_matrix)

def corr_calc (aligned_seqs):
    corr_calculator = DistanceCalculator('blosum62')
    correcting_matrix = corr_calculator.get_distance(aligned_seqs)
    return(correcting_matrix)

def sat_test (d_matrix, c_matrix):
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

def saturation_wrap (saturated):
    raw_mat = raw_calc(saturated)
    corr_mat = corr_calc(saturated)
    saturation = sat_test(raw_mat, corr_mat)
    return saturation

def hawk_wrap (file):
    merged_list = []
    raw_seqs = starting_pt(file)
    full_align = align(raw_seqs)
    whole_saturation = saturation_wrap(full_align)
    if whole_saturation > 0.01:
        output = sci_hier(raw_calc(full_align), corr_calc(full_align))
        saturated = True
        while saturated:
                num_clusters = len(output)
                for i in range(num_clusters):
                    if len(output[i]) > 1:
                        cluster = output[i]
                        cluster_seqs = [full_align[x] for x in cluster]
                        cluster_aligned = align(cluster_seqs)
                        cluster_sat = saturation_wrap(cluster_aligned)
                        if cluster_sat > 0.01:
                            saturated = True
                            output = sci_hier(raw_calc(cluster_aligned),
                                              corr_calc(cluster_aligned))
                            return output
                        elif cluster_sat < 0.01:
                            merged_list.append(cluster_seqs)
                    elif len(output[i]) <= 1:
                        cluster = output[i]
                        cluster_seqs = [full_align[x] for x in cluster]
                        merged_list.append(cluster_seqs)
                    saturated = False
    elif whole_saturation < 0.01:
        merged_list.append(full_align)
        return merged_list
    SeqIO.write(merged_list, "unlikely_merged_file.fasta", "fasta")
seqs = hawk_wrap("birds_shortened.fasta")




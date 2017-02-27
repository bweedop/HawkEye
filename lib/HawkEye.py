# Headers
from Bio import AlignIO, SeqIO, Phylo
from Bio.Align.Applications import MafftCommandline
from Bio.Seq import Seq
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor

import tempfile
import numpy as np
import dendropy
import matplotlib.pyplot as plt
import scipy.cluster.hierarchy as hac

# Functions

def starting_pt (file):
    records = list(SeqIO.parse(file, "fasta"))
    return records

def align(seqs):
    with tempfile.NamedTemporaryFile() as lower_file:
        SeqIO.write(seqs,lower_file.name,"fasta")
        mafft_cline = MafftCommandline(input=lower_file.name)
        stdout, stderr = mafft_cline()
        handle = open(lower_file.name, "w")
        handle.write(stdout)
        handle.close()
        
        #the mafft alignment put thesequences in lower case. For the rest of this
        #   program, upper case lettering is needed.
        with tempfile.NamedTemporaryFile() as upper_file:
            records = (rec.upper() for rec in SeqIO.parse(lower_file.name, "fasta"))
            SeqIO.write(records, upper_file.name, "fasta")
            #parse the sequences and save in a list
            aligned_list = AlignIO.read(open(upper_file.name), 'fasta')
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
    
    dend = hac.fcluster(hier,float(np.mean(hier[:,[2]])), criterion = 'distance')
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

def alignment_wrap(cluster):
    aligned_cluster = align(cluster)
    return aligned_cluster

def saturation_wrap (saturated):
    raw_mat = raw_calc(saturated)
    corr_mat = corr_calc(saturated)
    saturation = sat_test(raw_mat, corr_mat)
    return saturation

def hawk_wrap (file):
    merged_list = []
    raw_seqs = starting_pt(file)
    full_align = alignment_wrap(raw_seqs)
    saturation = saturation_wrap(full_align)
    if saturation > 0.01:
        output = sci_hier(raw_calc(full_align), corr_calc(full_align))
        saturated_list = output
        while saturated_list:
            current_cluster = saturated_list.pop(0)
            if len(current_cluster) > 1:
                cluster_seqs = [full_align[x] for x in current_cluster]
                cluster_aligned = alignment_wrap(cluster_seqs)
                cluster_sat = saturation_wrap(cluster_aligned)
                if cluster_sat > 0.01:
                        new_clusters = sci_hier(raw_calc(cluster_aligned),
                                          corr_calc(cluster_aligned))
                        if len(new_clusters) == 1:
                            problem_seqs = new_clusters[0]
                            while problem_seqs:
                                solution = []
                                solution.append(problem_seqs.pop(0))
                                single_seqs = [current_cluster[x]
                                               for x in solution]
                                merged_list.append(single_seqs)
                        elif len(new_clusters) >= 1:
                            while new_clusters:
                                separate = new_clusters.pop(0)
                                new_cluster_seqs = [current_cluster[x]
                                                    for x in separate]
                                saturated_list.append(new_cluster_seqs)
                elif cluster_sat < 0.01:
                        merged_list.append(current_cluster)
            elif len(current_cluster) <= 1:
                    merged_list.append(current_cluster)
        return merged_list
    elif whole_saturation < 0.01:
        merged_list.append(full_align)
        return merged_list



# Headers
from Bio import AlignIO, SeqIO, Phylo
from Bio.Align.Applications import MafftCommandline
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import SingleLetterAlphabet
from Bio.Emboss.Applications import FDNADistCommandline as dist
import os
import subprocess

import tempfile
import numpy as np
import dendropy
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
    with tempfile.NamedTemporaryFile() as align_file:
        SeqIO.write(aligned_seqs,align_file.name,"fasta")
        with tempfile.NamedTemporaryFile() as upper_file:
            
            dnadist_cline = dist(sequence=align_file.name,method = 's', outfile=
                                 upper_file.name)
            stdout, stderr = dnadist_cline()
            lines = upper_file.read()
            l = []
            for t in lines.split():
                try:
                    l.append(float(t))
                except ValueError:
                    pass
        array_dim = int(l.pop(0))
        sim_raw_array = np.array(l)
        distance_raw_array = 1 - sim_raw_array
        distance_raw_array.resize((array_dim, array_dim))
        return distance_raw_array

def corr_calc (aligned_seqs):
    with tempfile.NamedTemporaryFile() as align_file:
        SeqIO.write(aligned_seqs,align_file.name,"fasta")
        with tempfile.NamedTemporaryFile() as upper_file:
            
            dnadist_cline = dist(sequence=align_file.name ,method = 'j', outfile=
                                 upper_file.name)
            stdout, stderr = dnadist_cline()
            lines = upper_file.read()
            l = []
            for t in lines.split():
                try:
                    l.append(float(t))
                except ValueError:
                    pass
    array_dim = int(l.pop(0))
    corr_array = np.array(l)
    corr_array.resize((array_dim, array_dim))
    return corr_array

def sat_test (d_matrix, c_matrix):
    diff_matrix = d_matrix - c_matrix
    rows = diff_matrix.shape[0]
    cols = diff_matrix.shape[1]
    for i in range(0, rows):
        for j in range(0, cols):
            med = np.median(abs(diff_matrix[i,j]-np.median(diff_matrix)))
    mad = 1.4826 * med
    return mad

def sci_hier (d_matrix, c_matrix):
    diff_matrix = d_matrix - c_matrix
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
    whole_saturation = saturation_wrap(full_align)
    if whole_saturation > 0.01:
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
                            if len(new_clusters[0]) > 2:
                                print("ERROR: Issue with separating clusters")
                                break
                            elif len(new_clusters[0]) <= 2:
                                
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

def clusters_alignment (file):
    full_alignment = starting_pt(file)
    clusters = hawk_wrap(file)
    consensus_list = []
    cluster_index = []
    for x in range(len(clusters)):
        if len(clusters[x]) > 1:
            cluster_index.append(x)
    while clusters:
        current_cluster = clusters.pop(0)
        if len(current_cluster) > 1:
            multiple_seqs = [full_alignment[x] for x in current_cluster]
            with tempfile.NamedTemporaryFile() as alignment_file:
                SeqIO.write(multiple_seqs,alignment_file.name,"fasta")
                with tempfile.NamedTemporaryFile() as consensus_file:
                    subprocess.call(["em_cons", alignment_file.name,consensus_file.name])
                    seq = consensus_file.name
                    data = open(seq).read()
                    con = []
                    dash = 'n'
                    R_DNA = ['A','a','C','c','T','t','G','g','U','u']
                    for i in data.split():
                        for j in i:
                            if j in R_DNA:
                                con.append(j)
                            elif j not in R_DNA:
                                pass
                    con_str = ''.join(str(k) for k in con)
                simple_seq_r = SeqRecord(Seq(con_str, SingleLetterAlphabet()),
                                         id = "CLUSTER"+str(cluster_index.pop(0)))
                consensus_list.append(simple_seq_r)
        elif len(current_cluster) == 1:
            single_seq = [full_alignment[x] for x in current_cluster]
            consensus_list.append(single_seq.pop(0))
    with tempfile.NamedTemporaryFile() as consensus_alignment: 
        SeqIO.write(consensus_list,consensus_alignment.name,"fasta")
        file = consensus_alignment.name
        in_file = consensus_alignment.name
        mafft_cline = MafftCommandline(input=in_file)
        stdout, stderr = mafft_cline()
        handle = open(file, "w")
        handle.write(stdout)
        handle.close()
        path = consensus_alignment.name
        data = open(path).read()
    return data


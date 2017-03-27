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
    dend = hac.fcluster(hier,float(2*np.std(hier[:,[2]])), criterion = 'distance')
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
        onlyCluster = list(range(0, len(full_align)))
        return onlyCluster

def clusters_alignment (file):
    full_alignment = starting_pt(file)
    clusters = hawk_wrap(file)
    consensus_list = []
    cluster_index = []
    if type(clusters[0]) == list:
        for x in range(len(clusters)):
            if len(clusters[x]) > 1:
                cluster_index.append(x)
    while clusters:
        if type(clusters[0]) == list:
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
                cluster_seq_id = []
                while multiple_seqs:
                     seq_info = multiple_seqs.pop(0)
                     cluster_seq_id.append(seq_info.id)

                seqid_string = '|'.join(str(l) for l in cluster_seq_id)
                simple_seq_r = SeqRecord(Seq(con_str, SingleLetterAlphabet()),
                                         id = "CLUSTER_"+str(cluster_index.pop(0))
                                         +": "+seqid_string)
                consensus_list.append(simple_seq_r)
            elif len(current_cluster) == 1:
                single_seq = [full_alignment[x] for x in current_cluster]
                consensus_list.append(single_seq.pop(0))
        elif type(clusters[0]) == int:
            noSaturation = alignment_wrap(full_alignment)
            return noSaturation
    final_aligned_clusters = []
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
        with tempfile.NamedTemporaryFile() as clusters_file:
            records = (rec.upper() for rec in SeqIO.parse(consensus_alignment.name,
                                                          "fasta"))
            SeqIO.write(records, clusters_file.name, "fasta")
            aligned_list = AlignIO.read(open(clusters_file.name), 'fasta')
    return aligned_list

def grande_alignment (file):
    raw_seqs = starting_pt(file)
    before_segment = clusters_alignment(file)
    list_of_clusters = hawk_wrap(file)
    print("Clusters formed:")
    print(list_of_clusters)
    seqs_in_consensus = []
    allEqualSeqs = []
    position = 0
    if type(list_of_clusters[0]) == list:
        big_final_alignment = []
        while list_of_clusters:
            current_cluster = list_of_clusters.pop(0)
            if len(current_cluster) == 1:
                single_seq = before_segment[position]
                big_final_alignment.append(single_seq)
                position += 1
            elif len(current_cluster) > 1:
                seqs_in_the_cluster = [raw_seqs[x] for x in current_cluster]
                multiple_seq = before_segment[position]
                while seqs_in_the_cluster:
                    seqs_in_consensus.append(seqs_in_the_cluster.pop(0))
                    with tempfile.NamedTemporaryFile() as segment_align:
                        SeqIO.write(seqs_in_consensus,segment_align.name, "fasta")
                        original_seqs = list(SeqIO.parse(segment_align.name, "fasta"))
                        with tempfile.NamedTemporaryFile() as segmenter:
                            dash = '-'
                            dashes = []
                            consensus = multiple_seq.seq
                            seq = original_seqs
                            for n in range(len(seq)):
                                seq_str = seq[n].seq
                                seq_id = seq[n].id
                                able_to_insert_seq = seq_str.tomutable()
                            for x in consensus:
                                if x == dash:
                                    dashes = [y for y, x in enumerate(consensus) if x == dash]
                            while dashes:
                                dash_position = dashes.pop(0)
                                able_to_insert_seq.insert(dash_position, '-')
                            new_seq_record = SeqRecord(Seq(str(able_to_insert_seq),
                                                           SingleLetterAlphabet()),
                                                       id = seq_id)
                            big_final_alignment.append(new_seq_record)       
                position += 1
        with tempfile.NamedTemporaryFile() as unequalSeqs:
            dash = '-'
            SeqIO.write(big_final_alignment,unequalSeqs.name, "fasta")
            big_final_alignment = list(SeqIO.parse(unequalSeqs.name, "fasta"))
            largestSeq = len(max([big_final_alignment[ind].seq
                                  for ind in range(len(big_final_alignment))]))
            while big_final_alignment:
                checkSeq = big_final_alignment.pop(0)
                if len(checkSeq.seq) < largestSeq:
                    smallerLength = len(checkSeq.seq)
                    seqStr = checkSeq.seq
                    seqId = checkSeq.id
                    needsEndingFilled = seqStr.tomutable()
                    endingDashes = list(range(smallerLength+1, largestSeq+1))
                    while endingDashes:
                        j = endingDashes.pop(0)
                        needsEndingFilled.insert(j, '-')
                        nowEqual = SeqRecord(Seq(str(needsEndingFilled),
                                                 SingleLetterAlphabet()),
                                             id = seqId)
                    allEqualSeqs.append(nowEqual)
                elif len(checkSeq.seq) == largestSeq:
                    allEqualSeqs.append(checkSeq)
    elif type(list_of_clusters[0]) == int:
        allEqualSeqs = before_segment
    for seqs in range(len(allEqualSeqs)):
        print(">"+allEqualSeqs[seqs].id)
        print(allEqualSeqs[seqs].seq) 












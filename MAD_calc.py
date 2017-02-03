from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio import AlignIO
align = AlignIO.read(open('alignment.fasta'), 'fasta')


calculator = DistanceCalculator('identity')
distance_matrix = calculator.get_distance(align)


calculator = DistanceCalculator('blosum62')
correcting_matrix = calculator.get_distance(align)


def saturation (d_matrix, c_matrix):
    import numpy as np
    d_mat = np.array(distance_matrix)
    c_mat = np.array(correcting_matrix)
    diff_matrix = d_mat - c_mat
    rows = diff_matrix.shape[0]
    cols = diff_matrix.shape[1]
    for i in range(0, rows):
        for j in range(0, cols):
            med = np.median(abs(diff_matrix[i,j]-np.median(diff_matrix)))
    mad = 1.4826 * med
    print(mad)
saturation(distance_matrix, correcting_matrix)

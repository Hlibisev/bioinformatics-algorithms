import numpy as np
import pandas as pd


def value_from_matrix(protein1, protein2, matrix):
    ind1 = list(matrix.columns).index(protein1)
    ind2 = list(matrix.columns).index(protein2)

    return matrix.iloc[ind1, ind2]


def get_wunsch_matrix(protein1, protein2, delta=-4, BLOSUM_matrix=None):
    """
    :type protein1: list
    :type protein2: list
    :type delta: int
    :type mu: int
    """
    matrix = np.zeros((len(protein1) + 1, len(protein2) + 1))
    matrix[:, 0] = np.array([-delta * i for i in range(len(protein1) + 1)])
    matrix[0, :] = np.array([-delta * i for i in range(len(protein2) + 1)])

    for i in range(matrix.shape[0] - 1):
        for j in range(matrix.shape[1] - 1):
            comp = value_from_matrix(protein1[i], protein2[j], BLOSUM_matrix)
            matrix[i + 1, j + 1] = max(matrix[i, j] + comp,
                                       matrix[i, j + 1] + delta,
                                       matrix[i + 1, j] + delta)
    return matrix


def sequence_alignment(protein1, protein2, delta=-4, BLOSUM_matrix=None):
    """
    :type protein1: string
    :type protein2: string
    :type delta: int
    :type mu: int
    """
    alignment1, alignment2 = [], []
    protein1, protein2 = list(protein1), list(protein2)
    matrix = get_wunsch_matrix(protein1, protein2, delta, BLOSUM_matrix)

    i, j = len(protein1) - 1, len(protein2) - 1
    while i > -1 and j > -1:
        comp = value_from_matrix(protein1[i], protein2[j], BLOSUM_matrix)

        if matrix[i + 1, j + 1] == matrix[i, j + 1] + delta:
            alignment1.append(protein1[i])
            alignment2.append("_")
            i -= 1

        elif matrix[i + 1, j + 1] == matrix[i + 1, j] + delta:
            alignment1.append("_")
            alignment2.append(protein2[j])
            j -= 1

        elif matrix[i + 1, j + 1] == matrix[i, j] + comp:
            alignment1.append(protein1[i])
            alignment2.append(protein2[j])
            i -= 1
            j -= 1

    while i > -1:
        alignment1.append(protein1[i])
        alignment2.append("_")
        i -= 1

    while j > -1:
        alignment1.append("_")
        alignment2.append(protein2[j])
        j -= 1

    return "".join(alignment1)[::-1], "".join(alignment2)[::-1]

protein2 = "GCATGCV"
protein1 = "GATTACA"

matrix = pd.read_csv("/Needleman Wunsch Alrogithm/BLOSUM")
print(sequence_alignment(protein1, protein2, -4, matrix))
# Вывод: ('__GATTACA', 'GC_ATG_CV')

matrix.iloc[5,5] = 10 # меняем G -> G на 10, вместо 6
print(sequence_alignment(protein1, protein2, -4, matrix))
# Вывод: ('GATTACA', 'GCATGCV')

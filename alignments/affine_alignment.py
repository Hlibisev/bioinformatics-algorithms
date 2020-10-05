import numpy as np
import pandas as pd
from tools import value_from_matrix


def affine_alignment(protein1, protein2, open_gap_penalty, continue_gap_penalty, weight_matrix):
    """
    :param protein1: first list of strings which we will alignment
    :param protein2: second list of strings which we will alignment
    :param open_gap_penalty:
    :param continue_gap_penalty:
    :param weight_matrix: matrix of weights or list of weight_match and weight_mismatch
    :return: tuple with 2 strings
    """

    # Init 3-leveled Manhattan Grid
    s_down = np.zeros((len(protein1) + 1, len(protein2) + 1))
    s_right = np.zeros((len(protein1) + 1, len(protein2) + 1))
    s_diag = np.zeros((len(protein1) + 1, len(protein2) + 1))

    way = np.zeros((len(protein1) + 1, len(protein2) + 1))

    for matrix in [s_down, s_right, s_diag]:
        matrix[:, 0] = np.array([continue_gap_penalty * i for i in range(len(protein1) + 1)]) + open_gap_penalty
        matrix[0, :] = np.array([continue_gap_penalty * i for i in range(len(protein2) + 1)]) + open_gap_penalty

    # Filling in three matrix's
    for i in range(s_diag.shape[0] - 1):
        for j in range(s_diag.shape[1] - 1):
            comp = value_from_matrix(protein1[i], protein2[j], weight_matrix)

            s_down[i + 1, j + 1] = max(s_down[i, j + 1] + continue_gap_penalty,
                                       s_diag[i, j + 1] + open_gap_penalty + continue_gap_penalty)

            s_right[i + 1, j + 1] = max(s_right[i + 1, j] + continue_gap_penalty,
                                        s_diag[i + 1, j] + open_gap_penalty + continue_gap_penalty)

            s_diag[i + 1, j + 1] = max(s_diag[i, j] + comp,
                                       s_right[i + 1, j + 1],
                                       s_down[i + 1, j + 1])

            # memorizing the best path
            way[i + 1, j + 1] = np.array([s_diag[i, j] + comp, s_right[i + 1, j + 1], s_down[i + 1, j + 1]]).argmax()

    i, j = len(protein1) - 1, len(protein2) - 1
    alignment1, alignment2 = [], []

    while i > -1 and j > -1:
        if way[i + 1, j + 1] == 0:
            alignment1.append(protein1[i])
            alignment2.append(protein2[j])
            i -= 1
            j -= 1

        elif way[i + 1, j + 1] == 1:
            alignment1.append("_")
            alignment2.append(protein2[j])
            j -= 1

        else:
            alignment1.append(protein1[i])
            alignment2.append("_")
            i -= 1

    while i > -1:
        alignment1.append(protein1[i])
        alignment2.append("_")
        i -= 1

    while j > -1:
        alignment1.append("_")
        alignment2.append(protein2[j])
        j -= 1

    return "".join(alignment1)[::-1], "".join(alignment2)[::-1]


if __name__ == "__main__":

    protein1 = "GGGACTGAG"
    protein2 = "GACTA"

    matrix = pd.read_csv("/Users/hlibisev/Documents/GitHub/bioinformatics-algorithms/alignments/BLOSUM")
    print(affine_alignment(protein1, protein2, -4, -1, matrix))
    # Вывод: ('GGGACTGAG', '__GACT__A') - что интересно, потому что алгоритм не хочет открывать 2 раза _

    print(affine_alignment(protein1, protein2, -1, -4, matrix))
    # Вывод: ('GGGACTGAG', '__GACT_A_')


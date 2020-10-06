import numpy as np
import pandas as pd
from tools import get_wunsch_matrix, value_from_matrix


def local_alignment(protein1, protein2, delta=-4, BLOSUM_matrix=None):
    """
    :type protein1: list
    :type protein2: list
    :type delta: int
    :type BLOSUM_matrix: pandas.DataFrame
    :return tuple with 2 strings
    _____________________________________
    Algorithm of local alignment
    """
    alignment1, alignment2 = [], []
    protein1, protein2 = list(protein1), list(protein2)
    matrix = get_wunsch_matrix(protein1, protein2, delta, BLOSUM_matrix, True)

    end1, end2 = np.unravel_index(matrix.argmax(), matrix.shape)
    i, j = end1 - 1, end2 - 1
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

        else:
            break

    while i > -1 or j > -1:
        letter1 = "_" if i <= -1 else protein1[i]
        letter2 = "_" if j <= -1 else protein2[j]

        alignment1.append(letter1)
        alignment2.append(letter2)

        i -= 1
        j -= 1

    while end1 < len(protein1) or end2 < len(protein2):
        letter1 = "_" if end1 >= len(protein1) else protein1[end1]
        letter2 = "_" if end2 >= len(protein2) else protein2[end2]

        alignment1.insert(0, letter1)
        alignment2.insert(0, letter2)

        end1 += 1
        end2 += 1

    return "".join(alignment1)[::-1], "".join(alignment2)[::-1]


if __name__ == "__main__":
    protein1 = "GGGACTGAG"
    protein2 = "GACT"

    matrix = pd.read_csv("/Users/hlibisev/Documents/GitHub/bioinformatics-algorithms/alignments/BLOSUM")
    print(local_alignment(protein1, protein2, -4, matrix))

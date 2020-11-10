import pandas as pd
from tools import get_wunsch_matrix, value_from_matrix


def global_alignment(protein1, protein2, delta=-4, BLOSUM_matrix=None):
    """
    :type protein1: list
    :type protein2: list
    :type delta: int
    :type BLOSUM_matrix: pandas.DataFrame
    :return tuple with 2 strings
    _____________________________________
    Algorithm of global alignment
    """
    alignment1, alignment2 = [], []
    protein1, protein2 = list(protein1), list(protein2)
    matrix = get_wunsch_matrix(protein1, protein2, delta, BLOSUM_matrix, False)
    print(matrix)
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



if __name__ == "__main__":

    protein1 = "SENR"
    protein2 = "ANR"

    matrix = pd.read_csv("/Users/hlibisev/Documents/GitHub/bioinformatics-algorithms/alignments/BLOSUM")
    print(global_alignment(protein1, protein2, -10, matrix))

import numpy as np


def get_wunsch_matrix(protein1: list, protein2: list, delta: float = -4, BLOSUM_matrix = None, local: bool = False):
    """
    :param protein1: first list of strings which we will alignment
    :param protein2: second list of strings which we will alignment
    :param delta: penalty for mismatch
    :param BLOSUM_matrix:
    :param local: if local alignment then True else False
    :return: matrix used for algorithms of alignment
    """

    if local is False:
        local = float("-inf")
    else:
        local = 0

    matrix = np.zeros((len(protein1) + 1, len(protein2) + 1))
    matrix[:, 0] = np.array([delta * i for i in range(len(protein1) + 1)])
    matrix[0, :] = np.array([delta * i for i in range(len(protein2) + 1)])

    for i in range(matrix.shape[0] - 1):
        for j in range(matrix.shape[1] - 1):
            comp = value_from_matrix(protein1[i], protein2[j], BLOSUM_matrix)
            matrix[i + 1, j + 1] = max(local,
                                       matrix[i, j] + comp,
                                       matrix[i, j + 1] + delta,
                                       matrix[i + 1, j] + delta)
    return matrix


def value_from_matrix(protein1, protein2, matrix):
    """
    :param protein1: first letter
    :param protein2: second letter
    :param matrix: matrix of weights or list of weight_match and weight_mismatch
    :return: penalty value
    """

    # If we have weight_match and weight_mismatch
    if matrix is list:
        if protein1 == protein2:
            return matrix[0]
        else:
            return matrix[1]

    # If we have marix, BLOSUM, for example
    ind1 = list(matrix.columns).index(protein1)
    ind2 = list(matrix.columns).index(protein2)

    return matrix.iloc[ind1, ind2]

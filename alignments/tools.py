import numpy as np


def get_wunsch_matrix(protein1, protein2, delta=-4, BLOSUM_matrix=None, local=False):
    """
    ____________________________________
    :type protein1: list
    :type protein2: list
    :type delta: int
    :type BLOSUM_matrix: pandas.DataFrame
    :type local: boolean
    :return pandas.DataFrame
    _____________________________________
    Create matrix used for algorithms of alignment
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

    ind1 = list(matrix.columns).index(protein1)
    ind2 = list(matrix.columns).index(protein2)

    return matrix.iloc[ind1, ind2]

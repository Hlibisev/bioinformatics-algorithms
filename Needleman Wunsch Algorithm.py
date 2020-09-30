import numpy as np

vocab_proteins = list("A R N D C Q E G H I L K M F P S T W Y V B Z X".split())
protein2 = "GCATGCU"
protein1 = "GATTACA"

protein1, protein2 = input().split()
mu = 1
delta = 1


def get_wunsch_matrix(protein1, protein2, delta=1, mu=1):
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
            comp = 1 if protein1[i] == protein2[j] else -mu
            matrix[i + 1, j + 1] = max(matrix[i, j] + comp,
                                       matrix[i, j + 1] - delta,
                                       matrix[i + 1, j] - delta)
    return matrix


def sequence_alignment(protein1, protein2, delta=1, mu=1):
    """
    :type protein1: string
    :type protein2: string
    :type delta: int
    :type mu: int
    """
    alignment1, alignment2 = [], []
    protein1, protein2 = list(protein1), list(protein2)
    matrix = get_wunsch_matrix(protein1, protein2, delta, mu)

    i, j = len(protein1) - 1, len(protein2) - 1

    while i > -1 and j > -1:
        comp = 1 if protein1[i] == protein2[j] else -mu

        if matrix[i + 1, j + 1] == matrix[i, j + 1] - delta:
            alignment1.append(protein1[i])
            alignment2.append("_")
            i -= 1

        elif matrix[i + 1, j + 1] == matrix[i + 1, j] - delta:
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


anw1, anw2 = sequence_alignment(protein1, protein2, delta, mu)
print(anw1, anw2)

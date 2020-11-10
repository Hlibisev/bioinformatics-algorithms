import numpy as np


def def_pair(str1, str2):
    a = ["A", "U", "C", "G"]
    if a.index(str1) + a.index(str2) == 1:  # A, U
        return True
    if a.index(str1) + a.index(str2) == 5:  # C, G
        return True
    return False


def get_nussinov_matrix(RNA, min_loop=3):
    L = len(RNA)
    matrix = np.empty((L, L))
    trace = np.ones((L, L), dtype=int) * -2

    for t in range(0, min_loop):
        for i in range(L - t):
            j = i + t
            matrix[i, j] = 0

    for t in range(min_loop, L):
        for i in range(L - t):
            j = i + t

            jump = max(*[(matrix[i, k] + matrix[k + 1, j], k) for k in range(i + 1, j)], key=lambda x: x[0])

            value, step = max(
                (matrix[i + 1, j], 1),
                (matrix[i, j - 1], -1),
                (matrix[i + 1, j - 1] + 1 if def_pair(RNA[i], RNA[j]) else matrix[i + 1, j - 1], 0),
                jump,
                key=lambda x: x[0]
            )

            trace[i, j] = step
            matrix[i, j] = value
    return matrix, trace


def get_matches(trace, i, j):
    if trace[i, j] == -2:
        return []
    elif trace[i, j] == -1:
        match = get_matches(trace, i, j - 1)
    elif trace[i, j] == 0:
        match = [(i, j)] + get_matches(trace, i + 1, j - 1)
    elif trace[i, j] == 1:
        match = get_matches(trace, i + 1, j)
    elif trace[i, j] > 1:
        match = get_matches(trace, i + trace[i, j], j) + get_matches(trace, i, j - trace[i, j])
    return match


if __name__ == "__main__":
    RNA = "CCCUUUAG"
    matrix, trace = get_nussinov_matrix(RNA)

    print(trace, 0, len(matrix) - 1)




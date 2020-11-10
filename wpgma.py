import numpy as np


class Node():
    def __init__(self, left, right, distance, size):
        self.left = left
        self.right = right
        self.distance = distance
        self.size = size


class Leaf():
    def __init__(self, name, distance, size):
        self.name = name
        self.distance = distance
        self.size = size


class Cluster():
    def __init__(self, left, right, distance):
        self.left = left
        self.right = right
        self.distance = distance


def get_wpgma_string(leave):
    if leave.size == 1:
        return leave.name, 0
    else:
        left, l_dist = get_wpgma_string(leave.left)
        right, r_dist = get_wpgma_string(leave.right)

        return f'({left}:{leave.distance / 2 - l_dist},{right}:{leave.distance / 2 - r_dist})', leave.distance / 2


def wpgma(matrix, leaves, weighted=True):
    leaves = [Leaf(name, 0, 1) for name in leaves]

    while matrix is not None:
        matrix, leaves = wpgma_step(matrix, leaves, weighted)

    string, _ = get_wpgma_string(leaves[0])
    return string


def argmin(matrix):
    a = matrix.argmin()
    i, j = a // len(matrix), a % len(matrix)
    return (i, j)


def wpgma_step(matrix, leaves, weighted=True):
    i, j = argmin(matrix)
    D = matrix[i, j]

    if len(matrix) <= 2:
        leaves.append(Node(leaves[i], leaves[j], D, leaves[i].size + leaves[j].size))
        del leaves[i]
        del leaves[j - 1]
        return None, leaves

    w_i = leaves[i].size if weighted is False else 1
    w_j = leaves[j].size if weighted is False else 1

    a = [[(matrix[min(i, k), max(i, k)] * w_i + matrix[min(j, k), max(j, k)] * w_j) / (w_i + w_j)]
         for k in range(0, len(matrix))
         if k != i and k != j]

    for axis in (0, 1):
        matrix = np.delete(matrix, (j, i), axis=axis)

    matrix = np.append(matrix, a, axis=1)
    matrix = np.append(matrix, [[np.inf for _ in range(matrix.shape[1])]], axis=0)  # attention

    leaves.append(Node(leaves[i], leaves[j], D, leaves[i].size + leaves[j].size))
    del leaves[i]
    del leaves[j - 1]

    return matrix, leaves


if __name__ == "__main__":
    inf = np.inf
    leafs = ['A', 'B', 'C', 'D']
    distances = np.array([[inf, 16, 16, 10],
                          [inf, inf, 8, 8],
                          [inf, inf, inf, 4],
                          [inf, inf, inf, inf]], dtype=np.float32)

    print(wpgma(distances, leafs))
    print(wpgma(distances, leafs, weighted=False))



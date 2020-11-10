import numpy as np
from wpgma import Leaf, Cluster, argmin


def neighbor_joining(matrix, leaves):
    leaves = [Leaf(name, 0, 1) for name in leaves]

    while matrix is not None:
        matrix, leaves = nj_step(matrix, leaves)

    string = get_nj_string(leaves[0])
    return string


def get_nj_string(leaves):
    if type(leaves) is Leaf:
        return leaves.name
    else:
        left = get_nj_string(leaves.left)
        right = get_nj_string(leaves.right)
        l_dist, r_dist = leaves.distance
        return f'({left}:{l_dist},{right}:{r_dist})'


def nj_step(matrix, leaves):

    if len(matrix) == 2:
        return None, [Cluster(leaves[0], leaves[1], (matrix[0, 1] / 2, matrix[0, 1] / 2))]

    def mean_dist(i, j):
        distance = 0
        for k in range(len(leaves)):
            if k != i and k != j:
                value = matrix[min(i, k), max(i, k)]
                distance = distance + value if value != np.inf else distance
        return distance / (len(leaves) - 2)

    def nj_dist(i, j):
        distance = 0
        for k in range(len(leaves)):
            if k != i and k != j:
                value1 = matrix[min(i, k), max(i, k)]
                value2 = matrix[min(j, k), max(j, k)]

                distance = distance + value1 if value1 != np.inf else distance
                distance = distance + value2 if value2 != np.inf else distance

        return matrix[min(i, j), max(i, j)] - distance / (len(leaves) - 2)

    dists_matrix = np.zeros((len(leaves), len(leaves))) + np.inf

    for i in range(len(leaves)):
        for j in range(i, len(leaves)):
            dists_matrix[i, j] = nj_dist(i, j)


    i, j = argmin(matrix)
    distance = matrix[i, j]
    dist_U_A = mean_dist(i, j)
    dist_U_B = mean_dist(j, i)
    dist_U_X = []

    # counting new distance
    for k in range(len(leaves)):
        if k != i and k != j:
            dist_U_X.append([(matrix[min(i, k), max(i, k)] + matrix[min(j, k), max(j, k)] - distance) / 2])

    for axis in (0, 1):
        matrix = np.delete(matrix, (j, i), axis=axis)

    # print(matrix.shape, len(dist_U_X))
    matrix = np.append(matrix, dist_U_X, axis=1)
    matrix = np.append(matrix, [[np.inf for _ in range(matrix.shape[1])]], axis=0)  # attention

    leaves.append(Cluster(leaves[i], leaves[j],
                          ((distance + dist_U_A - dist_U_B) / 2,
                          (distance + dist_U_B - dist_U_A) / 2)))

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

    print(neighbor_joining(distances, leafs))

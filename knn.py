import matplotlib.pyplot as plt
import numpy as np
from sklearn.neighbors import NearestNeighbors

# X = np.array(
#     [
#         [-1, -1],
#         [-2, -1],
#         [-3, -2],
#         [1, 1],
#         [2, 1],
#         [3, 2],
#     ]
# )

X = np.array(
    [
        [0, 0],
        [1, 2],
        [5, 3],
        [2, 2],
        [15, 30],
        [1, 20],
    ]
)

nbrs = NearestNeighbors(n_neighbors=2, algorithm="ball_tree").fit(X)

distances, indices = nbrs.kneighbors(X)
print(indices)
print(distances)

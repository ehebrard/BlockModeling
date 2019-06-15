from __future__ import print_function
import torch
import  numpy as np
import torch.nn.functional as F


import random
from math import log

rand = random.Random()
rand.seed(42)

class Graph:
    def __init__(self, V):
        self.V = V
        self.adj_list = [[] for i in range(V)]

    def add_edge(self, i, j):
        self.adj_list[i].append(j)

    def adj(self, i):
        return self.adj_list[i]

    def adj_matrix(self):
        matrix = [[0] * self.V for i in range(self.V)]
        for i in range(self.V):
            for j in self.adj_list[i]:
                matrix[i][j] = 1
        return matrix
    def E(self):
        return sum([len(adj) for adj in self.adj_list])

class BlockModel:
    def __init__(self, K):
        self.K = K
        self.M = [[[0,0.7][rand.randint(0,1)]+rand.random()/4  for i in range(K)] for j in range(K)]
        #self.M = [[1.0 if i==j else 0.0 for i in range(K)] for j in range(K)]
    def random_graph(self,nodes):
        g = Graph(sum(nodes))
        # block of each node
        block = [e for i in range(len(nodes)) for e in [i]*nodes[i]]
        for i in range(g.V):
            for j in range(g.V):
                if rand.random() <= self.M[block[i]][block[j]]:
                    g.add_edge(i,j)
        return g


def bit_length(g,block):
    K = max(block)+1
    n = [0]*K
    for b in block:
        n[b] += 1
    d = [[0]*K for i in range(K)]
    for i in range(g.V):
        for j in g.adj(i):
            d[block[i]][block[j]] += 1
    tot_len = 0.0
    for b1 in range(K):
        for b2 in range(K):
            n1,n2 = n[b1],n[b2]
            n12 = d[b1][b2]
            if n12 != n1*n2 and n12 != 0:
                tot_len += - n12 * log(n12/(n1*n2)) - (n1*n2-n12) * log(1-n12/(n1*n2))
    return tot_len


nodes = [10,5,7]

n = sum(nodes)
K = len(nodes)
model = BlockModel(K)
g = model.random_graph(nodes)
block_opt = [e for i in range(K) for e in [i] * nodes[i]]
opt_length = bit_length(g, block_opt)


adj_matrix = torch.from_numpy(np.array(g.adj_matrix(),dtype=np.float32))

print(adj_matrix)

# initial best weights
#weights = torch.tensor(np.array([[1 if j==block_opt[i] else 0 for j in range(K)] for i in range(n)],dtype=np.float32),requires_grad=True)

# random initial weigths
weights = torch.randn(sum(nodes),K,requires_grad=True)

optimizer = torch.optim.Adam([weights], lr=0.0001)

for iter in range(30000):

    optimizer.zero_grad()
    weights_norm = F.softmax(weights, dim=1)
    # number of nodes in each block
    n_nodes = torch.sum(weights_norm, dim=0)
    # matrix with the number of edges between the blocks
    n_edges = [[torch.sum(torch.mul(torch.ger(weights_norm[:, i], weights_norm[:, j]), adj_matrix)) for j in range(K)] for i in range(K)]

    objective = torch.tensor(.0)
    for i in range(K):
        for j in range(K):
            # computes the number of bits between the block i and j
            objective -= n_edges[i][j] * torch.log(n_edges[i][j] / (n_nodes[i] * n_nodes[j]))
            objective -= (n_nodes[i] * n_nodes[j] - n_edges[i][j]) * torch.log(1 - n_edges[i][j] / (n_nodes[i] * n_nodes[j]))

    # gradient descent step
    objective.backward(retain_graph=True)
    optimizer.step()

    if iter%1000 == 0:
        print("Objective:", objective)
        print("assignment vector of nodes to blocks")
        print(weights_norm)
        print("block sizes",n_nodes)
        print("number of edges between the blocks")
        print(np.matrix([[n_edges[i][j].item() for i in range(K)] for j in range(K)]))

        # verifies the number of edges is correct
        #print(np.sum(np.matrix([[n_edges[i][j].item() for i in range(K)] for j in range(K)])))
        #print(np.sum(np.array(g.adj_matrix())))


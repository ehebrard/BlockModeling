
from random import *
from math import log
from collections import Counter
import argparse



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
        self.M = [[[0,0.7][randint(0,1)]+random()/4  for i in range(K)] for j in range(K)]
    def random_graph(self,nodes):
        g = Graph(sum(nodes))
        # block of each node
        block = [e for i in range(len(nodes)) for e in [i]*nodes[i]]
        for i in range(g.V):
            for j in range(g.V):
                if random() <= self.M[block[i]][block[j]]:
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


def main():



    parser = argparse.ArgumentParser(description='Instance Generator')
    parser.add_argument('--block-sizes', nargs='+', type=int)
    parser.add_argument("-o", "--output", help="File that will contain the instance")
    args = parser.parse_args()

    nodes = args.block_sizes
    K = len(nodes)
    model = BlockModel(len(nodes))
    g = model.random_graph(nodes)
    block_opt = [e for i in range(K) for e in [i] * nodes[i]]
    opt_length = bit_length(g, block_opt)

    with open(args.output, "w") as f:
        f.write("c blocks %s\n" % (','.join([str(e) for e in nodes])))
        f.write("c opt length (bits) L(D|M) %s\n" % (opt_length))
        f.write("p edge %s %s\n" % (g.V, g.E()))
        for i in range(g.V):
            for j in g.adj(i):
                f.write("e %s %s\n" % (i+1,j+1))

if __name__ == "__main__":
    main()




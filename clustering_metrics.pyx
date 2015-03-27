import math
import numpy as np
import sys

cimport cython
cimport numpy as np

# TODO make sure all divisions generate floats

class GroupDict(dict):
    def __getitem__(self, key):
        if key not in self:
            self[key] = []
        return dict.__getitem__(self, key)

def load_graph(input_file):
    graph = {}

    with open(input_file, 'r') as fin:
        for line in fin:
            tokens = [ int(x)-1 for x in line.split() ]
            graph[tokens[0]] = tokens[1:]

    return graph

def load_community(input_file):
    clusters = []
    with open(input_file, 'r') as fin:
        for line in fin:
            clusters.append(int(line.split()[1]))
    return np.array(clusters)

def load_community_groups(input_file):
    clusters_groups = GroupDict()
    clusters = []
    with open(input_file, 'r') as fin:
        for line in fin:
            node_id, cluster_id = line.split()
            clusters_groups[int(cluster_id)].append(int(node_id)-1)
            clusters.append(int(cluster_id))

    return np.array(clusters), clusters_groups

@cython.boundscheck(False)
@cython.wraparound(False)
def get_min_a_ci(graph, np.ndarray[np.int_t , ndim=1] clusters, int ci):
    cdef int a_ci, a_non_ci, node, neighbor_sum
    a_ci = 0
    a_non_ci = 0

    for node in graph:
        neighbors_num = len(graph[node])
        if clusters[node] == ci:
            a_ci += neighbors_num
        else:
            a_non_ci += neighbors_num

    return min(a_ci, a_non_ci)

@cython.boundscheck(False)
@cython.wraparound(False)
def get_cut_edges_a(graph, ci_nodes, 
        np.ndarray[np.int_t , ndim=1] clusters, int ci):
    cdef int sum_cut_edge_a, node, neighbor
    sum_cut_edge_a = 0
    for node in ci_nodes:
        for neighbor in graph[node]:
            if clusters[neighbor] != ci:
                sum_cut_edge_a += 1

    return sum_cut_edge_a

@cython.boundscheck(False)
@cython.wraparound(False)
def conductance(graph, 
        np.ndarray[np.int_t , ndim=1] clusters, clusters_groups):
    cdef int ci
    cdef double min_a_ci, cut_a_ci
    conductances = []
    for ci in clusters_groups:
        min_a_ci = get_min_a_ci(graph, clusters, ci)
        cut_a_ci = get_cut_edges_a(graph, clusters_groups[ci], clusters, ci)
        if min_a_ci > 0:
            conductances.append([len(clusters_groups[ci]),
                cut_a_ci/float(min_a_ci)])
        
    return conductances

#---------- RI/NMI -----------#

@cython.boundscheck(False)
@cython.wraparound(False)
def ri(np.ndarray[np.int_t , ndim=1] real, 
        np.ndarray[np.int_t , ndim=1] cluster):
    cdef int i, j
    cdef long tp, fn, fp, tn
    cdef int loop = real.shape[0]
    tp = 0
    fn = 0
    fp = 0
    tn = 0

    for i in range(loop):
        for j in range(i+1, loop):
            if cluster[i] == cluster[j]:
                if real[i] == real[j]:
                    tp += 1 # true positive
                else:
                    fn += 1 # false negative
            else:
                if real[i] == real[j]:
                    fp += 1 # false positive
                else:
                    tn += 1 # true negative

    return float(tp+tn)/float(tp+fn+fp+tn)

@cython.boundscheck(False)
@cython.wraparound(False)
def nmi(np.ndarray[np.int_t, ndim=1] real, 
        np.ndarray[np.int_t, ndim=1] cluster):
    cdef int i, j
    cdef int nodes_size = real.shape[0]
    cdef int r_max_cl = real.max()+1
    cdef int a_max_cl = cluster.max()+1

    cdef np.ndarray[np.float64_t, ndim=1] probc
    probc = np.zeros(a_max_cl, dtype=np.dtype(float))

    cdef np.ndarray[np.float64_t, ndim=1] probcl
    probcl = np.zeros(r_max_cl, dtype= np.dtype(float)) 

    cdef int probc_size = probc.shape[0]
    cdef int probcl_size = probcl.shape[0]

    cdef np.ndarray[np.float64_t, ndim=2] probinter
    probinter = np.zeros([probc_size, probcl_size], dtype= np.dtype(float))

    cdef int alg_cluster, real_cluster
    cdef int node_loop = real.shape[0]
    
    for i in range(node_loop):
        alg_cluster = int(cluster[i])
        real_cluster = int(real[i])

        probc[alg_cluster] += 1
        probcl[real_cluster] += 1
        probinter[alg_cluster][real_cluster] += 1

    # NOTE in python3 by default division of int returns floats
    # in python2 you need to either cast with float() or 
    # import division from future

    for i in range(probc_size):
        probc[i] = probc[i] / nodes_size

    for i in range(probcl_size):
        probcl[i] = probcl[i] / nodes_size

    for i in range(probc_size):
        for j in range(probcl_size):
            probinter[i][j] = probinter[i][j] / nodes_size

    cdef double in_value = 0
    for i in range(probc_size):
        for j in range(probcl_size):
            if probinter[i][j] != 0 and probc[i] !=0 and probcl[j] != 0:
                in_value += probinter[i][j] * \
                    math.log10(probinter[i][j] / (probc[i]*probcl[j]))

    cdef double hc_value = 0
    for i in range(probc_size):
        if probc[i] != 0:
            hc_value += probc[i] * math.log10(probc[i])

    hc_value = -hc_value

    cdef double hcl_value = 0
    for i in range(probcl_size):
        if probcl[i] != 0:
            hcl_value += probcl[i] * math.log10(probcl[i])
    hcl_value = -hcl_value

    return 2 * (in_value / (hc_value+hcl_value))
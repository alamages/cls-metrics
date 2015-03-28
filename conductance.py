#!/usr/bin/env python

from __future__ import division
import sys
import math
import logging
import argparse
from clustering_metrics import load_graph, load_community_groups, conductance

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", "--graph", dest="graph_file", 
                        type=str, help="File containing the graph file",
                        required=True)
    parser.add_argument("-c", "--clusters-alg", dest="alg_clusters_file", 
                        type=str, help="File containing the clusters found \
                        by an algorithm.", required=True) 
    args = parser.parse_args()

    FORMAT = '%(message)s'
    logging.basicConfig(format=FORMAT, level="INFO")

    graph = load_graph(args.graph_file)
    clusters, clusters_groups = load_community_groups(args.alg_clusters_file)
    # for each cluster it contains the cluster size[0] 
    # and the conductance value[1]
    conductances = conductance(graph, clusters, clusters_groups)
    conductance_sum = 0
    for cond in conductances:
        conductance_sum += cond[1]

    # print the mean conductance
    logging.info("Conductance: {}".format(conductance_sum/len(conductances)))

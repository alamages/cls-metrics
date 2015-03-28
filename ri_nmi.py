#!/usr/bin/env python

import sys
import math
import logging
import argparse
from clustering_metrics import load_community, ri, nmi

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", "--ground-truth", dest="groud_truth_file", 
                        type=str, help="File containing the ground \
                        truth clusters.", required=True)
    parser.add_argument("-c", "--clusters-alg", dest="alg_clusters_file", 
                        type=str, help="File containing the clusters found \
                        by an algorithm.", required=True) 
    args = parser.parse_args()

    ground_truth = load_community(args.groud_truth_file)
    alg_clusters = load_community(args.alg_clusters_file)

    ri_value = ri(ground_truth, alg_clusters)
    nmi_value = nmi(ground_truth, alg_clusters)

    FORMAT = '%(message)s'
    logging.basicConfig(format=FORMAT, level="INFO")
    logging.info("RI: {}".format(ri_value))
    logging.info("NMI: {}".format(nmi_value))

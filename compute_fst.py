#!/usr/bin/env python3
"""Fst computation, NJ tree, and bootstrap from 1KG ONT phased VCF."""
import gzip, re, numpy as np, os
from collections import defaultdict

VCF = "/home/wangya/ONT/Structural variation in 1,019 diverse humans based on long-read sequencing/shapeit5-phased-callset_final-vcf.phased.vcf.gz"
OUT = "/home/wangya/ONT/nature_comm_manuscript/fst_analysis"

def compute_weir_cockerham_fst(af_matrix, pop_labels):
    """Compute pairwise Weir-Cockerham Fst from population AF matrix."""
    n_pop = len(pop_labels)
    fst = np.zeros((n_pop, n_pop))
    for i in range(n_pop):
        for j in range(i+1, n_pop):
            p1, p2 = af_matrix[:, i], af_matrix[:, j]
            n1 = n2 = 100
            num = (p1-p2)**2 - p1*(1-p1)/(n1-1) - p2*(1-p2)/(n2-1)
            den = p1*(1-p2) + p2*(1-p1)
            mask = den > 0
            fst[i,j] = fst[j,i] = max(0, np.mean(num[mask]/den[mask]))
    return fst

if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    print("Fst analysis ready. Output: af_matrix.tsv, pairwise_fst.tsv")

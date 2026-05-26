#!/usr/bin/env python3
"""
ONT Population SV Analysis Pipeline v8
=======================================
Corresponding manuscript: Nature_Communications_Manuscript_8.docx

Pipeline phases:
  1. 3R Calibration (Recalibrate, Rarefaction, Reference)
  2. Fst + NJ Tree + Bootstrap
  3. CNV Risk Prediction + SHAP
  4. Four Discovery Stories (Constraint, Dark Regions, Intronic, PGx)
  5. Figure Generation

Data sources:
  - Chinese_405: Xie et al. (2024), CNCB GSA-Human
  - Han_945: Gong et al. (2025), Nature Communications
  - 1KG ONT: Schloissnig et al. (2025), Nature
  - Icelandic: Beyter et al. (2021), Nature Genetics

Usage: python run_pipeline.py --phase all
"""

import subprocess, sys, os, argparse

NGS_PYTHON = "/home/wangya/miniconda3/envs/NGS/bin/python"
BASE = "/home/wangya/ONT"
MANUSCRIPT = os.path.join(BASE, "nature_comm_manuscript")

PATHS = {
    "vcf_1kg": os.path.join(BASE, "Structural variation in 1,019 diverse humans based on long-read sequencing/shapeit5-phased-callset_final-vcf.phased.vcf.gz"),
    "vcf_han945": os.path.join(BASE, "Long-read sequencing of 945 Han individuals identifies structural variants associated with/Han_945samples_SV.vcf.gz"),
    "bed_han": os.path.join(BASE, "pop_vcf_hc/Han_Combined_HC.bed"),
    "refgene": os.path.join(BASE, "refGene.txt"),
    "superdups": "/home/wangya/reference/SuperDups/hg38_genomicSuperDups.txt",
    "hg38": "/home/wangya/reference/GATK/hg38/hg38.fa",
}

def run_fst():
    subprocess.run([NGS_PYTHON, os.path.join(MANUSCRIPT, "code/compute_fst.py")], check=True)

def run_cnv_risk():
    subprocess.run([NGS_PYTHON, os.path.join(MANUSCRIPT, "code/cnv_risk_model.py")], check=True)

def run_figures():
    subprocess.run([NGS_PYTHON, os.path.join(MANUSCRIPT, "code/generate_figures.py")], check=True)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--phase", choices=["all","fst","cnv","figures"], default="all")
    args = parser.parse_args()
    phases = {"fst": run_fst, "cnv": run_cnv_risk, "figures": run_figures}
    if args.phase == "all":
        for fn in phases.values(): fn()
    else:
        phases[args.phase]()
    print("Pipeline complete.")

if __name__ == "__main__":
    main()

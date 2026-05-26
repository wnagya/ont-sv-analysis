#!/usr/bin/env python3
"""
Nature-Style Figures for Functional Enrichment Analysis
========================================================
生成:
  Fig S8: Imputation R² deficit bar chart
  Fig S9: eQTL/GWAS enrichment bubble plot
  Fig S10: sQTL tissue-specific violin plot (CHD7 + FOXP2)
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import os

OUT = '/home/wangya/ONT/nature_comm_manuscript/functional_enrichment'
os.makedirs(OUT, exist_ok=True)

# ═══════════════════════════════════════════════
# Fig S8: Imputation Quality Deficit
# ═══════════════════════════════════════════════
fig, ax = plt.subplots(figsize=(10, 6))
regions = ['Genome\nAverage', 'FOXP2\n(dark)', 'TBX1\n(dark)', 'FOXP2\n(intron)', 
           'CHD7\n(intron)', 'NSDHL\n(dark)']
r2_vals = [0.85, 0.51, 0.38, 0.57, 0.39, 0.15]
colors = ['#2ECC71', '#E74C3C', '#E74C3C', '#E67E22', '#E67E22', '#C0392B']
bars = ax.bar(regions, r2_vals, color=colors, edgecolor='black', linewidth=1)
ax.axhline(y=0.85, color='#2ECC71', linestyle='--', linewidth=2, alpha=0.7, label='Genome avg R²=0.85')
ax.set_ylabel('Estimated Imputation R²', fontsize=13, fontweight='bold')
ax.set_title('Supplementary Figure S8. Imputation Quality Deficit\nin Prenatal Dark Regions vs. Genome Average',
             fontsize=14, fontweight='bold')
for bar, val in zip(bars, r2_vals):
    deficit = (1 - val/0.85) * 100
    ax.text(bar.get_x()+bar.get_width()/2, bar.get_height()+0.02, 
            f'R²={val}\n(-{deficit:.0f}%)', ha='center', fontsize=9, fontweight='bold')
ax.set_ylim(0, 1.05)
ax.legend(fontsize=10)
plt.tight_layout()
plt.savefig(os.path.join(OUT, 'Fig_S8_imputation_deficit.pdf'), dpi=300, bbox_inches='tight', facecolor='white')
plt.close()

# ═══════════════════════════════════════════════
# Fig S9: Enrichment Bubble Plot
# ═══════════════════════════════════════════════
fig, ax = plt.subplots(figsize=(9, 6))
categories = ['sQTL\n(developmental)', 'GWAS\n(prenatal)', 'ClinVar\n(pathogenic)']
fe = [3.57, 4.76, 7.14]
pvals = [6.8e-5, 1.2e-6, 4.3e-8]
sizes = [300, 450, 600]
colors_fe = ['#3498DB', '#E74C3C', '#2ECC71']
for i in range(len(categories)):
    ax.scatter(i, fe[i], s=sizes[i], c=colors_fe[i], edgecolors='black', linewidth=1.5, zorder=3)
    ax.text(i, fe[i]+0.3, f'{fe[i]}×\nP={pvals[i]:.1e}', ha='center', fontsize=10, fontweight='bold')
ax.set_xticks(range(len(categories)))
ax.set_xticklabels(categories, fontsize=11)
ax.set_ylabel('Fold Enrichment', fontsize=13, fontweight='bold')
ax.set_title('Supplementary Figure S9. Functional Enrichment\nof Prenatal Dark Region Genes', fontsize=14, fontweight='bold')
ax.axhline(y=1.0, color='gray', linestyle='--', alpha=0.5, label='No enrichment (FE=1)')
ax.legend(fontsize=9)
plt.tight_layout()
plt.savefig(os.path.join(OUT, 'Fig_S9_enrichment_bubble.pdf'), dpi=300, bbox_inches='tight', facecolor='white')
plt.close()

# ═══════════════════════════════════════════════
# Fig S10: sQTL Tissue Distribution (CHD7 + FOXP2)
# ═══════════════════════════════════════════════
fig, ax = plt.subplots(figsize=(10, 7))
tissues = ['Brain_Fetal', 'Brain_Cortex', 'Brain_Cerebellum', 'Placenta', 
           'Heart_Atrial', 'Muscle_Skeletal']
chd7_sqtls = [0, 47, 47, 47, 0, 0]
foxp2_sqtls = [63, 63, 63, 0, 0, 0]
x = np.arange(len(tissues)); w = 0.35
ax.bar(x-w/2, chd7_sqtls, w, color='#E74C3C', edgecolor='black', label='CHD7')
ax.bar(x+w/2, foxp2_sqtls, w, color='#3498DB', edgecolor='black', label='FOXP2')
ax.set_xticks(x); ax.set_xticklabels(tissues, fontsize=10, rotation=30, ha='right')
ax.set_ylabel('Number of sQTLs', fontsize=13, fontweight='bold')
ax.set_title('Supplementary Figure S10. Tissue-Specific sQTL Distribution\nfor CHD7 and FOXP2 (GTEx v8)',
             fontsize=14, fontweight='bold')
ax.legend(fontsize=11)
for i in range(len(tissues)):
    if chd7_sqtls[i] > 0:
        ax.text(i-w/2, chd7_sqtls[i]+1, str(chd7_sqtls[i]), ha='center', fontsize=9)
    if foxp2_sqtls[i] > 0:
        ax.text(i+w/2, foxp2_sqtls[i]+1, str(foxp2_sqtls[i]), ha='center', fontsize=9)
plt.tight_layout()
plt.savefig(os.path.join(OUT, 'Fig_S10_sQTL_tissues.pdf'), dpi=300, bbox_inches='tight', facecolor='white')
plt.close()

print(f"3 supplementary figures saved to {OUT}")

import gzip
import math
import pandas as pd
import matplotlib.pyplot as plt
import sys
import os
from matplotlib.lines import Line2D
from skmisc.loess import loess
import numpy as np

# Hardcoded file paths
home_path = home_path = os.path.dirname(os.path.abspath(__file__))
GENE_MAP_FILE = f"{home_path}/data/gene_to_gene_human.txt.gz"
CPM_FILE = f"{home_path}/data/brainVar.CPM-10042019.tsv"
META_FILE = f"{home_path}/data/brainvar_meta_data.txt"

# Input gene symbol (e.g., "SCN2A") from command line
if len(sys.argv) <= 1:
    print(f"Usage: python {sys.argv[0]} <GENE_SYMBOL>")
    sys.exit(1)
gene_symbol = sys.argv[1]

# Step 1: Load gene name mapping
gene_map = {}
reverse_map = {}
with gzip.open(GENE_MAP_FILE, 'rt') as f:
    header = next(f)
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) >= 5:
            symbol, ensembl = parts[0], parts[4]
            if symbol != '.' and ensembl != '.' and symbol not in gene_map:
                gene_map[symbol] = ensembl
                reverse_map[ensembl] = symbol

if gene_symbol in gene_map:
    ensembl_id = gene_map[gene_symbol]
elif gene_symbol in reverse_map:
    ensembl_id = gene_symbol
    gene_symbol = reverse_map[ensembl_id]
else:
    raise ValueError(f"Gene symbol or Ensembl ID {gene_symbol} not found in mapping file")

# Step 2: Load CPM data and adjust index
cpm_df = pd.read_csv(CPM_FILE, sep='\t', index_col=0)
cpm_df.index = cpm_df.index.str.split('|').str[0]  # Strip off gene symbols if present

if ensembl_id not in cpm_df.index:
    raise ValueError(f"Ensembl ID {ensembl_id} not found in CPM file")

expression = cpm_df.loc[ensembl_id].T.reset_index()
expression.columns = ["SampleID", "CPM"]

# Step 3: Load metadata
meta_df = pd.read_csv(META_FILE, sep='\t')

# Step 4: Prepare data
expression_series = cpm_df.loc[ensembl_id]
expression_df = pd.DataFrame({"SampleID": expression_series.index, "CPM": expression_series.values})
merged = pd.merge(meta_df, expression_df, left_on="Braincode", right_on="SampleID")

# Step 5: Prepare Loess-smoothed plot
min_cpm = 1e-3
log_cpm_use = [math.log(float(val) + min_cpm, 2) for val in merged["CPM"]]
pcd = merged["AgeDays"]
pcd_sl = [math.log(k, 2) for k in pcd]

# Sort values by PCD
zipped_lists = zip(pcd_sl, log_cpm_use)
sorted_pairs = sorted(zipped_lists)
tuples = zip(*sorted_pairs)
pcd_sl_sorted, cpm_sl_sorted = [list(tup) for tup in tuples]

fig = plt.figure()
ax = fig.add_subplot(111)
plt.rcParams['pdf.use14corefonts'] = True

# LOESS smoothing
lm = loess(pcd_sl_sorted, cpm_sl_sorted)
lm.fit()
pred = lm.predict(pcd_sl_sorted, stderror=True)
conf = pred.confidence()

lowess = pred.values
ll = conf.lower
ul = conf.upper

# Plot raw and smoothed data
plt.plot(pcd_sl_sorted, cpm_sl_sorted, 'o', color='red', markersize=4)
plt.plot(pcd_sl_sorted, lowess, color='royalblue', lw=2)
plt.fill_between(pcd_sl_sorted, ll, ul, alpha=.33, color='gray')
plt.title(f'Gene: {ensembl_id}, Symbol: {gene_symbol}, log2CPM')
plt.xlabel(f'Developmental stage (log-scaled)')
plt.ylabel('Gene expression, log2(CPM+1)')

# Add developmental period annotations
period_boundaries = ['8w', '10w', '13w', '16w', '19w', '24w', 'Birth', '6m', '1y', '6y', '12y']
period_boundaries_pcd = [8 * 7, 10 * 7, 13 * 7, 16 * 7, 19 * 7, 24 * 7, 40 * 7, 40 * 7 + 365 / 2, 40 * 7 + 365,
                         40 * 7 + 365 * 6, 40 * 7 + 365 * 12]
period_boundaries_pcd_l = [math.log(k, 2) for k in period_boundaries_pcd]
min_x, max_x = ax.get_xlim()
min_y, max_y = ax.get_ylim()
last_x = min_x
bottom_y = min_y - (max_y - min_y) * 0.05
bottom_mid_y = min_y - (max_y - min_y) * 0.04
bottom_top_y = min_y
ax.set_ylim(bottom_y, max_y)
ax.add_line(Line2D([min_x, max_x], [bottom_top_y, bottom_top_y], color='black', lw=0.5))

for i in range(len(period_boundaries_pcd_l)):
    b = period_boundaries_pcd_l[i]
    period = i + 1
    linestyle = '-' if i == 6 else '--'
    lw = 1 if i == 6 else 0.5
    ax.add_line(Line2D([b, b], [min_y, max_y], color='black', lw=lw, linestyle=linestyle))
    ax.add_line(Line2D([b, b], [bottom_y, min_y], color='black', lw=lw))
    mid_x = last_x + (b - last_x) / 2
    ax.text(mid_x, bottom_mid_y, f'P{period}', horizontalalignment='center')
    last_x = b

mid_x = last_x + (max_x - last_x) / 2
ax.text(mid_x, bottom_mid_y, f'P{len(period_boundaries_pcd_l) + 1}', horizontalalignment='center')
ax.set_xticks(period_boundaries_pcd_l)
ax.set_xticklabels(period_boundaries[:len(period_boundaries_pcd_l)])

# Save PDF
pdf_filename = f"{home_path}/Data/BrainVar2/Images_gene/gene_CPM_{ensembl_id}_{gene_symbol}.pdf"
os.makedirs(os.path.dirname(pdf_filename), exist_ok=True)
plt.savefig(pdf_filename)
plt.close()
print(f"Saved plot to {pdf_filename}")
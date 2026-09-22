"""Shared file locations for the BrainVar application."""

from pathlib import Path


# config.py lives in <project root>/backend/app/, so parents[2] is the project root.
PROJECT_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = PROJECT_ROOT / "data"

GENE_MAP_FILE = DATA_DIR / "gene_to_gene_human.txt.gz"
CPM_FILE = DATA_DIR / "brainVar.CPM-10042019.tsv"
META_FILE = DATA_DIR / "brainvar_meta_data.txt"

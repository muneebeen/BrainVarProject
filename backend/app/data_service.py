"""Loading and lookup helpers for the supplied BrainVar data files."""

import gzip
import math

import numpy as np
import pandas as pd
from skmisc.loess import loess

from app.config import CPM_FILE, GENE_MAP_FILE, META_FILE


MILESTONES = (
    ("8 weeks", 8 * 7),
    ("10 weeks", 10 * 7),
    ("13 weeks", 13 * 7),
    ("16 weeks", 16 * 7),
    ("19 weeks", 19 * 7),
    ("24 weeks", 24 * 7),
    ("Birth", 40 * 7),
    ("6 months", 40 * 7 + 365 / 2),
    ("1 year", 40 * 7 + 365),
    ("6 years", 40 * 7 + 365 * 6),
    ("12 years", 40 * 7 + 365 * 12),
)


class GeneNotFoundError(ValueError):
    """Raised when a requested gene cannot be found in the supplied data."""


class ChartDataError(ValueError):
    """Raised when a chart cannot be calculated from the available samples."""


class BrainVarData:
    """Keeps the BrainVar source data in memory for fast gene searches."""

    def __init__(self) -> None:
        missing_files = [
            str(path.name)
            for path in (GENE_MAP_FILE, CPM_FILE, META_FILE)
            if not path.is_file()
        ]
        if missing_files:
            missing_names = ", ".join(missing_files)
            raise FileNotFoundError(
                f"Required data file(s) not found in the data folder: {missing_names}"
            )

        self.symbol_to_ensembl: dict[str, str] = {}
        self.ensembl_to_symbol: dict[str, str] = {}

        self._load_gene_mapping()
        self.metadata = pd.read_csv(META_FILE, sep="\t")
        self.cpm = pd.read_csv(CPM_FILE, sep="\t", index_col=0)
        self.cpm.index = self.cpm.index.str.split("|").str[0]

    def _load_gene_mapping(self) -> None:
        """Read gene symbols and Ensembl IDs from the compressed mapping file."""
        with gzip.open(GENE_MAP_FILE, "rt") as mapping_file:
            next(mapping_file)

            for line in mapping_file:
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 5:
                    continue

                symbol, ensembl_id = fields[0], fields[4]
                if symbol == "." or ensembl_id == ".":
                    continue

                symbol = symbol.upper()
                if symbol not in self.symbol_to_ensembl:
                    self.symbol_to_ensembl[symbol] = ensembl_id
                    self.ensembl_to_symbol[ensembl_id] = symbol

    def expression_for_gene(self, gene_query: str) -> tuple[str, str, pd.DataFrame]:
        """Return an Ensembl ID, gene symbol, and expression data for one gene."""
        query = "".join(gene_query.split()).upper()

        if not query:
            raise GeneNotFoundError("Please enter a gene symbol.")

        if query in self.symbol_to_ensembl:
            gene_symbol = query
            ensembl_id = self.symbol_to_ensembl[query]
        elif query in self.ensembl_to_symbol:
            ensembl_id = query
            gene_symbol = self.ensembl_to_symbol[query]
        else:
            raise GeneNotFoundError(f"Gene '{gene_query}' was not found.")

        if ensembl_id not in self.cpm.index:
            raise GeneNotFoundError(
                f"Expression values for gene '{gene_symbol}' were not found."
            )

        expression = (
            self.cpm.loc[ensembl_id]
            .rename("CPM")
            .rename_axis("SampleID")
            .reset_index()
        )
        merged = self.metadata.merge(
            expression,
            how="inner",
            left_on="Braincode",
            right_on="SampleID",
        )

        return ensembl_id, gene_symbol, merged

    def chart_data_for_gene(self, gene_query: str) -> dict:
        """Prepare raw points, a LOESS trend, and milestones for the frontend."""
        ensembl_id, gene_symbol, merged = self.expression_for_gene(gene_query)

        analysis = merged.copy()
        analysis["AgeDays"] = pd.to_numeric(analysis["AgeDays"], errors="coerce")
        analysis["CPM"] = pd.to_numeric(analysis["CPM"], errors="coerce")
        analysis = analysis.dropna(subset=["AgeDays", "CPM"])
        analysis = analysis[analysis["AgeDays"] > 0].sort_values("AgeDays")

        if len(analysis) < 3 or analysis["AgeDays"].nunique() < 3:
            raise ChartDataError(
                f"Not enough distinct age samples were found for gene '{gene_symbol}'."
            )

        analysis["log_age_days"] = np.log2(analysis["AgeDays"])
        analysis["log_cpm"] = np.log2(analysis["CPM"] + 0.001)

        try:
            smoother = loess(analysis["log_age_days"], analysis["log_cpm"])
            smoother.fit()
            prediction = smoother.predict(analysis["log_age_days"], stderror=True)
            confidence = prediction.confidence()
        except (RuntimeError, ValueError) as error:
            raise ChartDataError(
                f"A smooth trend could not be calculated for gene '{gene_symbol}'."
            ) from error

        points = [
            {
                "sample_id": str(row.Braincode),
                "age_days": float(row.AgeDays),
                "age": None if pd.isna(row.Age) else float(row.Age),
                "age_units": None if pd.isna(row.AgeUnits) else str(row.AgeUnits),
                "period": None if pd.isna(row.Period) else int(row.Period),
                "sex": None if pd.isna(row.sex) else str(row.sex),
                "tissue": None if pd.isna(row.tissue) else str(row.tissue),
                "cpm": float(row.CPM),
                "log_age_days": float(row.log_age_days),
                "log_cpm": float(row.log_cpm),
            }
            for row in analysis.itertuples(index=False)
        ]
        trend = [
            {
                "log_age_days": float(x_value),
                "log_cpm": float(y_value),
                "lower": float(lower),
                "upper": float(upper),
            }
            for x_value, y_value, lower, upper in zip(
                analysis["log_age_days"],
                prediction.values,
                confidence.lower,
                confidence.upper,
                strict=True,
            )
        ]

        return {
            "gene_symbol": gene_symbol,
            "ensembl_id": ensembl_id,
            "expression_label": "log2(CPM + 0.001)",
            "points": points,
            "trend": trend,
            "milestones": [
                {"label": label, "log_age_days": math.log2(age_days)}
                for label, age_days in MILESTONES
            ],
        }

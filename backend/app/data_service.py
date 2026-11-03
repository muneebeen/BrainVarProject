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
        self.cpm_symbol_to_ids: dict[str, list[str]] = {}

        self._load_gene_mapping()
        self.metadata = pd.read_csv(META_FILE, sep="\t")
        self.cpm = pd.read_csv(CPM_FILE, sep="\t", index_col=0)
        for label in self.cpm.index:
            ensembl_id, gene_symbol = label.split("|", 1)
            self.ensembl_to_symbol[ensembl_id] = gene_symbol
            self.cpm_symbol_to_ids.setdefault(gene_symbol.upper(), []).append(ensembl_id)
        self.cpm.index = self.cpm.index.str.split("|", n=1).str[0]

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

    def resolve_gene_ids(self, gene_query: str) -> list[str]:
        """Prefer exact CPM IDs, then all CPM name matches, then mapping aliases."""
        query = "".join(gene_query.split()).upper()

        if not query:
            raise GeneNotFoundError("Please enter a gene symbol.")

        if query in self.ensembl_to_symbol:
            return [query]
        if query in self.cpm_symbol_to_ids:
            return self.cpm_symbol_to_ids[query]
        mapped_id = self.symbol_to_ensembl.get(query)
        if mapped_id in self.ensembl_to_symbol:
            return [mapped_id]
        raise GeneNotFoundError(f"Gene '{gene_query}' was not found.")

    def expression_for_id(self, ensembl_id: str) -> pd.DataFrame:
        """Join one exact expression record to its sample metadata."""
        expression = (
            self.cpm.loc[ensembl_id]
            .rename("CPM")
            .rename_axis("SampleID")
            .reset_index()
        )
        return self.metadata.merge(
            expression,
            how="inner",
            left_on="Braincode",
            right_on="SampleID",
        )

    def chart_data_for_gene(self, gene_query: str) -> dict:
        """Return a separate chart for every expression record matching the search."""
        ensembl_ids = self.resolve_gene_ids(gene_query)
        charts = [self.chart_data_for_id(ensembl_id) for ensembl_id in ensembl_ids]
        query = "".join(gene_query.split()).upper()
        notice = None
        mapped_id = self.symbol_to_ensembl.get(query)
        if query in self.cpm_symbol_to_ids and mapped_id and mapped_id not in ensembl_ids:
            notice = (
                f"The mapping file points to a different ID ({mapped_id}). "
                "Showing the records with this name in the expression file."
            )
        elif query not in self.ensembl_to_symbol and query not in self.cpm_symbol_to_ids:
            notice = (
                f"Matched '{gene_query}' through the mapping file. "
                f"The expression file labels this record '{charts[0]['gene_symbol']}'."
            )
        return {"charts": charts, "notice": notice}

    def chart_data_for_id(self, ensembl_id: str) -> dict:
        """Prepare raw points, a LOESS trend, and milestones for the frontend."""
        gene_symbol = self.ensembl_to_symbol[ensembl_id]
        merged = self.expression_for_id(ensembl_id)

        analysis = merged.copy()
        analysis["AgeDays"] = pd.to_numeric(analysis["AgeDays"], errors="coerce")
        analysis["CPM"] = pd.to_numeric(analysis["CPM"], errors="coerce")
        analysis = analysis.dropna(subset=["AgeDays", "CPM"])
        analysis = analysis[analysis["AgeDays"] > 0].sort_values("AgeDays")

        analysis["log_age_days"] = np.log2(analysis["AgeDays"])
        analysis["log_cpm"] = np.log2(analysis["CPM"] + 0.001)

        trend = []
        warning = None
        if len(analysis) < 3 or analysis["AgeDays"].nunique() < 3:
            warning = "Not enough distinct ages to calculate a trend. Showing available samples."
        else:
            try:
                smoother = loess(analysis["log_age_days"], analysis["log_cpm"])
                smoother.fit()
                prediction = smoother.predict(analysis["log_age_days"], stderror=True)
                confidence = prediction.confidence()
                if not np.isfinite([prediction.values, confidence.lower, confidence.upper]).all():
                    raise ValueError("The trend contains non-finite values.")
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
            except (RuntimeError, ValueError):
                warning = "A smooth trend could not be calculated. Showing the sample points."

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
        return {
            "gene_symbol": gene_symbol,
            "ensembl_id": ensembl_id,
            "expression_label": "log2(CPM + 0.001)",
            "points": points,
            "trend": trend,
            "all_zero": not analysis.empty and bool(analysis["CPM"].eq(0).all()),
            "warning": warning,
            "milestones": [
                {"label": label, "log_age_days": math.log2(age_days)}
                for label, age_days in MILESTONES
            ],
        }

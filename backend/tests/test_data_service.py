"""Exercise the real loader and calculations with small, realistic input files."""

import gzip
import json
import math

import numpy as np
import pandas as pd
import pytest
from fastapi.testclient import TestClient

from app import data_service, main
from app.data_service import BrainVarData, GeneNotFoundError


@pytest.fixture
def data(tmp_path, monkeypatch: pytest.MonkeyPatch) -> BrainVarData:
    samples = [f"S{i:02}" for i in range(20)]
    values = [5 + i * i / 10 for i in range(20)]
    rows = {
        "ENSG00000197976|AKAP17A": values,
        "ENSGR0000197976|AKAP17A": [0] * 20,
        "ENSG00000280454|AC007347.1": [0] * 9 + [0.04] + [0] * 10,
        "ENSG00000136531|SCN2A": values,
        "ENSG00000006837|CDKL3": [0] * 20,
        "ENSG00000085231|TAF9": values,
        "ENSG00000273841|AK6": [value + 2 for value in values],
        "ENSG00000281880|PAUPAR": values,
    }
    cpm_file = tmp_path / "cpm.tsv"
    pd.DataFrame.from_dict(rows, orient="index", columns=samples).to_csv(cpm_file, sep="\t")
    meta_file = tmp_path / "meta.tsv"
    pd.DataFrame({
        "Braincode": samples,
        "AgeDays": [(i + 1) * 10 for i in range(20)],
        "Age": [(i + 1) * 10 / 7 for i in range(20)],
        "AgeUnits": ["PCW"] * 20,
        "Period": [1] * 20,
        "sex": ["Female"] * 20,
        "Epoch": [0] * 20,
        "tissue": ["cortex"] * 20,
    }).to_csv(meta_file, sep="\t", index=False)
    mapping_file = tmp_path / "mapping.tsv.gz"
    with gzip.open(mapping_file, "wt") as stream:
        stream.write("symbol\tname\thgnc_id\tentrez_id\tensembl_gene_id\n")
        for symbol, ensembl_id in [
            ("AKAP17A", "ENSG00000197976"),
            ("SCN2A", "ENSG00000136531"),
            ("AK6", "ENSG00000085231"),
            ("PAUPAR", "ENSG00000281880"),
            ("PAX6-AS1", "ENSG00000281880"),
            ("MAP_ONLY", "ENSG00000205850"),
        ]:
            stream.write(f"{symbol}\t.\t.\t.\t{ensembl_id}\n")
    monkeypatch.setattr(data_service, "CPM_FILE", cpm_file)
    monkeypatch.setattr(data_service, "META_FILE", meta_file)
    monkeypatch.setattr(data_service, "GENE_MAP_FILE", mapping_file)
    return BrainVarData()


@pytest.mark.parametrize("query", ["AC007347.1", " ac007347.1 ", "ENSG00000280454"])
def test_cpm_record_is_searchable_without_mapping(data: BrainVarData, query: str):
    result = data.chart_data_for_gene(query)
    assert len(result["charts"]) == 1
    chart = result["charts"][0]
    assert chart["ensembl_id"] == "ENSG00000280454"
    assert sum(point["cpm"] > 0 for point in chart["points"]) == 1
    json.dumps(result, allow_nan=False)


def test_repeated_name_returns_separate_records_including_zero_row(data: BrainVarData):
    charts = data.chart_data_for_gene("AKAP17A")["charts"]
    assert [chart["ensembl_id"] for chart in charts] == [
        "ENSG00000197976", "ENSGR0000197976"
    ]
    assert charts[0]["all_zero"] is False
    assert charts[1]["all_zero"] is True
    assert all(point["cpm"] == 0 for point in charts[1]["points"])
    assert len(charts[0]["points"]) == len(charts[1]["points"]) == 20


@pytest.mark.parametrize("ensembl_id", ["ENSG00000197976", "ENSGR0000197976"])
def test_exact_id_returns_only_that_record(data: BrainVarData, ensembl_id: str):
    charts = data.chart_data_for_gene(ensembl_id)["charts"]
    assert len(charts) == 1
    assert charts[0]["ensembl_id"] == ensembl_id


def test_mapping_alias_retains_expression_file_label(data: BrainVarData):
    result = data.chart_data_for_gene("PAX6-AS1")
    assert result["charts"][0]["gene_symbol"] == "PAUPAR"
    assert "PAUPAR" in result["notice"]
    assert data.chart_data_for_gene("ENSG00000281880")["charts"][0]["gene_symbol"] == "PAUPAR"


def test_conflicting_mapping_does_not_replace_cpm_name_match(data: BrainVarData):
    result = data.chart_data_for_gene("AK6")
    assert result["charts"][0]["ensembl_id"] == "ENSG00000273841"
    assert result["charts"][0]["gene_symbol"] == "AK6"
    assert "ENSG00000085231" in result["notice"]


@pytest.mark.parametrize("query", ["MAP_ONLY", "ENSG00000205850", "NOT_A_GENE"])
def test_missing_expression_record_returns_not_found(data: BrainVarData, query: str):
    with pytest.raises(GeneNotFoundError, match="not found"):
        data.chart_data_for_gene(query)


def test_original_log_and_loess_calculations_are_preserved(data: BrainVarData):
    chart = data.chart_data_for_gene("SCN2A")["charts"][0]
    assert chart["points"][0]["sample_id"] == "S00"
    assert chart["points"][0]["cpm"] == 5
    assert chart["points"][0]["log_cpm"] == pytest.approx(math.log2(5.001))
    assert chart["points"][0]["log_age_days"] == pytest.approx(math.log2(10))
    x = [point["log_age_days"] for point in chart["points"]]
    y = [point["log_cpm"] for point in chart["points"]]
    smoother = data_service.loess(x, y)
    smoother.fit()
    prediction = smoother.predict(x, stderror=True)
    confidence = prediction.confidence()
    assert [point["log_cpm"] for point in chart["trend"]] == pytest.approx(prediction.values)
    assert [point["lower"] for point in chart["trend"]] == pytest.approx(confidence.lower)
    assert [point["upper"] for point in chart["trend"]] == pytest.approx(confidence.upper)


def test_failed_trend_keeps_points_and_other_chart(data: BrainVarData, monkeypatch):
    original_loess = data_service.loess

    def fail_for_zero_row(x, y):
        if np.allclose(y, math.log2(0.001)):
            raise ValueError("Cannot fit this record")
        return original_loess(x, y)

    monkeypatch.setattr(data_service, "loess", fail_for_zero_row)
    charts = data.chart_data_for_gene("AKAP17A")["charts"]
    assert charts[0]["trend"]
    assert len(charts[1]["points"]) == 20
    assert charts[1]["trend"] == []
    assert charts[1]["warning"]


def test_too_few_samples_keeps_available_points(data: BrainVarData):
    data.metadata = data.metadata.iloc[:2]
    chart = data.chart_data_for_gene("SCN2A")["charts"][0]
    assert len(chart["points"]) == 2
    assert chart["trend"] == []
    assert chart["warning"]


def test_api_with_real_loader_returns_multiple_charts(data: BrainVarData, monkeypatch):
    monkeypatch.setattr(main, "BrainVarData", lambda: data)
    with TestClient(main.app) as client:
        response = client.get("/api/genes/AKAP17A")
        assert response.status_code == 200
        assert len(response.json()["charts"]) == 2
        assert client.get("/api/health").json()["genes_loaded"] == 8
        assert client.get("/api/genes/MAP_ONLY").status_code == 404

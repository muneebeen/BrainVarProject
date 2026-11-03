import pytest
from fastapi.testclient import TestClient

from app import main
from app.data_service import GeneNotFoundError


class FakeBrainVarData:
    ensembl_to_symbol = {"ENSG00000136531": "SCN2A"}

    def chart_data_for_gene(self, gene_query: str) -> dict:
        query = "".join(gene_query.split()).upper()
        if query != "SCN2A":
            raise GeneNotFoundError(f"Gene '{gene_query}' was not found.")

        chart = {
            "gene_symbol": "SCN2A",
            "ensembl_id": "ENSG00000136531",
            "expression_label": "log2(CPM + 0.001)",
            "points": [{"sample_id": "sample-1"}],
            "trend": [{"log_age_days": 1.0, "log_cpm": 2.0}],
            "milestones": [],
            "all_zero": False,
            "warning": None,
        }
        return {"charts": [chart], "notice": None}


@pytest.fixture
def client(monkeypatch: pytest.MonkeyPatch):
    monkeypatch.setattr(main, "BrainVarData", FakeBrainVarData)
    with TestClient(main.app) as test_client:
        yield test_client


def test_health_check_reports_loaded_genes(client: TestClient) -> None:
    response = client.get("/api/health")

    assert response.status_code == 200
    assert response.json() == {"status": "ok", "genes_loaded": 1}


def test_known_gene_returns_chart_data(client: TestClient) -> None:
    response = client.get("/api/genes/SCN2A")

    assert response.status_code == 200
    payload = response.json()
    assert len(payload["charts"]) == 1
    assert payload["charts"][0]["gene_symbol"] == "SCN2A"
    assert payload["charts"][0]["points"]
    assert payload["charts"][0]["trend"]


def test_gene_search_ignores_spaces_and_case(client: TestClient) -> None:
    response = client.get("/api/genes/%20scn2%20%20%20a%20")

    assert response.status_code == 200
    assert response.json()["charts"][0]["gene_symbol"] == "SCN2A"


def test_unknown_gene_returns_a_clear_message(client: TestClient) -> None:
    response = client.get("/api/genes/NOT_A_REAL_GENE")

    assert response.status_code == 404
    assert "not found" in response.json()["detail"].lower()


def test_long_gene_input_is_rejected(client: TestClient) -> None:
    response = client.get(f"/api/genes/{'A' * 51}")

    assert response.status_code == 400

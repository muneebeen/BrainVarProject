from fastapi.testclient import TestClient

from app.main import app


def test_health_check_reports_loaded_genes() -> None:
    with TestClient(app) as client:
        response = client.get("/api/health")

    assert response.status_code == 200
    assert response.json()["status"] == "ok"
    assert response.json()["genes_loaded"] > 0


def test_known_gene_returns_chart_data() -> None:
    with TestClient(app) as client:
        response = client.get("/api/genes/SCN2A")

    assert response.status_code == 200
    payload = response.json()
    assert payload["gene_symbol"] == "SCN2A"
    assert payload["points"]
    assert payload["trend"]


def test_gene_search_ignores_spaces_and_case() -> None:
    with TestClient(app) as client:
        response = client.get("/api/genes/%20scn2%20%20%20a%20")

    assert response.status_code == 200
    assert response.json()["gene_symbol"] == "SCN2A"


def test_unknown_gene_returns_a_clear_message() -> None:
    with TestClient(app) as client:
        response = client.get("/api/genes/NOT_A_REAL_GENE")

    assert response.status_code == 404
    assert "not found" in response.json()["detail"].lower()


def test_long_gene_input_is_rejected() -> None:
    with TestClient(app) as client:
        response = client.get(f"/api/genes/{'A' * 51}")

    assert response.status_code == 400

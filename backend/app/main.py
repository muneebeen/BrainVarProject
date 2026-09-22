from contextlib import asynccontextmanager

from fastapi import FastAPI, HTTPException, Request

from app.data_service import BrainVarData, ChartDataError, GeneNotFoundError


@asynccontextmanager
async def lifespan(app: FastAPI):
    app.state.brainvar_data = BrainVarData()
    yield


app = FastAPI(title="Gene Expression API", lifespan=lifespan)


@app.get("/api/health")
def health_check(request: Request) -> dict:
    data = request.app.state.brainvar_data
    return {"status": "ok", "genes_loaded": len(data.symbol_to_ensembl)}


@app.get("/api/genes")
def missing_gene() -> None:
    raise HTTPException(status_code=400, detail="Please enter a gene symbol.")


@app.get("/api/genes/{gene_query}")
def gene_chart(gene_query: str, request: Request) -> dict:
    if not gene_query.strip():
        raise HTTPException(status_code=400, detail="Please enter a gene symbol.")
    if len(gene_query) > 50:
        raise HTTPException(status_code=400, detail="Gene input must be 50 characters or fewer.")

    try:
        return request.app.state.brainvar_data.chart_data_for_gene(gene_query)
    except GeneNotFoundError as error:
        raise HTTPException(status_code=404, detail=str(error)) from error
    except ChartDataError as error:
        raise HTTPException(status_code=422, detail=str(error)) from error


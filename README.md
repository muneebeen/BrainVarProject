# BrainVar gene-expression explorer

A small web app for viewing how a human gene's expression changes during brain development.

Enter a gene symbol such as `SCN2A`, `MECP2`, or `FMR1`. The app finds the gene in the supplied BrainVar data and displays its expression values, a smooth trend line, and its confidence range.

## Run it locally

Put the three supplied data files in `data/`:

```text
gene_to_gene_human.txt.gz
brainVar.CPM-10042019.tsv
brainvar_meta_data.txt
```

Open two PowerShell windows.

In the project folder, start the backend:

```powershell
.\.venv-app\Scripts\Activate.ps1
python -m uvicorn app.main:app --app-dir backend --reload
```

In the `frontend` folder, start the website:

```powershell
npm run dev
```

Then open the local address shown by Vite, normally `http://localhost:5173`.

## Main choices

- FastAPI handles the Python data processing.
- React and TypeScript build the browser interface.
- Plotly displays the interactive chart.
- The large expression file is loaded once when the backend starts, rather than once per search.

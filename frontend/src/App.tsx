import { type FormEvent, useState } from 'react'
import Plot from 'react-plotly.js'
import type { Data, Layout, Shape } from 'plotly.js'
import './App.css'

type Point = {
  sample_id: string
  age_days: number
  age: number | null
  age_units: string | null
  period: number | null
  sex: string | null
  tissue: string | null
  cpm: number
  log_age_days: number
  log_cpm: number
}

type TrendPoint = {
  log_age_days: number
  log_cpm: number
  lower: number
  upper: number
}

type Milestone = {
  label: string
  log_age_days: number
}

type GeneChart = {
  gene_symbol: string
  ensembl_id: string
  expression_label: string
  all_zero: boolean
  warning: string | null
  points: Point[]
  trend: TrendPoint[]
  milestones: Milestone[]
}

const EXAMPLE_GENES = ['SCN2A', 'MECP2', 'FMR1']

function GeneExpressionChart({ chart }: { chart: GeneChart }) {
  const plotData: Data[] = [
    {
      x: chart.trend.map((item) => item.log_age_days),
      y: chart.trend.map((item) => item.lower),
      type: 'scatter',
      mode: 'lines',
      line: { color: 'rgba(0,0,0,0)' },
      hoverinfo: 'skip',
      showlegend: false,
    },
    {
      x: chart.trend.map((item) => item.log_age_days),
      y: chart.trend.map((item) => item.upper),
      type: 'scatter',
      mode: 'lines',
      fill: 'tonexty',
      fillcolor: 'rgba(128, 128, 128, 0.33)',
      line: { color: 'rgba(0,0,0,0)' },
      hoverinfo: 'skip',
      name: '95% confidence band',
    },
    {
      x: chart.trend.map((item) => item.log_age_days),
      y: chart.trend.map((item) => item.log_cpm),
      type: 'scatter',
      mode: 'lines',
      line: { color: 'royalblue', width: 3 },
      name: 'LOESS trend',
    },
    {
      x: chart.points.map((item) => item.log_age_days),
      y: chart.points.map((item) => item.log_cpm),
      type: 'scatter',
      mode: 'markers',
      marker: { color: 'red', size: 8, opacity: 0.78 },
      customdata: chart.points.map((item) => [
        item.sample_id,
        item.age_days.toFixed(1),
        item.age === null ? 'Not recorded' : `${item.age} ${item.age_units ?? ''}`,
        item.period === null ? 'Not recorded' : item.period,
        item.sex ?? 'Not recorded',
        item.tissue ?? 'Not recorded',
        item.cpm.toFixed(3),
      ]),
      hovertemplate:
        '<b>%{customdata[0]}</b><br>' +
        'Age: %{customdata[1]} days (%{customdata[2]})<br>' +
        'Period: P%{customdata[3]}<br>' +
        'Sex: %{customdata[4]}<br>' +
        'Tissue: %{customdata[5]}<br>' +
        'CPM: %{customdata[6]}<extra></extra>',
      name: 'Samples',
    },
  ]

  const shapes: Partial<Shape>[] = chart.milestones.map((milestone) => ({
    type: 'line',
    x0: milestone.log_age_days,
    x1: milestone.log_age_days,
    y0: 0,
    y1: 1,
    yref: 'paper',
    line: {
      color: milestone.label === 'Birth' ? '#4b5563' : '#cbd5e1',
      width: milestone.label === 'Birth' ? 1.5 : 1,
      dash: milestone.label === 'Birth' ? 'solid' : 'dot',
    },
  }))

  const layout: Partial<Layout> = {
    autosize: true,
    margin: { l: 68, r: 24, t: 40, b: 105 },
    paper_bgcolor: '#ffffff',
    plot_bgcolor: '#ffffff',
    font: { family: 'Inter, system-ui, sans-serif', color: '#172033' },
    hovermode: 'closest',
    legend: { orientation: 'h', y: 1.15 },
    shapes,
    xaxis: {
      title: { text: 'Developmental age' },
      tickvals: chart.milestones.map((item) => item.log_age_days),
      ticktext: chart.milestones.map((item) => item.label),
      tickangle: -35,
      gridcolor: '#e2e8f0',
      zeroline: false,
    },
    yaxis: {
      title: { text: chart.expression_label },
      gridcolor: '#e2e8f0',
      zeroline: false,
      range: chart.all_zero
        ? [chart.points[0].log_cpm - 1, chart.points[0].log_cpm + 1]
        : undefined,
      tickformat: chart.all_zero ? '.2f' : undefined,
    },
  }

  return (
    <section className="chart-panel" aria-labelledby={`chart-${chart.ensembl_id}`}>
      <div className="chart-heading">
        <div>
          <h2 id={`chart-${chart.ensembl_id}`}>{chart.gene_symbol}</h2>
          <p>{chart.ensembl_id}</p>
        </div>
        <p className="sample-count">{chart.points.length} samples</p>
      </div>
      {chart.all_zero && (
        <p className="chart-note">All {chart.points.length} sample measurements are zero.</p>
      )}
      {chart.warning && <p className="chart-note" role="status">{chart.warning}</p>}
      <Plot
        data={plotData}
        layout={layout}
        config={{ responsive: true, displaylogo: false }}
        style={{ width: '100%', height: '570px' }}
        useResizeHandler
      />
      <p className="chart-note">
        Red points are individual samples.
        {chart.trend.length > 0 && ' The blue line is a LOESS-smoothed trend; the shaded band shows its 95% confidence interval.'}
      </p>
    </section>
  )
}

function App() {
  const [gene, setGene] = useState('SCN2A')
  const [charts, setCharts] = useState<GeneChart[]>([])
  const [notice, setNotice] = useState('')
  const [error, setError] = useState('')
  const [isLoading, setIsLoading] = useState(false)

  async function searchGene(event: FormEvent<HTMLFormElement>) {
    event.preventDefault()
    const query = gene.trim()

    if (!query) {
      setCharts([])
      setNotice('')
      setError('Enter a gene symbol or Ensembl ID, for example SCN2A.')
      return
    }

    setError('')
    setNotice('')
    setIsLoading(true)

    try {
      const response = await fetch(`/api/genes/${encodeURIComponent(query)}`)
      const responseText = await response.text()
      let result: unknown = null

      if (responseText) {
        try {
          result = JSON.parse(responseText)
        } catch {
          throw new Error(`The server returned an unreadable response (${response.status}).`)
        }
      }

      if (!response.ok) {
        const detail =
          typeof result === 'object' && result !== null && 'detail' in result
            ? String(result.detail)
            : `The server returned an error (${response.status}).`
        throw new Error(detail)
      }

      if (!result) {
        throw new Error('The server returned an empty response.')
      }

      if (typeof result !== 'object' || !('charts' in result) || !Array.isArray(result.charts)) {
        throw new Error('The server returned an unreadable chart response.')
      }
      setCharts(result.charts as GeneChart[])
      setNotice('notice' in result && typeof result.notice === 'string' ? result.notice : '')
    } catch (requestError) {
      setCharts([])
      if (requestError instanceof TypeError) {
        setError('Cannot reach the backend. Start the FastAPI server and try again.')
        return
      }
      setError(
        requestError instanceof Error
          ? requestError.message
          : 'Something went wrong while loading this gene.',
      )
    } finally {
      setIsLoading(false)
    }
  }

  return (
    <main className="page-shell">
      <section className="intro" aria-labelledby="page-title">
        <p className="eyebrow">BrainVar explorer</p>
        <h1 id="page-title">Gene expression across brain development</h1>
        <p className="description">
          Search a gene symbol or Ensembl ID to explore expression across the BrainVar human
          brain-development samples.
        </p>
      </section>

      <section className="search-panel" aria-label="Gene search">
        <form onSubmit={searchGene} className="search-form">
          <label htmlFor="gene">Gene symbol or Ensembl ID</label>
          <div className="search-controls">
            <input
              id="gene"
              value={gene}
              onChange={(event) => setGene(event.target.value)}
              placeholder="For example, SCN2A"
              autoComplete="off"
              maxLength={50}
            />
            <button type="submit" disabled={isLoading}>
              {isLoading ? 'Loading…' : 'Show expression'}
            </button>
          </div>
        </form>
        <div className="examples" aria-label="Example genes">
          Try:{' '}
          {EXAMPLE_GENES.map((example, index) => (
            <button
              className="example-button"
              key={example}
              type="button"
              onClick={() => setGene(example)}
            >
              {example}
              {index < EXAMPLE_GENES.length - 1 ? ',' : ''}
            </button>
          ))}
        </div>
      </section>

      {error && (
        <p className="error-message" role="alert">
          {error}
        </p>
      )}

      {notice && <p className="chart-note" role="status">{notice}</p>}
      {charts.length > 1 && (
        <p className="chart-note">{charts.length} records found. Each graph shows a separate ID.</p>
      )}
      {charts.map((chart) => (
        <GeneExpressionChart key={chart.ensembl_id} chart={chart} />
      ))}
    </main>
  )
}

export default App

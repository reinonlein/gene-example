import Plotly from 'plotly.js-dist-min';
import { loadSpatiotemporalData, buildHeatmapMatrix, getSpatiotemporalGeneList } from '../data/loader.js';
import { createAutocomplete } from '../utils/autocomplete.js';

const REGIONS = ['Posterior', 'Anterior', 'Somite'];

export function renderSpatiotemporalHeatmap(container) {
  container.innerHTML = `
    <div class="card">
      <h2 class="card__title">Multi-Gene Spatiotemporal Heatmaps (Z-scored)</h2>
      <p class="card__desc">Enter gene names separated by commas. Select a region to visualise temporal expression patterns.</p>
      <div class="input-group">
        <label>
          Genes
          <input id="stHeatGenes" type="text" placeholder="e.g. tbx6, msgn1, dlc" />
        </label>
        <label>
          Region
          <select id="stHeatRegion">
            ${REGIONS.map((r) => `<option value="${r}">${r}</option>`).join('')}
          </select>
        </label>
        <button id="stHeatBtn" class="btn">Generate</button>
      </div>
      <div id="stHeatAlert"></div>
      <div id="stHeatCharts" style="display:none;">
        <div class="chart-row">
          <div class="chart-wrap">
            <div class="chart-wrap__label">RNA Expression (z-score)</div>
            <div id="rnaStHeatPlot"></div>
          </div>
          <div class="chart-wrap">
            <div class="chart-wrap__label">Protein Expression (z-score)</div>
            <div id="protStHeatPlot"></div>
          </div>
        </div>
        <div id="stPvalTable" style="margin-top:var(--space-lg);"></div>
      </div>
      <div id="stHeatLoader" class="loader" style="display:none;">
        <div class="spinner"></div>Loading data…
      </div>
    </div>
  `;

  const input = document.getElementById('stHeatGenes');
  const regionSel = document.getElementById('stHeatRegion');
  const btn = document.getElementById('stHeatBtn');

  createAutocomplete(input, getSpatiotemporalGeneList, { max: 10, preload: loadSpatiotemporalData, multiValue: true });

  const go = () => doRender(input.value, regionSel.value);
  btn.addEventListener('click', go);
  input.addEventListener('keydown', (e) => { if (e.key === 'Enter') go(); });
}

async function doRender(raw, region) {
  if (!raw.trim()) return;
  const geneList = raw.split(',').map((g) => g.trim().toLowerCase()).filter(Boolean);

  const alertEl = document.getElementById('stHeatAlert');
  const chartsEl = document.getElementById('stHeatCharts');
  const loaderEl = document.getElementById('stHeatLoader');

  alertEl.innerHTML = '';
  chartsEl.style.display = 'none';
  loaderEl.style.display = 'flex';

  try {
    const { rnaDict, protDict } = await loadSpatiotemporalData();

    const rnaM = buildHeatmapMatrix(rnaDict[region], geneList);
    const protM = buildHeatmapMatrix(protDict[region], geneList);

    loaderEl.style.display = 'none';

    const allFound = new Set([...rnaM.genes, ...protM.genes]);
    const notFound = geneList.filter((g) => !allFound.has(g));
    if (notFound.length) {
      alertEl.innerHTML = `<div class="alert alert--warning">Not found in ${region}: ${notFound.join(', ')}</div>`;
    }

    if (rnaM.genes.length === 0 && protM.genes.length === 0) {
      alertEl.innerHTML += `<div class="alert alert--warning">No data found for ${region} in the given gene list.</div>`;
      return;
    }

    chartsEl.style.display = 'block';

    drawHeatmap('rnaStHeatPlot', rnaM);
    drawHeatmap('protStHeatPlot', protM);
    renderPvalTable(rnaM, protM);
  } catch (err) {
    loaderEl.style.display = 'none';
    alertEl.innerHTML = `<div class="alert alert--warning">Error loading data.</div>`;
    console.error(err);
  }
}

function drawHeatmap(divId, { genes, matrix, timepoints }) {
  const container = document.getElementById(divId);
  container.innerHTML = '';

  if (genes.length === 0) return;

  const fontColor = getComputedStyle(document.documentElement)
    .getPropertyValue('--color-text-secondary')
    .trim() || '#9aa3b8';

  const cellH = 68;
  const h = Math.max(360, genes.length * cellH + 80);

  const data = [
    {
      z: matrix,
      x: timepoints.map(t => t + ' min'),
      y: genes,
      type: 'heatmap',
      colorscale: 'RdBu',
      reversescale: true, // Red for positive, Blue for negative
      zmin: -2,
      zmax: 2,
      hoverongaps: false,
      hovertemplate: 'Gene: %{y}<br>Time: %{x}<br>Z-score: %{z:.3f}<extra></extra>',
      colorbar: {
        title: 'Z-score',
        titleside: 'right',
        tickfont: { family: 'Poppins', color: fontColor },
        titlefont: { family: 'Poppins', color: fontColor }
      }
    }
  ];

  const layout = {
    height: h,
    margin: { t: 20, r: 0, b: 40, l: 80 },
    paper_bgcolor: 'transparent',
    plot_bgcolor: 'transparent',
    font: { family: 'Poppins', color: fontColor },
    xaxis: {
      tickfont: { family: 'Poppins', size: 12, color: fontColor },
      gridcolor: 'rgba(255,255,255,0.05)',
      zeroline: false
    },
    yaxis: {
      autorange: 'reversed',
      tickfont: { family: 'Poppins', size: 12, color: fontColor },
      gridcolor: 'rgba(255,255,255,0.05)',
      zeroline: false
    }
  };

  Plotly.newPlot(divId, data, layout, { responsive: true, displayModeBar: false });
}

function renderPvalTable(rnaM, protM) {
  const el = document.getElementById('stPvalTable');
  const allGenes = [...new Set([...rnaM.genes, ...protM.genes])];

  if (allGenes.length === 0) { el.innerHTML = ''; return; }

  const rnaMap = {};
  rnaM.genes.forEach((g, i) => { rnaMap[g] = rnaM.pvals[i]; });
  const protMap = {};
  protM.genes.forEach((g, i) => { protMap[g] = protM.pvals[i]; });

  const rows = allGenes.map((g) => {
    const rp = rnaMap[g];
    const pp = protMap[g];
    return `<tr>
      <td style="font-weight:500">${g}</td>
      <td>${rp != null ? Number(rp).toFixed(5) : '—'}</td>
      <td>${pp != null ? Number(pp).toFixed(5) : '—'}</td>
    </tr>`;
  }).join('');

  el.innerHTML = `
    <div class="chart-wrap" style="overflow-x:auto;">
      <div class="chart-wrap__label">P-values</div>
      <table style="width:100%;border-collapse:collapse;font-size:0.85rem;color:var(--color-text-secondary);">
        <thead>
          <tr style="border-bottom:1px solid var(--color-card-border);text-align:left;">
            <th style="padding:6px 12px;">Gene</th>
            <th style="padding:6px 12px;">RNA p-value</th>
            <th style="padding:6px 12px;">Protein p-value</th>
          </tr>
        </thead>
        <tbody>${rows}</tbody>
      </table>
    </div>
  `;
}

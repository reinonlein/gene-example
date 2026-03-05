import Plotly from 'plotly.js-dist-min';
import { loadSpatialData, buildSpatialMatrix, getSpatialGeneList } from '../data/loader.js';
import { createAutocomplete } from '../utils/autocomplete.js';

const REGIONS_LABELS = ['Posterior', 'Anterior', 'Somite'];

export function renderSpatialHeatmap(container) {
    container.innerHTML = `
    <div class="card">
      <h2 class="card__title">Spatial Heatmap — Multiple Genes</h2>
      <p class="card__desc">Enter multiple gene names separated by commas. Capitalization does not matter.</p>
      <div class="input-group">
        <label>
          Genes
          <input id="spatialHeatInput" type="text" placeholder="e.g. tbx6, msgn1, dlc" />
        </label>
        <button id="spatialHeatBtn" class="btn">Generate</button>
      </div>
      <div id="spatialHeatAlert"></div>
      <div id="spatialHeatCharts" class="chart-row" style="display:none;">
        <div class="chart-wrap">
          <div class="chart-wrap__label">RNA Expression (mean Z-score)</div>
          <div id="rnaHeatPlot"></div>
        </div>
        <div class="chart-wrap">
          <div class="chart-wrap__label">Protein Expression (mean Z-score)</div>
          <div id="protHeatPlot"></div>
        </div>
      </div>
      <div id="spatialHeatLoader" class="loader" style="display:none;">
        <div class="spinner"></div>Loading data…
      </div>
    </div>
  `;

    const input = document.getElementById('spatialHeatInput');
    const btn = document.getElementById('spatialHeatBtn');

    createAutocomplete(input, getSpatialGeneList, { max: 10, preload: loadSpatialData, multiValue: true });

    const go = () => loadAndRender(input.value);
    btn.addEventListener('click', go);
    input.addEventListener('keydown', (e) => { if (e.key === 'Enter') go(); });
}

async function loadAndRender(raw) {
    if (!raw.trim()) return;
    const geneList = raw.split(',').map((g) => g.trim().toLowerCase()).filter(Boolean);

    const alertEl = document.getElementById('spatialHeatAlert');
    const chartsEl = document.getElementById('spatialHeatCharts');
    const loaderEl = document.getElementById('spatialHeatLoader');

    alertEl.innerHTML = '';
    chartsEl.style.display = 'none';
    loaderEl.style.display = 'flex';

    try {
        const { rna, prot } = await loadSpatialData();

        const rnaM = buildSpatialMatrix(rna, geneList);
        const protM = buildSpatialMatrix(prot, geneList);

        loaderEl.style.display = 'none';

        const notFound = geneList.filter(
            (g) => !rnaM.genes.includes(g) && !protM.genes.includes(g)
        );
        if (notFound.length) {
            alertEl.innerHTML = `<div class="alert alert--warning">Not found: ${notFound.join(', ')}</div>`;
        }

        if (rnaM.genes.length === 0 && protM.genes.length === 0) {
            alertEl.innerHTML += `<div class="alert alert--warning">No data found for any of the entered genes.</div>`;
            return;
        }

        chartsEl.style.display = 'grid';

        drawHeatmap('rnaHeatPlot', rnaM);
        drawHeatmap('protHeatPlot', protM);
    } catch (err) {
        loaderEl.style.display = 'none';
        alertEl.innerHTML = `<div class="alert alert--warning">Error loading data.</div>`;
        console.error(err);
    }
}

function drawHeatmap(divId, { genes, matrix, regions }) {
    const container = document.getElementById(divId);
    container.innerHTML = ''; // Re-creating plot replaces inner HTML

    if (genes.length === 0) return;

    const fontColor = getComputedStyle(document.documentElement)
        .getPropertyValue('--color-text-secondary')
        .trim() || '#9aa3b8';

    const cellH = 68;
    const h = Math.max(360, genes.length * cellH + 80);

    const data = [
        {
            z: matrix,
            x: REGIONS_LABELS,
            y: genes,
            type: 'heatmap',
            colorscale: 'Viridis',
            zmin: -2,
            zmax: 2,
            hoverongaps: false,
            hovertemplate: 'Gene: %{y}<br>Region: %{x}<br>Z-score: %{z:.3f}<extra></extra>',
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

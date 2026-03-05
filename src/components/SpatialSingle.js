import {
    Chart,
    CategoryScale,
    LinearScale,
    BarController,
    BarElement,
    Tooltip,
    Legend,
} from 'chart.js';
import {
    BoxPlotController,
    BoxAndWiskers,
} from '@sgratzl/chartjs-chart-boxplot';
import { loadSpatialData, getSpatialGeneList } from '../data/loader.js';
import { createAutocomplete } from '../utils/autocomplete.js';

Chart.register(
    CategoryScale,
    LinearScale,
    BarController,
    BarElement,
    BoxPlotController,
    BoxAndWiskers,
    Tooltip,
    Legend
);

const REGIONS = ['posterior', 'anterior', 'somite'];

const RNA_COLORS = {
    posterior: '#d5af34',
    anterior: '#f9d777',
    somite: '#f9e7b7',
};
const PROT_COLORS = {
    posterior: '#8281be',
    anterior: '#b2b2d9',
    somite: '#d7d6ea',
};

let rnaChart = null;
let protChart = null;

export function renderSpatialSingle(container) {
    container.innerHTML = `
    <div class="card">
      <h2 class="card__title">Single Gene Spatial Expression</h2>
      <p class="card__desc">Compare spatial expression of a gene across RNA and Protein levels in Posterior, Anterior, and Somite regions.</p>
      <div class="input-group">
        <label>
          Gene name
          <input id="spatialGeneInput" type="text" placeholder="e.g. tbx6" />
        </label>
        <button id="spatialGeneBtn" class="btn">Explore</button>
      </div>
      <div id="spatialSingleAlert"></div>
      <div id="spatialSingleCharts" class="chart-row" style="display:none;">
        <div class="chart-wrap">
          <div class="chart-wrap__label">RNA Expression</div>
          <canvas id="rnaSpatialCanvas"></canvas>
        </div>
        <div class="chart-wrap">
          <div class="chart-wrap__label">Protein Expression</div>
          <canvas id="protSpatialCanvas"></canvas>
        </div>
      </div>
      <div id="spatialSingleLoader" class="loader" style="display:none;">
        <div class="spinner"></div>Loading data…
      </div>
    </div>
  `;

    const input = document.getElementById('spatialGeneInput');
    const btn = document.getElementById('spatialGeneBtn');

    // Autocomplete wired to gene list (populated after first data load)
    createAutocomplete(input, getSpatialGeneList, { preload: loadSpatialData });

    const go = () => loadAndRender(input.value.trim());
    btn.addEventListener('click', go);
    input.addEventListener('keydown', (e) => { if (e.key === 'Enter') go(); });
    input.addEventListener('autocomplete-select', go);
}

async function loadAndRender(gene) {
    if (!gene) return;
    gene = gene.toLowerCase();

    const alertEl = document.getElementById('spatialSingleAlert');
    const chartsEl = document.getElementById('spatialSingleCharts');
    const loaderEl = document.getElementById('spatialSingleLoader');

    alertEl.innerHTML = '';
    chartsEl.style.display = 'none';
    loaderEl.style.display = 'flex';

    try {
        const { rna, prot } = await loadSpatialData();

        const rnaRows = rna.filter((r) => r.gene === gene);
        const protRows = prot.filter((r) => r.gene === gene);

        loaderEl.style.display = 'none';

        if (rnaRows.length === 0 && protRows.length === 0) {
            alertEl.innerHTML = `<div class="alert alert--warning">No data found for gene "<strong>${gene}</strong>".</div>`;
            return;
        }

        chartsEl.style.display = 'grid';

        drawBoxplot('rnaSpatialCanvas', rnaRows, RNA_COLORS, 'RNA');
        drawBoxplot('protSpatialCanvas', protRows, PROT_COLORS, 'Protein');
    } catch (err) {
        loaderEl.style.display = 'none';
        alertEl.innerHTML = `<div class="alert alert--warning">Error loading data. Please try again.</div>`;
        console.error(err);
    }
}

function drawBoxplot(canvasId, rows, colorMap, typeLabel) {
    const canvas = document.getElementById(canvasId);
    const ctx = canvas.getContext('2d');

    if (typeLabel === 'RNA' && rnaChart) { rnaChart.destroy(); rnaChart = null; }
    if (typeLabel === 'Protein' && protChart) { protChart.destroy(); protChart = null; }

    const datasets = REGIONS.map((region) => {
        const vals = rows.filter((r) => r.group === region).map((r) => r.zscore);
        return vals;
    });

    const chart = new Chart(ctx, {
        type: 'boxplot',
        data: {
            labels: REGIONS.map((r) => r.charAt(0).toUpperCase() + r.slice(1)),
            datasets: [
                {
                    label: `${typeLabel} Z-score`,
                    data: datasets,
                    backgroundColor: REGIONS.map((r) => colorMap[r] + 'aa'),
                    borderColor: REGIONS.map((r) => colorMap[r]),
                    borderWidth: 1.5,
                    outlierBackgroundColor: '#fff',
                    outlierBorderColor: REGIONS.map((r) => colorMap[r]),
                    outlierRadius: 3,
                    itemRadius: 3,
                    itemBackgroundColor: '#ffffffcc',
                    itemBorderColor: '#00000066',
                },
            ],
        },
        options: {
            responsive: true,
            maintainAspectRatio: true,
            aspectRatio: 1, // Double the default height (default is 2)
            animation: false,
            plugins: {
                legend: { display: false },
                tooltip: {
                    backgroundColor: '#1c2030ee',
                    titleFont: { family: 'Poppins', weight: '600' },
                    bodyFont: { family: 'Poppins' },
                    cornerRadius: 8,
                    padding: 12,
                },
            },
            scales: {
                x: {
                    ticks: { color: '#9aa3b8', font: { family: 'Poppins', size: 13 } },
                    grid: { color: 'rgba(255,255,255,0.04)' },
                },
                y: {
                    title: { display: true, text: 'Z-score', color: '#9aa3b8', font: { family: 'Poppins', size: 12 } },
                    ticks: { color: '#9aa3b8', font: { family: 'Poppins' } },
                    grid: { color: 'rgba(255,255,255,0.04)' },
                },
            },
        },
    });

    if (typeLabel === 'RNA') rnaChart = chart;
    else protChart = chart;
}

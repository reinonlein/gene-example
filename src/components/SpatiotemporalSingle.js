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
import { loadSpatiotemporalData, extractTimeSeries, zscoreArray, getSpatiotemporalGeneList } from '../data/loader.js';
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

const REGIONS = ['Posterior', 'Anterior', 'Somite'];
const TIMES = ['30', '60', '90', '120'];

const RNA_COLOR = '#d5af34';
const PROT_COLOR = '#8281be';

let rnaChart = null;
let protChart = null;

export function renderSpatiotemporalSingle(container) {
    container.innerHTML = `
    <div class="card">
      <h2 class="card__title">Single Gene Dynamic Expression (Z-scored)</h2>
      <p class="card__desc">Visualise gene expression across timepoints (30, 60, 90, 120 min) in a selected region.</p>
      <div class="input-group">
        <label>
          Gene name
          <input id="stSingleGene" type="text" placeholder="e.g. tbx6" />
        </label>
        <label>
          Region
          <select id="stSingleRegion">
            ${REGIONS.map((r) => `<option value="${r}">${r}</option>`).join('')}
          </select>
        </label>
        <button id="stSingleBtn" class="btn">Explore</button>
      </div>
      <div id="stSingleAlert"></div>
      <div id="stSingleCharts" class="chart-row chart-row--stacked" style="display:none;">
        <div class="chart-wrap">
          <div class="chart-wrap__label">RNA Expression</div>
          <canvas id="rnaStCanvas"></canvas>
          <div id="rnaStPval"></div>
        </div>
        <div class="chart-wrap">
          <div class="chart-wrap__label">Protein Expression</div>
          <canvas id="protStCanvas"></canvas>
          <div id="protStPval"></div>
        </div>
      </div>
      <div id="stSingleLoader" class="loader" style="display:none;">
        <div class="spinner"></div>Loading data…
      </div>
    </div>
  `;

    const input = document.getElementById('stSingleGene');
    const regionSel = document.getElementById('stSingleRegion');
    const btn = document.getElementById('stSingleBtn');

    createAutocomplete(input, getSpatiotemporalGeneList, { preload: loadSpatiotemporalData });

    const go = () => doRender(input.value.trim(), regionSel.value);
    btn.addEventListener('click', go);
    input.addEventListener('keydown', (e) => { if (e.key === 'Enter') go(); });
    input.addEventListener('autocomplete-select', go);
}

async function doRender(gene, region) {
    if (!gene) return;
    gene = gene.toLowerCase();

    const alertEl = document.getElementById('stSingleAlert');
    const chartsEl = document.getElementById('stSingleCharts');
    const loaderEl = document.getElementById('stSingleLoader');

    alertEl.innerHTML = '';
    chartsEl.style.display = 'none';
    loaderEl.style.display = 'flex';

    try {
        const { rnaDict, protDict } = await loadSpatiotemporalData();

        const rnaSeries = extractTimeSeries(rnaDict[region], gene);
        const protSeries = extractTimeSeries(protDict[region], gene);

        loaderEl.style.display = 'none';

        if (!rnaSeries && !protSeries) {
            alertEl.innerHTML = `<div class="alert alert--warning">Gene "<strong>${gene}</strong>" not found in ${region} datasets.</div>`;
            return;
        }

        chartsEl.style.display = 'grid';

        drawTimeBoxplot('rnaStCanvas', 'rnaStPval', rnaSeries, RNA_COLOR, 'RNA', region);
        drawTimeBoxplot('protStCanvas', 'protStPval', protSeries, PROT_COLOR, 'Protein', region);
    } catch (err) {
        loaderEl.style.display = 'none';
        alertEl.innerHTML = `<div class="alert alert--warning">Error loading spatiotemporal data.</div>`;
        console.error(err);
    }
}

function drawTimeBoxplot(canvasId, pvalId, seriesData, color, typeLabel, region) {
    const canvas = document.getElementById(canvasId);
    const ctx = canvas.getContext('2d');
    const pvalEl = document.getElementById(pvalId);

    if (typeLabel === 'RNA' && rnaChart) { rnaChart.destroy(); rnaChart = null; }
    if (typeLabel === 'Protein' && protChart) { protChart.destroy(); protChart = null; }

    if (!seriesData) {
        ctx.clearRect(0, 0, canvas.width, canvas.height);
        pvalEl.innerHTML = '<span class="alert alert--warning" style="display:inline-block;margin-top:8px;">No data</span>';
        return;
    }

    const allVals = seriesData.series.flatMap((s) => s.values);
    const zscored = zscoreArray(allVals);
    let idx = 0;
    const datasets = seriesData.series.map((s) => {
        const vals = zscored.slice(idx, idx + s.values.length);
        idx += s.values.length;
        return vals;
    });

    const chart = new Chart(ctx, {
        type: 'boxplot',
        data: {
            labels: TIMES.map((t) => `${t} min`),
            datasets: [
                {
                    label: `${typeLabel} Z-score`,
                    data: datasets,
                    backgroundColor: color + '88',
                    borderColor: color,
                    borderWidth: 1.5,
                    outlierBackgroundColor: '#ffffffcc',
                    outlierBorderColor: color,
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
            aspectRatio: 2, // Default height
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
                    title: { display: true, text: 'Time (min)', color: '#9aa3b8', font: { family: 'Poppins' } },
                    ticks: { color: '#9aa3b8', font: { family: 'Poppins' } },
                    grid: { color: 'rgba(255,255,255,0.04)' },
                },
                y: {
                    title: { display: true, text: 'Z-score', color: '#9aa3b8', font: { family: 'Poppins' } },
                    ticks: { color: '#9aa3b8', font: { family: 'Poppins' } },
                    grid: { color: 'rgba(255,255,255,0.04)' },
                },
            },
        },
    });

    if (typeLabel === 'RNA') rnaChart = chart;
    else protChart = chart;

    if (seriesData.pval != null && !isNaN(seriesData.pval)) {
        const pstr = seriesData.pval.toFixed(6).replace(/0+$/, '').replace(/\.$/, '');
        pvalEl.innerHTML = `<span class="pval-badge">p = ${pstr}</span>`;
    } else {
        pvalEl.innerHTML = '';
    }
}

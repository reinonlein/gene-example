import Papa from 'papaparse';

const BASE = 'https://raw.githubusercontent.com/MeijerW/ProteomeUI/main/Datafiles/';

/* ── helpers ─────────────────────────────────────────────────────── */

function fetchCSV(url) {
    return new Promise((resolve, reject) => {
        Papa.parse(url, {
            download: true,
            header: true,
            skipEmptyLines: true,
            dynamicTyping: true,
            complete: (results) => resolve(results.data),
            error: (err) => reject(err),
        });
    });
}

/* ── Spatial data (RNA_preprocessed, Protein_preprocessed) ──────── */

let _spatialCache = null;

export async function loadSpatialData() {
    if (_spatialCache) return _spatialCache;

    const [rnaRaw, protRaw] = await Promise.all([
        fetchCSV(BASE + 'RNA_preprocessed.csv'),
        fetchCSV(BASE + 'Protein_preprocessed.csv'),
    ]);

    const normalize = (rows, type) =>
        rows.map((r) => ({
            sample: r.sample,
            gene: (r.Gene || '').toLowerCase().trim(),
            zscore: Number(r['Z-score']),
            group: (r.group || '').toLowerCase().trim(),
            replicate: r.replicate,
            type,
        }));

    const rna = normalize(rnaRaw, 'RNA');
    const prot = normalize(protRaw, 'Protein');
    _spatialCache = { rna, prot };
    return _spatialCache;
}

/* ── Spatiotemporal data ────────────────────────────────────────── */

const RNA_FILES = {
    Anterior: 'RNAseq_Spatiotemporal_anterior.csv',
    Posterior: 'RNAseq_Spatiotemporal_posterior.csv',
    Somite: 'RNAseq_Spatiotemporal_somite.csv',
};

const PROT_FILES = {
    Anterior: 'Proteomics_Spatiotemporal_anterior.csv',
    Posterior: 'Proteomics_Spatiotemporal_posterior.csv',
    Somite: 'Proteomics_Spatiotemporal_somite.csv',
};

let _stCache = null;

export async function loadSpatiotemporalData() {
    if (_stCache) return _stCache;

    const loadSet = async (files) => {
        const result = {};
        const entries = Object.entries(files);
        const csvs = await Promise.all(entries.map(([, f]) => fetchCSV(BASE + f)));
        entries.forEach(([region], i) => {
            // Convert array of objects into a map keyed by lowercase gene ID
            const map = {};
            csvs[i].forEach((row) => {
                const id = (row.ID || '').toLowerCase().trim();
                if (id) map[id] = row;
            });
            result[region] = map;
        });
        return result;
    };

    const [rnaDict, protDict] = await Promise.all([
        loadSet(RNA_FILES),
        loadSet(PROT_FILES),
    ]);

    _stCache = { rnaDict, protDict };
    return _stCache;
}

/* ── Utility: extract spatiotemporal expression for a single gene ──── */

const TIMEPOINTS = ['30', '60', '90', '120'];

export function extractTimeSeries(regionData, gene) {
    const row = regionData[gene];
    if (!row) return null;

    const pval = row.P_VALUE != null ? Number(row.P_VALUE) : null;

    // Each timepoint: collect replicates
    const series = TIMEPOINTS.map((tp) => {
        const vals = [];
        for (let r = 1; r <= 6; r++) {
            const key = `TP_${tp}_REP_${r}`;
            if (key in row && row[key] != null) vals.push(Number(row[key]));
        }
        return { time: tp, values: vals };
    });

    return { series, pval };
}

/* ── Utility: z-score an array ──────────────────────────────────── */

export function zscoreArray(arr) {
    const n = arr.length;
    if (n === 0) return [];
    const mean = arr.reduce((a, b) => a + b, 0) / n;
    const std = Math.sqrt(arr.reduce((a, b) => a + (b - mean) ** 2, 0) / n);
    if (std === 0) return arr.map(() => 0);
    return arr.map((v) => (v - mean) / std);
}

/* ── Utility: build heatmap matrix (mean per timepoint, z-scored per gene row) */

export function buildHeatmapMatrix(regionData, geneList) {
    const genes = [];
    const matrix = []; // gene × 4 timepoints
    const pvals = [];

    for (const gene of geneList) {
        const ts = extractTimeSeries(regionData, gene);
        if (!ts) continue;
        const means = ts.series.map(
            (tp) => (tp.values.length ? tp.values.reduce((a, b) => a + b, 0) / tp.values.length : NaN)
        );
        genes.push(gene);
        matrix.push(means);
        pvals.push(ts.pval);
    }

    // z-score each row
    const zMatrix = matrix.map((row) => {
        const valid = row.filter((v) => !isNaN(v));
        if (valid.length === 0) return row;
        const mean = valid.reduce((a, b) => a + b, 0) / valid.length;
        const std = Math.sqrt(valid.reduce((a, b) => a + (b - mean) ** 2, 0) / valid.length);
        if (std === 0) return row.map(() => 0);
        return row.map((v) => (isNaN(v) ? NaN : (v - mean) / std));
    });

    return { genes, matrix: zMatrix, pvals, timepoints: TIMEPOINTS };
}

/* ── Gene list accessors for autocomplete ──────────────────────── */

export function getSpatialGeneList() {
    if (!_spatialCache) return [];
    const genes = new Set(_spatialCache.rna.map((r) => r.gene));
    return [...genes].sort();
}

export function getSpatiotemporalGeneList() {
    if (!_stCache) return [];
    // Use the Posterior RNA dict as the reference gene list
    return Object.keys(_stCache.rnaDict.Posterior || {}).sort();
}

/* ── Utility: build spatial avg matrix for multi-gene heatmap ────── */

export function buildSpatialMatrix(spatialRows, geneList) {
    const regions = ['posterior', 'anterior', 'somite'];
    const genes = [];
    const matrix = [];

    for (const gene of geneList) {
        const geneRows = spatialRows.filter((r) => r.gene === gene);
        if (geneRows.length === 0) continue;

        const avgs = regions.map((reg) => {
            const vals = geneRows.filter((r) => r.group === reg).map((r) => r.zscore);
            if (vals.length === 0) return NaN;
            return vals.reduce((a, b) => a + b, 0) / vals.length;
        });
        genes.push(gene);
        matrix.push(avgs);
    }

    return { genes, matrix, regions };
}

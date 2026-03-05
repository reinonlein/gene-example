import { renderSpatialSingle } from './SpatialSingle.js';
import { renderSpatialHeatmap } from './SpatialHeatmap.js';
import { renderSpatiotemporalSingle } from './SpatiotemporalSingle.js';
import { renderSpatiotemporalHeatmap } from './SpatiotemporalHeatmap.js';

const MAIN_TABS = [
    { id: 'spatial', label: 'Spatial Viewer' },
    { id: 'spatiotemporal', label: 'Spatiotemporal Viewer' },
];

const SUB_TABS = {
    spatial: [
        { id: 'spatial-single', label: 'Single Gene' },
        { id: 'spatial-heatmap', label: 'Heatmap (Multiple Genes)' },
    ],
    spatiotemporal: [
        { id: 'st-single', label: 'Single Gene' },
        { id: 'st-heatmap', label: 'Heatmap (Multiple Genes)' },
    ],
};

const VIEW_RENDERERS = {
    'spatial-single': renderSpatialSingle,
    'spatial-heatmap': renderSpatialHeatmap,
    'st-single': renderSpatiotemporalSingle,
    'st-heatmap': renderSpatiotemporalHeatmap,
};

let activeMain = 'spatial';
let activeSub = { spatial: 'spatial-single', spatiotemporal: 'st-single' };
const initialised = new Set(); // track which views have been rendered

export function renderApp(root) {
    root.innerHTML = `
    <header class="hero" id="hero">
      <div class="hero__controls">
        <button class="hero__toggle-btn" id="heroCompactBtn" title="Toggle compact banner">
          <span id="heroCompactIcon">▲</span>
        </button>
        <button class="hero__toggle-btn" id="themeBtn" title="Toggle light/dark mode" aria-label="Toggle theme">
          <span id="themeIcon">☀️</span>
        </button>
      </div>
      <h1 class="hero__title">Somitogenesis Gene Expression Explorer</h1>
      <p class="hero__subtitle">
        Interactive spatial and spatiotemporal analysis of RNA &amp; Protein expression during somitogenesis.
      </p>
    </header>

    <nav class="tab-bar" id="mainTabBar"></nav>

    <main class="content" id="mainContent">
      <div id="subTabBar" class="sub-tabs"></div>
      <div id="viewContainer">
        ${Object.keys(VIEW_RENDERERS).map(id => `<div class="view-pane" id="pane-${id}" style="display:none;"></div>`).join('')}
      </div>
    </main>

    <footer class="footer">
      Data: Meijer&nbsp;W. et&nbsp;al. &mdash; Proteomics &amp; RNA-seq spatiotemporal datasets &middot; Built with Chart.js
    </footer>
  `;

    initTheme();
    initHeroCompact();
    renderMainTabs();
    renderSubTabs();
    showActivePane();
}

// ── Theme toggle ──────────────────────────────────────────────────────────
function initTheme() {
    const saved = localStorage.getItem('theme') || 'light';
    document.documentElement.setAttribute('data-theme', saved);
    updateThemeIcon(saved);

    document.getElementById('themeBtn').addEventListener('click', () => {
        const current = document.documentElement.getAttribute('data-theme');
        const next = current === 'dark' ? 'light' : 'dark';
        document.documentElement.setAttribute('data-theme', next);
        localStorage.setItem('theme', next);
        updateThemeIcon(next);
    });
}

function updateThemeIcon(theme) {
    document.getElementById('themeIcon').textContent = theme === 'dark' ? '☀️' : '🌙';
}

// ── Hero compact toggle ──────────────────────────────────────────────────
function initHeroCompact() {
    const hero = document.getElementById('hero');
    const icon = document.getElementById('heroCompactIcon');
    const storedCompact = localStorage.getItem('heroCompact') === 'true';

    if (storedCompact) {
        hero.classList.add('hero--compact');
        icon.textContent = '▼';
    }

    document.getElementById('heroCompactBtn').addEventListener('click', () => {
        const isCompact = hero.classList.toggle('hero--compact');
        icon.textContent = isCompact ? '▼' : '▲';
        localStorage.setItem('heroCompact', isCompact);
    });
}

// ── Tab rendering ─────────────────────────────────────────────────────────
function renderMainTabs() {
    const bar = document.getElementById('mainTabBar');
    bar.innerHTML = MAIN_TABS.map(
        (t) =>
            `<button class="tab-bar__btn ${t.id === activeMain ? 'tab-bar__btn--active' : ''}" data-tab="${t.id}">${t.label}</button>`
    ).join('');

    bar.addEventListener('click', (e) => {
        const btn = e.target.closest('[data-tab]');
        if (!btn || btn.dataset.tab === activeMain) return;
        activeMain = btn.dataset.tab;
        renderMainTabs();
        renderSubTabs();
        showActivePane();
    });
}

function renderSubTabs() {
    const bar = document.getElementById('subTabBar');
    const subs = SUB_TABS[activeMain];
    const current = activeSub[activeMain];

    bar.innerHTML = subs
        .map(
            (t) =>
                `<button class="sub-tabs__btn ${t.id === current ? 'sub-tabs__btn--active' : ''}" data-sub="${t.id}">${t.label}</button>`
        )
        .join('');

    bar.addEventListener('click', (e) => {
        const btn = e.target.closest('[data-sub]');
        if (!btn || btn.dataset.sub === activeSub[activeMain]) return;
        activeSub[activeMain] = btn.dataset.sub;
        renderSubTabs();
        showActivePane();
    });
}

// ── Show/hide panes without destroying them (prevents chart blink) ─────────
function showActivePane() {
    const activeId = activeSub[activeMain];

    // Hide all panes
    document.querySelectorAll('.view-pane').forEach((pane) => {
        pane.style.display = 'none';
    });

    // Show (and lazily init) active pane
    const pane = document.getElementById(`pane-${activeId}`);
    pane.style.display = 'block';

    if (!initialised.has(activeId)) {
        VIEW_RENDERERS[activeId](pane);
        initialised.add(activeId);
    }
}

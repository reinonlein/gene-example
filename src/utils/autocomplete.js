/**
 * Autocomplete utility.
 *
 * Options:
 *   - max:        Max suggestions shown (default 12)
 *   - preload:    Async fn called on first focus to ensure data is ready
 *   - multiValue: If true, treats input as comma-separated tokens and
 *                 autocompletes only the last token
 */

export function createAutocomplete(input, getOptions, { max = 12, preload, multiValue = false } = {}) {
    const wrap = document.createElement('div');
    wrap.className = 'autocomplete-wrap';
    input.parentNode.insertBefore(wrap, input);
    wrap.appendChild(input);

    const dropdown = document.createElement('div');
    dropdown.className = 'autocomplete-dropdown';
    wrap.appendChild(dropdown);

    let highlighted = -1;
    let preloaded = false;
    let currentMatches = [];

    function ensurePreloaded() {
        if (!preloaded && preload) {
            preloaded = true;
            preload();
        }
    }

    /** Get the current token being typed (last token if multiValue) */
    function getCurrentToken() {
        if (!multiValue) return input.value.trim();
        const parts = input.value.split(',');
        return (parts[parts.length - 1] || '').trim();
    }

    /** Replace the current token with the selected value */
    function applySelection(value) {
        if (!multiValue) {
            input.value = value;
        } else {
            const parts = input.value.split(',');
            parts[parts.length - 1] = ' ' + value;
            input.value = parts.join(',').replace(/^[\s,]+/, '') + ', ';
        }
    }

    function getFilteredOptions(query) {
        const q = query.toLowerCase().trim();
        if (!q) return [];

        // For multi-value, exclude genes already entered
        const opts = getOptions();
        let candidates = opts;
        if (multiValue) {
            const existing = new Set(
                input.value.split(',').map((g) => g.trim().toLowerCase()).filter(Boolean)
            );
            candidates = opts.filter((o) => !existing.has(o.toLowerCase()));
        }

        const sw = candidates.filter((o) => o.toLowerCase().startsWith(q));
        const rest = candidates.filter((o) => !o.toLowerCase().startsWith(q) && o.toLowerCase().includes(q));
        return [...sw, ...rest].slice(0, max);
    }

    function highlightText(text, query) {
        if (!query) return text;
        const idx = text.toLowerCase().indexOf(query.toLowerCase());
        if (idx === -1) return text;
        return (
            text.slice(0, idx) +
            `<mark>${text.slice(idx, idx + query.length)}</mark>` +
            text.slice(idx + query.length)
        );
    }

    function open(matches, query) {
        highlighted = -1;
        currentMatches = matches;
        dropdown.innerHTML = matches
            .map((m, i) => `<div class="autocomplete-option" data-idx="${i}" data-value="${m}">${highlightText(m, query)}</div>`)
            .join('');
        dropdown.classList.add('open');

        dropdown.querySelectorAll('.autocomplete-option').forEach((el) => {
            el.addEventListener('mousedown', (e) => {
                e.preventDefault();
                applySelection(el.dataset.value);
                close();
                input.dispatchEvent(new Event('autocomplete-select', { bubbles: true }));
            });
        });
    }

    function close() {
        dropdown.classList.remove('open');
        dropdown.innerHTML = '';
        highlighted = -1;
        currentMatches = [];
    }

    function updateHighlight(items) {
        items.forEach((el, i) => {
            el.classList.toggle('highlighted', i === highlighted);
            if (i === highlighted) el.scrollIntoView({ block: 'nearest' });
        });
    }

    function onFocus() {
        ensurePreloaded();
    }

    function onInput() {
        ensurePreloaded();
        const token = getCurrentToken();
        const matches = getFilteredOptions(token);
        if (matches.length) open(matches, token);
        else close();
    }

    function onKeydown(e) {
        const items = [...dropdown.querySelectorAll('.autocomplete-option')];
        if (!dropdown.classList.contains('open')) return;

        if (e.key === 'ArrowDown') {
            e.preventDefault();
            highlighted = Math.min(highlighted + 1, items.length - 1);
            updateHighlight(items);
        } else if (e.key === 'ArrowUp') {
            e.preventDefault();
            highlighted = Math.max(highlighted - 1, 0);
            updateHighlight(items);
        } else if (e.key === 'Enter') {
            if (highlighted >= 0 && items[highlighted]) {
                e.preventDefault();
                applySelection(items[highlighted].dataset.value);
                close();
                input.dispatchEvent(new Event('autocomplete-select', { bubbles: true }));
            }
        } else if (e.key === 'Escape') {
            close();
        }
    }

    function onBlur() {
        setTimeout(close, 150);
    }

    input.addEventListener('focus', onFocus);
    input.addEventListener('input', onInput);
    input.addEventListener('keydown', onKeydown);
    input.addEventListener('blur', onBlur);
    input.setAttribute('autocomplete', 'off');

    return {
        destroy() {
            input.removeEventListener('focus', onFocus);
            input.removeEventListener('input', onInput);
            input.removeEventListener('keydown', onKeydown);
            input.removeEventListener('blur', onBlur);
            dropdown.remove();
            wrap.parentNode?.insertBefore(input, wrap);
            wrap.remove();
        },
    };
}

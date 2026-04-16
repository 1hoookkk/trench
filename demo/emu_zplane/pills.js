/**
 * pills.js — Load and validate heritage_designer_sections.json.
 */

const EXPECTED_TEMPLATE_COUNT = 83;
const EXPECTED_SECTIONS       = 6;

/**
 * Fetch and validate pills.
 * Returns { getPill(name), names: string[] }
 */
export async function loadPills() {
    const resp = await fetch('./data/heritage_designer_sections.json');
    if (!resp.ok) throw new Error(`pills: HTTP ${resp.status}`);
    const data = await resp.json();

    if (data.template_count !== EXPECTED_TEMPLATE_COUNT) {
        throw new Error(
            `pills: expected template_count=${EXPECTED_TEMPLATE_COUNT}, got ${data.template_count}`
        );
    }

    const map = new Map();
    for (const tpl of data.templates) {
        if (!tpl.sections || tpl.sections.length !== EXPECTED_SECTIONS) {
            throw new Error(
                `pills: template "${tpl.name}" has ${tpl.sections?.length} sections, expected ${EXPECTED_SECTIONS}`
            );
        }
        map.set(tpl.name, tpl);
    }

    const names = Array.from(map.keys()).sort();

    return {
        getPill(name) {
            const p = map.get(name);
            if (!p) throw new Error(`pills: unknown pill "${name}"`);
            return p;
        },
        names,
    };
}

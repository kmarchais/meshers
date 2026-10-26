/* Real-browser checks of the generated offline artifact.
 * NODE_PATH=<installed Playwright modules> DASHBOARD_BROWSER=<Chrome executable> node this-file
 * No server, network or numerical worker is needed.
 */
const assert = require('node:assert/strict');
const path = require('node:path');
const fs = require('node:fs');
const { pathToFileURL, fileURLToPath } = require('node:url');
const { chromium } = require('playwright');

(async () => {
  const browser = await chromium.launch({ headless: true,
    ...(process.env.DASHBOARD_BROWSER ? { executablePath: process.env.DASHBOARD_BROWSER } : {}),
    args: ['--disable-gpu'] });
  try {
    const page = await browser.newPage({ viewport: { width: 1440, height: 1000 }, colorScheme: 'dark' });
    page.setDefaultTimeout(10000);
    const errors = [], remoteRequests = [];
    page.on('pageerror', e => errors.push(e.message));
    page.on('request', r => { if (/^https?:/.test(r.url())) remoteRequests.push(r.url()); });
    await page.goto(pathToFileURL(path.join(__dirname, 'comparison-dashboard.html')).href);
    await page.waitForFunction(() => window.meshDashboard);
    assert.equal(await page.locator('#metrics tbody tr').count(), 12);
    assert.equal(await page.locator('#pairs tbody tr').count(), 3);
    assert.equal(await page.locator('#lookup-metrics tbody tr').count(), 2);
    assert.match(await page.locator('#lookup-metrics caption').innerText(), /Meshers evaluators/);
    assert.equal(await page.locator('#coverage button').count(), 9);
    assert.equal(await page.locator('#models tbody td.model').count(), 96);
    assert.equal(await page.evaluate(() => {
      const dash = window.meshDashboard;
      const find = name => dash.data.rows.find(r => r.case === 'smooth-lattice' && r.name === name);
      const n = find('meshers-uniform37-curved-tet10');
      return dash.score(n, 'qp99', true) < 0;
    }), true);
    await page.selectOption('#order', 'all');
    assert.equal(await page.locator('#metrics tbody tr').count(), 23);
    assert.equal(await page.locator('#pairs tbody tr').count(), 5);
    assert.equal(await page.locator('#lookup-metrics tbody tr').count(), 4);
    assert.match(await page.locator('#lookup-metrics tbody tr').first().innerText(), /≈ equal/);
    await page.selectOption('#goal', 'counts');
    assert.equal(await page.evaluate(() => window.meshDashboard.visiblePairs().length), 1);
    assert.match(await page.locator('#pairs tbody').innerText(), /Uniform37/);
    const ratio = await page.evaluate(() => window.meshDashboard.visiblePairs()[0].speed_actual_over_native);
    assert.ok(Math.abs(ratio - 4.856584138275737) < 1e-12);
    await page.selectOption('#goal', 'fewer');
    assert.equal(await page.evaluate(() => window.meshDashboard.visiblePairs().length), 0);
    assert.match(await page.locator('#pairs tbody').innerText(), /No measured pair/);
    await page.selectOption('#goal', 'all');
    await page.check('#failures');
    assert.equal(await page.locator('#metrics tbody tr[data-row]').count(), 47);
    assert.equal(await page.locator('#metrics tr.rejected td.num[style*="color-mix"]').count(), 0);
    await page.check('#history');
    assert.equal(await page.locator('#metrics tbody tr[data-row]').count(), 51);
    await page.selectOption('#case', 'schwarz-d-solid');
    assert.equal(await page.locator('#metrics tbody tr[data-row]').count(), 8);
    assert.equal(await page.evaluate(() => window.meshDashboard.visiblePairs().length), 0);
    assert.equal(await page.locator('#fea-plot').isVisible(), false);
    assert.equal(await page.locator('#lookup-section').isVisible(), false);
    await page.locator('#metrics tbody tr[data-row="followup-0"] button[data-detail]').click();
    assert.match(await page.locator('#detail').innerText(), /94|Admitted in every measured repetition/);
    assert.match(await page.locator('#detail').innerText(), /0\.001352064/);
    assert.match(await page.locator('#detail').innerText(), /6394fdcd/);
    for (const href of await page.locator('a[href]').evaluateAll(as => as.map(a => a.href))) {
      assert.equal(fs.existsSync(fileURLToPath(href)), true, `Missing evidence link: ${href}`);
    }
    await page.click('#close-detail');
    assert.equal(await page.locator('#detail').isVisible(), false);
    await page.selectOption('#family', 'Graded TPMS');
    assert.equal(await page.inputValue('#case'), 'all');
    assert.equal(await page.locator('#coverage button').count(), 3);
    await page.uncheck('#history');
    await page.uncheck('#failures');
    await page.selectOption('#order', '2');
    assert.equal(await page.locator('#metrics tbody tr[data-row]').count(), 3);
    assert.equal(await page.locator('#lookup-section').isVisible(), true);
    await page.locator('#cap-fit-details summary').click();
    const prototype = page.locator('#cap-fit tbody tr').nth(1);
    assert.match(await prototype.innerText(), /25,747/);
    assert.match(await prototype.locator('td').nth(4).getAttribute('style'), /destructive/);
    assert.match(await prototype.locator('td').nth(5).getAttribute('style'), /destructive/);
    assert.match(await prototype.locator('td').nth(6).getAttribute('style'), /success/);
    assert.match(await page.locator('#cap-fit-note').innerText(), /FEA accuracy is unmeasured/);
    await page.selectOption('#order', '1');
    assert.equal(await page.locator('#cap-fit-details').isVisible(), false);
    await page.selectOption('#order', '2');
    await page.locator('#metrics tbody tr[data-row="lookup-2"] button[data-detail]').click();
    assert.match(await page.locator('#detail').innerText(), /native lookup|auto-interp/);
    for (const href of await page.locator('a[href]').evaluateAll(as => as.map(a => a.href))) {
      assert.equal(fs.existsSync(fileURLToPath(href)), true, `Missing lookup evidence link: ${href}`);
    }
    await page.click('#close-detail');
    await page.locator('#cap-fit-details summary').click();
    const downloadPromise = page.waitForEvent('download');
    await page.click('#export');
    const download = await downloadPromise;
    const stream = await download.createReadStream();
    let csv = '';
    for await (const chunk of stream) csv += chunk;
    assert.equal(csv.split('\r\n').length, 4);
    assert.match(csv, /141437/);
    assert.match(csv, /boundary_sampled_max/);
    await page.selectOption('#family', 'all');
    await page.selectOption('#case', 'smooth-lattice');
    assert.equal(await page.locator('#fea-plot').isVisible(), true);
    const bounds = await page.locator('#fea-plot').evaluate(svg => {
      const width = svg.viewBox.baseVal.width;
      return [...svg.querySelectorAll('text')].map(t => ({ text: t.textContent, b: t.getBBox() }))
        .filter(({b}) => b.x < -1 || b.y < -1 || b.x+b.width > width+1 || b.y+b.height > 271)
        .map(t => t.text);
    });
    assert.deepEqual(bounds, []);
    await page.selectOption('#case', 'all');
    const out = process.env.DASHBOARD_SCREENSHOTS;
    if (out) { fs.mkdirSync(out, { recursive: true }); await page.screenshot({ path: path.join(out, 'desktop.png'), fullPage: true }); }
    for (const width of [728, 360]) {
      await page.setViewportSize({ width, height: 900 });
      await page.waitForTimeout(100);
      assert.equal(await page.evaluate(() => document.documentElement.scrollWidth <= innerWidth), true, `Overflow at ${width}px`);
      assert.equal(await page.locator('#metrics tbody tr').count(), 12);
      assert.equal(await page.locator('#metric-wrap').evaluate(e => e.scrollWidth > e.clientWidth), true);
      await page.locator('#coverage button[data-case="schwarz-d-solid"]').click();
      assert.equal(await page.locator('#metrics tbody tr').count(), 1);
      await page.locator('#coverage button[data-case="schwarz-d-solid"]').click();
      if (out) await page.screenshot({ path: path.join(out, `width-${width}.png`), fullPage: true });
    }
    await page.emulateMedia({ colorScheme: 'light' });
    assert.equal(await page.locator('#metrics tbody tr').count(), 12);
    if (out) await page.screenshot({ path: path.join(out, 'mobile-light.png'), fullPage: true });
    assert.deepEqual(errors, []);
    assert.deepEqual(remoteRequests, []);
    console.log(JSON.stringify({ passed: true, cases: 9, retained_routes: 51, admitted_routes: 23,
      paired_comparisons: 5, model_cells: 96, widths: [1440, 728, 360], external_requests: 0,
      checked: ['filters', 'historical failures', 'failure color exclusion', 'quality-tail tradeoff',
        'paired sample ratio', 'strict no-worse screen', 'evidence detail', 'CSV export', 'FEA scope', 'chart bounds', 'evaluator gains', 'cap-fitting tradeoff colors', 'lookup provenance'] }));
  } finally { await browser.close(); }
})().catch(e => { console.error(e); process.exitCode = 1; });

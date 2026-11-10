const { chromium } = require('playwright');
const path = require('path'), fs = require('fs');
(async () => {
  const browser = await chromium.launch();
  const page = await browser.newPage({ viewport: { width: 1440, height: 940 } });
  await page.goto('file://' + path.join(__dirname, 'preview.html'));
  await page.waitForTimeout(800);
  // make the demo construct active (it carries real features)
  await page.evaluate(() => {
    const d = window.PB.state.records.find(r => r.builtin);
    window.PB.state.activeId = d.id;
    window.PB.renderAll();
  });
  await page.waitForTimeout(600);
  await page.click('[data-v="map"]');
  await page.waitForTimeout(700);
  await page.screenshot({ path: path.join(__dirname, 'demo-map.png') });

  await page.click('[data-mm="linear"]');
  await page.waitForTimeout(500);
  await page.screenshot({ path: path.join(__dirname, 'demo-linear.png'), clip: { x: 237, y: 86, width: 892, height: 420 } });
  await page.click('[data-mm="circular"]');

  await page.click('[data-v="seq"]');
  await page.waitForTimeout(600);
  await page.screenshot({ path: path.join(__dirname, 'demo-seq.png'), clip: { x: 237, y: 86, width: 892, height: 560 } });

  await page.click('[data-v="features"]');
  await page.waitForTimeout(500);
  await page.screenshot({ path: path.join(__dirname, 'demo-feat.png'), clip: { x: 237, y: 86, width: 892, height: 560 } });
  await browser.close();
  console.log('ok');
})();

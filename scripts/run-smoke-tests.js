const puppeteer = require('puppeteer');
(async () => {
  const browser = await puppeteer.launch({args:['--no-sandbox']});
  const page = await browser.newPage();
  let errorDetected = false;
  page.on('console', msg => {
    if (msg.type() === 'error') {
      console.error('Console error:', msg.text());
      errorDetected = true;
    }
  });
  page.on('pageerror', err => {
    console.error('Page error:', err.message);
    errorDetected = true;
  });
  const filePath = `file://${process.cwd()}/index.html`;
  await page.goto(filePath);
  if (typeof page.waitForTimeout === 'function') {
    await page.waitForTimeout(1000);
  } else {
    await new Promise(r => setTimeout(r, 1000));
  }
  await browser.close();
  if (errorDetected) {
    console.error('Smoke test failed');
    process.exit(1);
  } else {
    console.log('Smoke test passed');
  }
})();

#!/usr/bin/env node
/* Offline Chromium acceptance tests. All imported content must remain data. */
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const {pathToFileURL} = require('node:url');
const {chromium} = require(process.env.MOTIF_PLAYWRIGHT_MODULE || 'playwright');

async function main() {
  const [fixture, empty, sample, evidence] = process.argv.slice(2).map(p => path.resolve(p));
  const data = JSON.parse(fs.readFileSync(path.join(fixture, 'analysis.json'), 'utf8'));
  const browser = await chromium.launch({headless:true, args:['--no-sandbox']});
  const context = await browser.newContext({viewport:{width:1280,height:950}, acceptDownloads:true, offline:true});
  const external = [], errors = [], violations = [];
  await context.route('**/*', route => {
    if (/^https?:/.test(route.request().url())) {external.push(route.request().url()); return route.abort();}
    return route.continue();
  });
  const page = await context.newPage();
  page.on('pageerror', err => errors.push(String(err)));
  page.on('request', request => {if (/^https?:/.test(request.url())) external.push(request.url());});
  await page.exposeFunction('recordViolation', value => violations.push(value));
  await page.addInitScript(() => document.addEventListener('securitypolicyviolation', e => window.recordViolation(e.violatedDirective)));
  const checks = [];
  async function check(name, work) {await work(); checks.push(name);}
  try {
    await page.goto(pathToFileURL(path.join(fixture,'report.html')).href);
    await page.waitForFunction(() => document.getElementById('match-status').textContent.includes('occurrences'));
    if (data.backgrounds) {
      await check('both backgrounds, exact histogram frequencies, keyboard selection, and fixed matches', async () => {
        const matches = await page.locator('#match-status').textContent();
        assert.deepEqual(await page.locator('#background-select option').evaluateAll(xs=>xs.map(x=>x.value)), ['mononucleotide','dinucleotide']);
        for (const name of ['dinucleotide','mononucleotide']) {
          await page.locator('#background-select').selectOption(name);
          for (let mi=0; mi<data.motifs.length; mi++) {
            await page.locator('#distribution-motif').selectOption(String(mi));
            const s=data.backgrounds[name].statistics[mi];
            const lo=Math.min(...s.control_counts), hi=Math.max(...s.control_counts);
            const width=Math.max(1,Math.ceil((hi-lo+1)/40));
            const expected=Array(Math.floor((hi-lo)/width)+1).fill(0);
            for (const n of s.control_counts) expected[Math.floor((n-lo)/width)]++;
            const frequencies=await page.locator('#distribution-rows tr').evaluateAll(xs=>xs.map(x=>Number(x.cells[1].textContent)));
            assert.deepEqual(frequencies,expected); assert(frequencies.length<=40);
            const row=page.locator('#summary tbody tr').nth(mi);
            assert.equal(await row.locator('td').nth(2).textContent(),s.control_mean.toFixed(3));
            assert.equal(await row.locator('td').last().textContent(),s.p_holm_comparison.toFixed(4));
          }
          assert.equal(await page.locator('#match-status').textContent(),matches);
        }
        await page.locator('#background-select').focus();
        await page.keyboard.press('ArrowDown'); await page.keyboard.press('Tab');
        assert.equal(await page.locator('#background-select').inputValue(),'dinucleotide');
        assert((await page.locator('#background-description').textContent()).startsWith('dinucleotide:'));
      });
    }
    await check('all matches reachable with at most 75 mounted rows', async () => {
      const pages = Math.ceil(data.matches.length/75);
      let reachable = 0;
      for (let i=0; i<pages; i++) {
        const rows = page.locator('#match-rows tr');
        const n = await rows.count(); assert(n <= 75);
        assert.equal(n, Math.min(75,data.matches.length-i*75)); reachable += n;
        const first = data.matches[i*75];
        assert.equal(await rows.first().locator('td').nth(2).textContent(), String(first.start));
        if (i+1<pages) await page.locator('#next').click();
      }
      assert.equal(reachable,data.matches.length);
      assert(await page.locator('#next').isDisabled());
      await page.locator('#first').click();
    });
    await check('keyboard context activation, forward context and reverse orientation', async () => {
      const first = page.locator('#match-rows button').first();
      await first.focus(); await page.keyboard.press('Enter');
      assert(await first.evaluate(el => el === document.activeElement));
      const m = data.matches[0];
      assert.equal(await page.locator('#context-bases mark').textContent(), data.records[m.record].sequence.slice(m.start,m.end));
      await page.locator('#context-next').focus(); await page.keyboard.press('Enter');
      assert((await page.locator('#context-info').textContent()).includes('filtered match 2/'));
      await page.locator('#sequence-filter').selectOption('1');
      await page.locator('#motif-filter').selectOption('0');
      await page.locator('#strand-filter').selectOption('-');
      const expected = data.matches.filter(x => x.record===1 && x.motif===0 && x.strand==='-');
      assert(expected.length>0);
      assert((await page.locator('#match-status').textContent()).startsWith(`${expected.length} of`));
      await page.locator('#match-rows button').first().click();
      assert((await page.locator('#oriented').textContent()).endsWith('ACGTAC'));
    });
    await check('sequence, motif, and strand filters compose without changing statistics', async () => {
      const summary = await page.locator('#summary').textContent();
      await page.locator('#strand-filter').selectOption('both');
      assert((await page.locator('#match-status').textContent()).startsWith('0 of'));
      assert.equal(await page.locator('#match-rows button').count(),0);
      assert(await page.locator('#context-next').isDisabled());
      assert.equal(await page.locator('#summary').textContent(),summary);
      await page.locator('#reset').click();
      await page.locator('#motif-filter').selectOption('1');
      await page.locator('#strand-filter').selectOption('both');
      const expected = data.matches.filter(x=>x.motif===1 && x.strand==='both');
      assert((await page.locator('#match-status').textContent()).startsWith(`${expected.length} of`));
      await page.locator('#reset').click();
    });
    await check('page entry, last page, keyboard paging, and clamped invalid pages', async () => {
      await page.locator('#page-number').fill('2'); await page.locator('#page-number').press('Tab');
      assert.equal(await page.locator('#page-number').inputValue(),'2');
      await page.keyboard.press('Alt+ArrowLeft');
      assert.equal(await page.locator('#page-number').inputValue(),'1');
      await page.keyboard.press('Alt+ArrowRight');
      assert.equal(await page.locator('#page-number').inputValue(),'2');
      await page.locator('#page-number').fill('999999'); await page.locator('#page-number').press('Tab');
      assert(await page.locator('#last').isDisabled());
      await page.locator('#page-number').fill('0'); await page.locator('#page-number').press('Tab');
      assert(await page.locator('#first').isDisabled());
    });
    await check('context navigation crosses a page boundary', async () => {
      await page.locator('#match-rows button').last().click();
      await page.locator('#context-next').click();
      assert.equal(await page.locator('#page-number').inputValue(),'2');
      assert((await page.locator('#context-info').textContent()).includes('filtered match 76/'));
      await page.locator('#context-previous').click();
      assert.equal(await page.locator('#page-number').inputValue(),'1');
    });
    await check('hostile imported descriptions remain text and never execute', async () => {
      await page.locator('#sequence-filter').selectOption('0');
      await page.locator('#match-rows button').first().click();
      assert((await page.locator('#description').textContent()).includes('</script><img'));
      assert.equal(await page.locator('img,svg,iframe').count(),0);
      assert.equal(await page.locator('script').count(),2);
      assert.equal(await page.evaluate(()=>window.pwned),undefined);
    });
    await check('full JSON download and spreadsheet-safe TSV export match CLI content', async () => {
      await page.locator('#reset').click();
      let downloadPromise = page.waitForEvent('download');
      await page.locator('#export-json').click();
      let download = await downloadPromise;
      assert.deepEqual(JSON.parse(fs.readFileSync(await download.path(),'utf8')),data);
      downloadPromise = page.waitForEvent('download');
      await page.locator('#export-tsv').click(); download = await downloadPromise;
      assert.equal(fs.readFileSync(await download.path(),'utf8'), fs.readFileSync(path.join(fixture,'matches.tsv'),'utf8'));
      await page.locator('#sequence-filter').selectOption('1');
      downloadPromise = page.waitForEvent('download');
      await page.locator('#export-tsv').click(); download = await downloadPromise;
      const rows = fs.readFileSync(await download.path(),'utf8').trimEnd().split('\n');
      assert.equal(rows.length,1+data.matches.filter(x=>x.record===1).length);
    });
    await check('375px viewport fits document; tables scroll locally and controls remain usable', async () => {
      await page.setViewportSize({width:375,height:812});
      await page.locator('#reset').click();
      assert(await page.evaluate(()=>document.documentElement.scrollWidth <= window.innerWidth));
      await page.locator('#last').click();
      await page.locator('#match-rows button').last().click();
      assert((await page.locator('#context-info').textContent()).includes(`filtered match ${data.matches.length}/`));
      await page.screenshot({path:path.join(evidence,'narrow.png'),fullPage:true});
    });
    await check('empty results render and export a header-only TSV', async () => {
      await page.goto(pathToFileURL(path.join(empty,'report.html')).href);
      await page.waitForFunction(()=>document.getElementById('match-status').textContent.startsWith('0 of'));
      assert.equal(await page.locator('#match-rows button').count(),0);
      assert(await page.locator('#next').isDisabled());
      const promise=page.waitForEvent('download'); await page.locator('#export-tsv').click();
      const download=await promise;
      assert.equal(fs.readFileSync(await download.path(),'utf8').trimEnd().split('\n').length,1);
    });
    await check('sample report renders at desktop size', async () => {
      await page.setViewportSize({width:1280,height:950});
      await page.goto(pathToFileURL(path.join(sample,'report.html')).href);
      await page.locator('#match-rows button').first().click();
      await page.screenshot({path:path.join(evidence,'report.png'),fullPage:true});
    });
    await check('no JavaScript errors, CSP violations, or external requests', async () => {
      assert.deepEqual(errors,[]); assert.deepEqual(violations,[]); assert.deepEqual(external,[]);
    });
    const result={complete:true,chromium:browser.version(),playwright:require(path.join(process.env.MOTIF_PLAYWRIGHT_MODULE || 'playwright','package.json')).version,
      node:process.version,offline:true,external_requests:external,errors,csp_violations:violations,checks,
      max_mounted_match_rows:75,fixture_occurrences:data.matches.length,viewports:[[1280,950],[375,812]]};
    fs.writeFileSync(path.join(evidence,'browser.json'),JSON.stringify(result,null,2)+'\n');
    console.log(JSON.stringify(result));
  } finally {await browser.close();}
}
main().catch(error=>{console.error(error);process.exitCode=1;});

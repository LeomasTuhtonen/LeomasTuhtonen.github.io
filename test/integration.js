const fs = require('fs');
const assert = require('assert');

const html = fs.readFileSync('index.html', 'utf8');

assert(/id=["']addSpanBtn["']/.test(html), 'Add Span button not found');
assert(/id=["']exportBtn["']/.test(html), 'Export button not found');
assert(/id=["']frameSaveBtn["']/.test(html), 'Frame Save button not found');
assert(/id=["']frameLoadBtn["']/.test(html), 'Frame Load button not found');

console.log('Integration tests passed');

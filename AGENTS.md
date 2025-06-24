# Repository Guidance for Codex

The following checks must be run from the repository root before committing any changes:

- **Unit tests:** `npm test`
- **Linting:** `npm run lint`
- **Integration tests:** `npm run test:integration`

Codex should ensure these commands complete successfully. If they fail due to missing dependencies or other issues, report the failure in the PR description.

---

## Static analysis for temporal-dead-zone errors

- Ensure ESLint is configured with the following rules set to “error”:
  - `no-use-before-define`
  - `no-undef`
- Run a strict lint step before committing:
  ```bash
  npm run lint:strict
  ```
- Example `package.json` entry:
  ```json
  {
    "scripts": {
      "lint:strict": "eslint \"src/**/*.{js,ts}\" --rule \"no-use-before-define: error\" --rule \"no-undef: error\""
    }
  }
  ```

## E2E smoke tests with console-error checking

- Add a headless-browser smoke-test script (using Puppeteer or Playwright) that:
  1. Launches your app (e.g. `npm start`)
  2. Navigates to each critical page
  3. Listens for `console.error` and uncaught exceptions
  4. Exits with a non-zero code if any error is detected
- Run this as part of your validation pipeline:
  ```bash
  npm run test:smoke
  ```
- Example `package.json` entry:
  ```json
  {
    "scripts": {
      "test:smoke": "node ./scripts/run-smoke-tests.js"
    }
  }
  ```
- (In `./scripts/run-smoke-tests.js`, hook `page.on('console', msg => { if (msg.type() === 'error') process.exit(1) })`.)

---

By enforcing these additional steps, Codex will catch “Cannot access ‘sectionChart’ before initialization” and other runtime reference errors at build or PR-validation time, rather than in production.

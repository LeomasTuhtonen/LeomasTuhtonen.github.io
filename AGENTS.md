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


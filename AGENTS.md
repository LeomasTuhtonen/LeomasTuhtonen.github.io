# Repository Guidance for Codex

The following checks must be run from the repository root before committing any changes:

1. **Unit tests**: `npm test`
2. **Linting**: `npm run lint`
3. **Integration tests**: `npm run test:integration`

Codex should ensure these commands complete successfully. If they fail due to missing dependencies or other issues, report the failure in the PR description.

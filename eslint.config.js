export default [
  {
    files: ['**/*.js'],
    languageOptions: {
      ecmaVersion: 2020,
      sourceType: 'module',
      globals: { window: 'readonly', document: 'readonly' }
    },
    rules: {}
  }
];

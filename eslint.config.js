export default [
  {
    ignores: ['libs/**']
  },
  {
    files: ['**/*.js'],
    languageOptions: {
      ecmaVersion: 2020,
      sourceType: 'module',
      globals: {
        window: 'readonly',
        document: 'readonly',
        require: 'readonly',
        module: 'readonly',
        process: 'readonly',
        console: 'readonly',
        setTimeout: 'readonly'
      }
    },
    rules: {
      'no-use-before-define': 'error',
      'no-undef': 'error'
    }
  }
];

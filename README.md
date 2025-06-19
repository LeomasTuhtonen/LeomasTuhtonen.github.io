# structural-calc

A simple browser-based beam calculator using 2D finite element analysis. No installation is required; just open `index.html` in a modern web browser.

Features:
- Define any number of spans and lengths.
- Choose from a list of common European steel cross-sections (IPE/HEA) or enter your own moment of inertia.
- Add point and line loads at arbitrary locations.
- Results update automatically showing support reactions and diagrams for shear force, bending moment and deflection using modern charts.

The solver uses Euler-Bernoulli beam elements with adjustable element divisions per span.

## Running tests

```
npm test
```


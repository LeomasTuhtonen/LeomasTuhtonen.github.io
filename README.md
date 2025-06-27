# Beam Calculator

A browser-based beam calculator using 2D finite element analysis. Simply open
`index.html` in a modern web browser to get started.

Features:
- Define any number of spans with adjustable node spacing.
- Built-in steel and timber cross-section libraries with grade selection.
- Automatic self weight plus customizable point and line loads.
- User-defined load combinations.
- Live charts for loads, shear, bending moment and deflection.
- Design tab showing section properties and checks for bending, shear and lateral-torsional buckling.
- Export analysis results to PDF.

The solver uses Euler-Bernoulli beam elements with adjustable element divisions per span.

See [beam_buckling_fixed.md](beam_buckling_fixed.md) for the lateral-torsional buckling equations used when calculating design resistance.

## Running tests

```
npm test
```


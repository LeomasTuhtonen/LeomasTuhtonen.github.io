# ThunderStruct

A browser-based beam calculator using 2D finite element analysis. Simply open
`index.html` in a modern web browser to get started.

**Disclaimer:** This tool is provided for educational and preliminary design purposes only. The results should not be used for final structural design without independent verification by a qualified professional engineer. The author assumes no responsibility or liability for any errors or omissions in the content of this site or for any actions taken in reliance on the information provided.

Features:
- Define any number of spans with adjustable node spacing.
- Built-in steel and timber cross-section libraries with grade selection.
- Automatic self weight plus customizable point and line loads.
- User-defined load combinations.
- Live charts for loads, shear, bending moment and deflection.
- Frame tab with beams, supports and node loads plus new member point and line loads with on-diagram illustrations.
- Design tab showing section properties and checks for bending, shear and lateral-torsional buckling.
- Export analysis results to PDF.

The solver uses Euler-Bernoulli beam elements with adjustable element divisions per span.

See [beam_buckling_fixed.md](beam_buckling_fixed.md) for the lateral-torsional buckling equations used when calculating design resistance.

## Running tests

```
npm test
```


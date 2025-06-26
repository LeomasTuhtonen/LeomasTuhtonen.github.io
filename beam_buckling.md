# Lateral-Torsional Buckling Procedures

## Steel Beams (Eurocode 3: EN 1993-1-1)

1. **Determine the true unbraced length,** \(L_b\), measured between restraints of the shear centre or warping (not just flange braces).
2. **Compute the elastic critical moment,** \(M_{cr}\):
   \[
   M_{cr}
   = \frac{\pi}{L_b} \sqrt{E I_w \; G I_t} 
     \sqrt{1 + \dfrac{\pi^2 E I_w}{G I_t L_b^2}}
   \]
   - \(E\): Young’s modulus  
   - \(G\): Shear modulus  
   - \(I_w\): Warping constant  
   - \(I_t\): Torsional constant  
   - \(L_b\): Unbraced length (shear-centre restraints)

3. **Calculate non-dimensional slenderness and reduction factor:**
   \[
   \overline{\lambda}_{LT}
   = \sqrt{\frac{M_y / \gamma_{M1}}{M_{cr}}}
   \quad\longrightarrow\quad
   \chi_{LT}
   = \frac{1}{\varphi + \sqrt{\varphi^2 + \overline{\lambda}_{LT}^2}}
   \]
4. **Design buckling resistance:**
   \[
   M_{b,Rd} = \chi_{LT} \; \frac{M_y}{\gamma_{M1}}
   \]
5. **Notes:**
   - If braces only restrain lateral movement (no warping restraint), use the General Method (Clause 6.3.4).
   - For full restraint (lateral + warping), use the formula above (Clause 6.3.2).

---

## Timber Beams (Eurocode 5: EN 1995-1-1)

1. **Effective buckling length,** \(L_{cr}\), between lateral restraints of the compression flange.
2. **Elastic critical bending stress,** \(\sigma_{m,\mathrm{crit}}\):
   - **General section (Clause 6.3.3 eq. 6.31):**
     \[
     \sigma_{m,\mathrm{crit}}
     = \frac{\pi}{L_{cr} W_y} \sqrt{E I_z \; G I_t}
     \]
   - **Solid rectangular section (eq. 6.32):**
     \[
     \sigma_{m,\mathrm{crit}}
     = \frac{0.78 \, b^2 \, E}{h \, L_{cr}}
     \]
   - \(E\): Modulus of elasticity  
   - \(G\): Shear modulus  
   - \(I_z\): Second moment of area (weak axis)  
   - \(I_t\): Torsional constant  
   - \(W_y\): Section modulus (strong axis)  
   - \(b,h\): Width and height (rectangle)  
   - \(L_{cr}\): Unbraced length

3. **Non-dimensional slenderness and reduction factor (eq. 6.30–6.34):**
   \[
   \lambda_{rel,m}
   = \sqrt{\frac{f_{m,y,d}}{\sigma_{m,\mathrm{crit}}}}
   \quad\longrightarrow\quad
   k_{crit}
   = \frac{1}{\lambda_{rel,m} + \sqrt{\lambda_{rel,m}^2 + 0.25}}
   \]
4. **Design moment resistance:**
   \[
   M_{b,Rd} = k_{crit} \; W_y \; f_{m,y,d}
   \]
5. **Notes:**
   - \(L_{cr}\) must be between actual lateral supports of the compression flange.
   - Use the simplified eq. 6.32 only for solid rectangular cross-sections.


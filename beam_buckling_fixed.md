# Lateral‑Torsional Buckling Procedures – **Corrected & Conservative**

This note replaces the previous *beam_buckling.md*.  The steel‑beam section now uses the **complete linear‑buckling (LBA) expression** with axis nomenclature that matches EN 1993‑1‑1 (weak‑axis inertia **I<sub>z</sub>**, strong‑axis inertia **I<sub>y</sub>**).  All moment‑gradient factors *C*<sub>i</sub> and the warping‑fixity factor *k*<sub>w</sub> are preset to **the most unfavourable values** so the resistance you obtain is always on the safe side.

---

## 1 Steel beams (EN 1993‑1‑1 Annex A & NCCI SN003)

### 1.1 Elastic critical moment \(M_{cr}\)
For doubly‑symmetric I‑sections subjected to equal end moments:

\[
M_{cr}
=
C_1\,\frac{\pi}{L_b}\,\sqrt{E\,I_z\,G\,I_t}
\;\sqrt{\,k_w
  + C_2\,\frac{\pi^2 E I_w}{C_1^{2} G I_t L_b^{2}}
  + C_3\,\frac{\pi^4 E^{2} I_z I_w}{C_1^{2} G^{2} I_t^{2} L_b^{4}}\,}\tag{1}
\]

| symbol | description | units |
|--------|-------------|-------|
| \(L_b\) | unbraced length between points that restrain **both** lateral translation *and* warping | m |
| \(E\) | Young’s modulus (≈ 210 GPa for steel) | Pa |
| \(G\) | shear modulus (≈ 81 GPa) | Pa |
| \(I_z\) | weak‑axis second moment of area (smaller) | m<sup>4</sup> |
| \(I_t\) | St Venant torsion constant | m<sup>4</sup> |
| \(I_w\) | warping constant | m<sup>6</sup> |
| \(C_1,C_2,C_3\) | moment‑gradient factors (see §1.2) | – |
| \(k_w\) | warping‑fixity factor (see §1.2) | – |

> **Axis convention** – *y*‑axis is the **strong** bending axis; *z*‑axis is the **weak** axis.  Thus **I<sub>y</sub> > I<sub>z</sub>** for rolled I‑shapes.

### 1.2 Conservative factor set
The factors below give the **minimum possible \(M_{cr}\)** for typical geometries and loading, hence the lowest (conservative) buckling resistance.  Replace them with project‑specific values if you wish to reduce conservatism.

| factor | conservative value | reason |
|--------|--------------------|--------|
| \(C_1\) | **1.00** | lowest value (uniform moment diagram) |
| \(C_2\) | **0.00** | omits favourable moment gradient term |
| \(C_3\) | **0.00** | omits higher‑order favourable term |
| \(k_w\) | **1.00** | assumes warping is *free* at both ends |

### 1.3 Non‑dimensional slenderness & reduction factor

\[
\overline{\lambda}_{LT} = \sqrt{\frac{M_y/\gamma_{M1}}{M_{cr}}}
\quad\Longrightarrow\quad
\varphi = 0.5\,[1 + \alpha_{LT}(\overline{\lambda}_{LT} - 0.2) + \overline{\lambda}_{LT}^2]
\]
\[
\chi_{LT} = \frac{1}{\varphi + \sqrt{\varphi^{2} + \overline{\lambda}_{LT}^{2}}}
\tag{2}
\]
with the recommended imperfection factor **\(\alpha_{LT} = 0.34\)**.

### 1.4 Design moment resistance

\[
M_{b,Rd} = \chi_{LT} \; \frac{M_y}{\gamma_{M1}}
\tag{3}
\]

---

## 2 Timber beams (EN 1995‑1‑1, Clause 6.3.3) – unchanged

### 2.1 Elastic critical bending stress
For solid or glulam rectangular sections the Eurocode offers a closed form:

\[
\sigma_{m,\,crit} = 0.78\;\frac{b^{2}E}{h\,L_{cr}}\tag{4}
\]

### 2.2 Reduction factor

\[
\lambda_{rel,m} = \sqrt{\frac{f_{m,y,d}}{\sigma_{m,\,crit}}}
\quad\Longrightarrow\quad
k_{crit} = \frac{1}{\lambda_{rel,m} + \sqrt{\lambda_{rel,m}^{2} + 0.25}}
\tag{5}
\]

### 2.3 Design moment resistance

\[
M_{b,Rd} = k_{crit}\;W_y\;f_{m,y,d}\tag{6}
\]

> *Note* For non‑rectangular or curved timber members, use the more general procedure of Clause 6.3.3(2) with the elastic critical stress from Annex B.

---

### 3 References
* EN 1993‑1‑1:2005 §6.3.2 & Annex A (2018 Consolidated Edition)
* NCCI SN003b: "Buckling curves for LTB" (SCI, 2019)
* EN 1995‑1‑1:2014 §6.3.3


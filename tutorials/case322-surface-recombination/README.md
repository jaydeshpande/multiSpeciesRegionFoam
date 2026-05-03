# case322 — Surface Recombination / Dissociation Kinetics

**Solver:** `foamMultiRun` with `speciesSolid` region module

---

## What this case covers

When a metal is exposed to a molecular gas (H₂, D₂, T₂), two kinetic processes compete at the surface:

1. **Dissociation (absorption):** H₂ → 2H,  flux = Kd(T) · p_gas  [mol/(m²·s)]
2. **Recombination (desorption):** 2H → H₂,  flux = Kr(T) · C²  [mol/(m²·s)]

The net flux into the solid is:

```
D · ∂C/∂n  =  Kd(T) · pGas  −  Kr(T) · C²
```

At equilibrium (no net flux): C_eq = √(Kd · pGas / Kr) — which is Sieverts' law.

The **Damköhler number** Da = Kr · C_eq · L / D quantifies whether kinetics or diffusion limits the permeation:

| Da | Limiting step | Surface concentration |
|----|--------------|----------------------|
| Da ≫ 1 | Diffusion limited | C_s ≈ C_eq (Sieverts equilibrium) |
| Da ≈ 1 | Mixed | C_s = (√5−1)/2 · C_eq ≈ 0.618 C_eq |
| Da ≪ 1 | Kinetically limited | C_s ≪ C_eq |

This case demonstrates:

1. That the `surfaceRecombination` BC recovers the Sieverts `fixedValue` limit as kinetics become fast (Da → ∞).
2. The magnitude of flux reduction when surface kinetics are rate-limiting (Da = 1 here → 38.2% reduction).
3. The Robin-type BC formulation and its Picard convergence via PIMPLE outer iterations.

---

## Geometry

Single metal slab, 0.5 mm thick:

```
gas (pGas = 1 Pa)            purge (C = 0)
     │← metal: D=5e-11, L=0.5mm, 50 cells →│
  x=0  (surfaceRecombination BC)          x=0.5 mm
```

---

## Sub-cases

### sieverts_bc (reference — Da → ∞)

Left face: `fixedValue C = C_eq = 1.0 mol/m³`.

This models **instantaneous** surface kinetics — the Sieverts equilibrium is maintained at every instant. It provides the upper bound on permeation flux for the given equilibrium concentration.

```
J_ref = D · C_eq / L = 5×10⁻¹¹ × 1.0 / 5×10⁻⁴ = 1.0×10⁻⁷ mol/(m²·s)
```

### surface_recombination (Da = 1)

Left face: `type surfaceRecombination` with Kd = Kr = 1×10⁻⁷.

Parameters are chosen so that:
- Equilibrium: C_eq = √(Kd · pGas / Kr) = 1.0 mol/m³  (same as reference)
- Damköhler: Da = Kr · C_eq · L / D = 1.0  (kinetics and diffusion equally limiting)

Steady-state surface concentration (solving Kr·Cₛ² + (D/L)·Cₛ − Kd·pGas = 0):
```
1×10⁻⁷·Cₛ² + 1×10⁻⁷·Cₛ − 1×10⁻⁷ = 0
Cₛ² + Cₛ − 1 = 0
Cₛ = (√5 − 1)/2 ≈ 0.6180 mol/m³
```

The golden-ratio conjugate appears naturally when Da = 1 with Kd = Kr and pGas = 1.

```
J_surf = D · Cₛ / L ≈ 6.18×10⁻⁸ mol/(m²·s)
J_surf / J_ref ≈ 0.618  (38.2% reduction from surface kinetic resistance)
```

---

## Physical parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| D | 5×10⁻¹¹ m²/s | Diffusivity (constant, Ea=0) |
| L | 0.5 mm | Slab thickness |
| pGas | 1.0 Pa | Gas-phase H₂ partial pressure |
| Kd | 1×10⁻⁷ mol/(m²·s·Pa) | Dissociation rate coefficient |
| Kr | 1×10⁻⁷ m⁴/(mol·s) | Recombination rate coefficient |
| C_eq | 1.0 mol/m³ | Sieverts equilibrium concentration |
| T | 300 K | Isothermal |
| Da | 1.0 | Damköhler number |

---

## Analytical steady state

At steady state, the slab has a linear concentration profile (constant D) with C_s at the gas face and 0 at the purge face.

The surface concentration Cₛ satisfies the Robin balance:

```
Kd·pGas − Kr·Cₛ² = D·Cₛ/L        (absorption = permeation at SS)
```

Rearranging: Kr·Cₛ² + (D/L)·Cₛ − Kd·pGas = 0

With the chosen parameters (Kd=Kr=1e-7, D/L=1e-7, Kd·pGas=1e-7):
```
Cₛ² + Cₛ − 1 = 0  →  Cₛ = (√5−1)/2 ≈ 0.6180 mol/m³
```

General case (varying Da):
```
Cₛ = [−(D/L) + √((D/L)² + 4·Kr·Kd·pGas)] / (2·Kr)
```

As Kr → ∞ (keeping Kd/Kr = C_eq²/pGas fixed): Cₛ → C_eq (Sieverts limit).

---

## Flux reduction vs Damköhler number

The surface kinetic resistance adds in series with the diffusion resistance. For general Da:

```
1/J_total = 1/J_kinetic + 1/J_diffusion
```

| Da | Cₛ/C_eq | J/J_ref |
|----|---------|---------|
| 0.01 | ≈ 0.01 | ≈ 0.01 |
| 0.1  | ≈ 0.09 | ≈ 0.09 |
| 1.0  | 0.618  | 0.618  |
| 10   | ≈ 0.95 | ≈ 0.95 |
| 100  | ≈ 0.995| ≈ 0.995|

The `surfaceRecombination` BC spans the full range from kinetically limited to diffusion limited within a single Robin condition.

---

## How to run

```bash
cd tutorials/case322-surface-recombination
./Allrun
```

### Validate

```bash
python3 validation/recombination_compare.py
python3 validation/recombination_compare.py --plot
```

PASS criterion: permeation flux at the last write time within 5% of the analytical value for both sub-cases.

### Clean

```bash
./Allclean
```

---

## Implementation notes

### Picard convergence

The C² term is nonlinear. The BC uses Picard (fixed-point) iteration: at each PIMPLE outer iteration, C² is evaluated at the **previous-iterate** cell-centre concentration, making the gradient explicit:

```
gradient() = [Kd(T)·pGas − Kr(T)·C_cell²] / D_patch
```

After the species equation is solved, C_cell is updated, and the next outer iteration re-evaluates the gradient. Convergence is reached when the residual criterion in `fvSolution` is satisfied.

### Temperature dependence

Kd and Kr follow the Arrhenius form:
```
Kd(T) = Kd0 · exp(−Ea_d / (R·T))
Kr(T) = Kr0 · exp(−Ea_r / (R·T))
```
Setting Ea = 0 gives temperature-independent kinetics (isothermal, as in this case).
For real materials (e.g., Pd, SS316) the activation energies set recombination to dominate at low T and dissociation at high T.

---

## Key references

1. Pick, M.A. & Sonnenberg, K., *J. Nucl. Mater.* **131**, 208–220 (1985) — surface recombination model.
2. Anderl, R.A. *et al.*, *J. Nucl. Mater.* **196–198**, 986–991 (1992) — Kd, Kr for stainless steels.
3. Hattab, N. *et al.*, *Fusion Engineering and Design* **202**, 114362 (2024).

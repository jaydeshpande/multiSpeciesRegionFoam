# case319 — Hydrogen Permeation Barrier: WC Coating on SS316 Tube

**Solver:** `foamMultiRun` with `speciesSolid` region module

---

## What this case covers

A hydrogen-carrying tube is modelled as a solid wall consisting of two concentric
cylindrical layers.  The planar (1-D Cartesian) approximation is used, valid when
the wall thickness is small relative to the tube radius ($t/r \ll 1$).

Two configurations are compared:

| Sub-case | Geometry | Solver regions |
|----------|----------|----------------|
| `no_coating/` | SS316 wall only, 0.5 mm | 1 region: `steel` |
| `with_coating/` | WC coating 0.1 mm + SS316 0.5 mm | 2 regions: `wc`, `steel` |

The key output is the **permeation flux** at the outer wall vs time.  The WC
coating dramatically reduces hydrogen permeation — this tutorial demonstrates and
quantifies that effect.

---

## Governing equations

Each solid region solves:

$$\frac{\partial C}{\partial t} = \nabla \cdot (D(T)\,\nabla C)$$

with temperature-dependent Arrhenius diffusivity:

$$D(T) = D_0 \exp\!\left(-\frac{E_a}{RT}\right)$$

At the WC/SS316 interface, thermodynamic equilibrium is enforced by the
`sievertsCoupledMixed` boundary condition (Sieverts law):

$$\frac{C_\text{SS316}}{C_\text{WC}} = \frac{K_{s,\text{SS316}}}{K_{s,\text{WC}}} = 10$$

The concentration **jumps** by a factor of 10 when crossing from WC (low
solubility) into SS316 (high solubility), at equal hydrogen chemical potential.

---

## Material properties

### Stainless Steel 316
*Source: Forcey et al. 1988, J. Nucl. Mater. 160, 153–159 (Table 1, 316L SS)*

| Parameter | Value | Description |
|-----------|-------|-------------|
| $D_0$ | $5.08 \times 10^{-7}$ m²/s | Pre-exponential diffusivity |
| $E_{a,D}$ | 71.3 kJ/mol (= 8580 K × R) | Diffusion activation energy |
| $D(1073\,\text{K})$ | $1.73 \times 10^{-10}$ m²/s | Value at operating temperature |
| $K_{s,0}$ | $3.89 \times 10^{-4}$ mol/(m³·Pa^0.5) | Sieverts pre-exponential |
| $K_s(1073\,\text{K})$ | $2.05 \times 10^{-3}$ mol/(m³·Pa^0.5) | Endothermic dissolution |
| $C_\text{inner}$ | **0.651 mol/m³** | $K_s \times \sqrt{P}$ at 1 bar, 1073 K |

### Tungsten Carbide (WC)
*Source: Diffusivity estimated to reproduce PRF ≈ 200 vs SS316 at 1073 K,
consistent with experimentally reported WC coating PRF ranges of 10²–10³
(Hollenberg et al. 1995, Fusion Eng. Des. 28, 190–208).*

| Parameter | Value | Description |
|-----------|-------|-------------|
| $D_0$ | $1.00 \times 10^{-7}$ m²/s | Pre-exponential diffusivity |
| $E_{a,D}$ | 97.8 kJ/mol | Higher barrier than SS316 |
| $D(1073\,\text{K})$ | $1.73 \times 10^{-12}$ m²/s | **100× lower than SS316** |
| $K_{s,\text{WC}}$ | $K_{s,\text{SS316}} / 10$ | Lower H solubility in ceramic |
| $C_\text{inner}$ | **0.0651 mol/m³** | $K_s \times \sqrt{P}$ at 1 bar, 1073 K |

---

## Analytical steady-state solution

At steady state the flux $J$ is constant through both layers.  Letting $r = K_{s,\text{SS}} / K_{s,\text{WC}} = 10$:

$$C_\text{WC}^\text{if} = \frac{(D_\text{WC}/L_\text{WC})\,C_\text{in,WC}}{D_\text{WC}/L_\text{WC} + r\,D_\text{SS}/L_\text{SS}}$$

$$C_\text{SS}^\text{if} = r\,C_\text{WC}^\text{if}$$

$$J_\text{bare} = \frac{D_\text{SS}\,C_\text{in,SS}}{L_\text{SS}} \approx 2.25 \times 10^{-7}\ \text{mol/(m}^2\text{·s)}$$

$$J_\text{coat} = \frac{D_\text{SS}\,C_\text{SS}^\text{if}}{L_\text{SS}} \approx 1.12 \times 10^{-9}\ \text{mol/(m}^2\text{·s)}$$

$$\text{PRF} = \frac{J_\text{bare}}{J_\text{coat}} \approx \mathbf{201}$$

---

## Mesh and time parameters

| Parameter | no_coating | with_coating |
|-----------|-----------|--------------|
| WC cells | — | 20 (Δx = 5 μm) |
| SS316 cells | 100 (Δx = 5 μm) | 100 (Δx = 5 μm) |
| $\Delta t$ | 1 s | 1 s |
| End time | 4000 s | 4000 s |
| $\tau_{SS}$ (SS316) | 1446 s (= 2.77 end times) | 1446 s |
| $\tau_{WC}$ (WC) | — | 5780 s |
| Write interval | 400 s | 400 s |

At $t = 4000$ s, the uncoated SS316 is at $\approx 97\%$ of its steady-state
flux.  The coated case is still in transient but has converged enough to
demonstrate the large PRF.

---

## Interface BC: `sievertsCoupledMixed`

Both sides of the WC/SS316 interface must specify each material's Sieverts
constant so the solver can enforce chemical potential continuity:

```
# 0/wc/C_H2  —  wc_to_steel patch
type        sievertsCoupledMixed;
Ks    { X0 100;  Ea 0; }     # Ks_WC (this side)
KsNbr { X0 1000; Ea 0; }    # Ks_SS316 (neighbour)

# 0/steel/C_H2  —  steel_to_wc patch
type        sievertsCoupledMixed;
Ks    { X0 1000; Ea 0; }    # Ks_SS316 (this side)
KsNbr { X0 100;  Ea 0; }    # Ks_WC (neighbour)
```

The ratio $K_{s,\text{SS316}} / K_{s,\text{WC}} = 1000/100 = 10$ sets the
concentration jump at the interface.

---

## How to run

```bash
cd tutorials/case319-permeation-barrier
./Allrun
```

Or run each sub-case individually:

```bash
cd no_coating   && ./Allrun && cd ..
cd with_coating && ./Allrun && cd ..
```

### Validate and plot

```bash
# Text report with PASS/FAIL on no_coating steady-state flux
python3 validation/permeation_compare.py

# Save permeation_comparison.png (flux vs time + concentration profiles)
python3 validation/permeation_compare.py --plot
```

### Clean

```bash
./Allclean
```

---

## Expected results

| Metric | Value |
|--------|-------|
| $J_\text{bare}$ at $t=4000$ s | $\approx 2.2 \times 10^{-7}$ mol/(m²·s) |
| $J_\text{coat}$ at $t=4000$ s | $\approx 10^{-9}$ mol/(m²·s) |
| Simulated PRF | $\approx 200$ |
| Analytical PRF | 201 |

The `no_coating` case passes if the outer-wall flux at the last write time is
within 5% of the analytical steady-state value.

---

## Physical context

Permeation barriers are important in:
- **Hydrogen storage / fuel cells** — limiting hydrogen loss through metal tanks
- **High-temperature processing** — protecting structural metals from hydrogen
  embrittlement
- **Fusion energy** — preventing tritium permeation through first-wall and
  blanket structures into the coolant
- **Fission energy** — preventing tritium permeation through structural materials in Molten Salt Reactors

WC is an attractive candidate due to its high hardness, chemical stability at
elevated temperature, and low hydrogen permeability compared to metals.  Real
coating PRF depends strongly on coating density, grain boundaries, and defects —
dense, well-adhered coatings achieve PRF > 10³.

---

## Key references

1. Forcey K.S., Ross D.K., Simpson J.C.B., Evans D.S., "Hydrogen transport and
   solubility in 316L and 1.4914 steels for fusion reactor applications,"
   *J. Nucl. Mater.* **160**, 153–159 (1988).
2. Hollenberg G.W. *et al.*, "The effect of irradiation on the stability and
   properties of tritium permeation barriers,"
   *Fusion Eng. Des.* **28**, 190–208 (1995).
3. Pierson H.O., *Handbook of Refractory Carbides and Nitrides*, Noyes
   Publications, 1996.

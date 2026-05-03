# case320 — Shell-and-Tube Heat Exchanger: Coupled Heat and H₂ Permeation

**Solver:** `foamMultiRun` with `speciesSolid` region module

---

## What this case covers

A shell-and-tube intermediate heat exchanger (IHX) transfers heat from a hot primary molten salt (FLiBe) to a cooler secondary salt (Hitec).  The primary salt carries dissolved hydrogen (a tritium surrogate in fusion-relevant contexts), which permeates through the metal tube wall into the secondary loop.

The case models a 1-D cross-section through the tube wall — the dominant direction for both heat conduction and species diffusion — using **three coupled solid regions**:

| Region | Material | Thickness | Role |
|--------|----------|-----------|------|
| `flibe` | FLiBe proxy | 0.1 mm | Primary molten salt boundary layer |
| `wall`  | SS316 | 0.5 mm | Tube wall — rate-limiting permeation barrier |
| `hitec` | Hitec proxy | 0.1 mm | Secondary salt boundary layer |

FLiBe (66.7% LiF + 33.3% BeF₂) is a candidate primary coolant for molten salt reactors and fusion tritium breeding blankets.  Hitec (50 wt% NaNO₃ + 50 wt% KNO₃) is a common secondary heat transfer salt.

The physical counterflow operation (FLiBe in +x, Hitec in −x along the tube axis) determines the bulk temperatures imposed as boundary conditions at the two outer faces.  The axial variation of temperature along the 1 m tube length would require a fluid solver module; this case captures the cross-sectional transport at a representative axial location.

---

## Key physics demonstrated

### 1 — Conjugate heat transfer (three-region)

The energy equation is solved simultaneously in all three solid regions.  Heat flux is continuous at both interfaces via `coupledTemperature` boundary conditions.  At steady state the temperature profile is piecewise linear in each layer with slope $-J_q/k_i$:

$$J_q = \frac{T_{\text{FLiBe}} - T_{\text{Hitec}}}{\sum_i L_i / k_i} \approx 1.34 \times 10^6 \ \text{W/m}^2$$

| Location | T (steady state) |
|----------|-----------------|
| FLiBe outer face | 1073 K (800°C) |
| FLiBe/wall interface | ≈ 951 K |
| wall/Hitec interface | ≈ 907 K |
| Hitec outer face | 673 K (400°C) |

The wall carries only a modest temperature drop (44 K over 0.5 mm) because its thermal conductivity ($k = 15$ W/(m·K)) is much higher than the fluid layers.  Most of the 400 K temperature difference falls across the low-conductivity FLiBe and Hitec proxy layers.

### 2 — Arrhenius temperature-dependent permeation

Hydrogen diffusivity in SS316 follows the Arrhenius law:

$$D(T) = D_0 \exp\!\left(-\frac{E_a}{RT}\right) = 5.08 \times 10^{-7} \exp\!\left(-\frac{8580}{T}\right) \ \text{m}^2/\text{s}$$

The 44 K gradient across the wall produces a **56% variation** in $D(T)$:

| Wall position | T (K) | D (m²/s) |
|---------------|-------|----------|
| FLiBe side | 951 | $6.16 \times 10^{-11}$ |
| Hitec side | 907 | $3.94 \times 10^{-11}$ |

Because $D$ is higher on the hot (FLiBe) side, concentration decreases more steeply near the FLiBe/wall interface and more gradually near the Hitec side — the steady-state concentration profile **bows toward the Hitec side** relative to the straight line that would result from a constant $D$.

### 3 — Multi-region Sieverts coupling

The `sievertsCoupledMixed` boundary condition at both interfaces enforces flux continuity and thermodynamic equilibrium.  Because $K_{s,\text{FLiBe}} = K_{s,\text{wall}} = K_{s,\text{Hitec}} = 1000$ (equal Sieverts constants), the condition reduces to simple concentration continuity ($C$ and $J_C$ continuous), appropriate for the proxy-solid approximation.

---

## Physical parameters

### FLiBe primary salt
*Source: Romatoski & Hu, Nucl. Technol. 205 (2019) 1367*

| Property | Value |
|----------|-------|
| $\rho$ | 1940 kg/m³ |
| $C_p$ | 2386 J/(kg·K) |
| $k$ | 1.1 W/(m·K) |
| $T_{\text{bulk}}$ | 1073 K (800°C) |
| $C_{H_2}$ | **1.0 mol/m³** (dissolved tritium/H₂) |

### SS316 tube wall
*Source: Forcey et al. 1988, J. Nucl. Mater. 160, 153–159*

| Property | Value |
|----------|-------|
| $\rho$ | 7960 kg/m³ |
| $C_v$ | 530 J/(kg·K) |
| $k$ | 15 W/(m·K) |
| $D_0$ | $5.08 \times 10^{-7}$ m²/s |
| $E_{a,D}$ | 71.3 kJ/mol (= 8580 K × R) |
| Thickness | 0.5 mm |

### Hitec secondary salt
*Source: Bradshaw & Siegel, ASME ESDA2008-59283*

| Property | Value |
|----------|-------|
| $\rho$ | 1870 kg/m³ |
| $C_p$ | 1560 J/(kg·K) |
| $k$ | 0.57 W/(m·K) |
| $T_{\text{bulk}}$ | 673 K (400°C) |
| $C_{H_2}$ | **0** (continuously removed by H₂ purge) |

---

## Analytical steady-state solution

At steady state the species flux $J_C$ is constant through all layers.  The wall dominates the transport resistance because $D_{\text{wall}} \ll D_{\text{proxy}}$:

$$J_C = \frac{C_{\text{in}}}{R_{C,\text{total}}} \approx \frac{C_{\text{in}}}{\displaystyle\int_0^{L_{\text{wall}}} \frac{\mathrm{d}x}{D(T(x))}}$$

where $T(x)$ is the local temperature (linear across the wall).  This integral has no closed form when $T(x)$ is linear, but is easily evaluated numerically.  The validation script computes it with 5 000-point trapezoidal quadrature.

The steady-state concentration in the wall satisfies:

$$C(x) = C_{\text{fw}} - J_C \int_0^x \frac{\mathrm{d}s}{D(T(s))}$$

where $C_{\text{fw}} \approx C_{\text{in}}$ (proxy resistance negligible).

---

## Mesh and time parameters

| Parameter | Value |
|-----------|-------|
| FLiBe proxy cells | 5 ($\Delta x = 20\ \mu\text{m}$) |
| Wall cells | 50 ($\Delta x = 10\ \mu\text{m}$) |
| Hitec proxy cells | 5 ($\Delta x = 20\ \mu\text{m}$) |
| $\Delta t$ | 100 s |
| End time | 25 000 s |
| $\tau_{\text{wall}}$ (species) | $\approx$ 4 500 s |
| Write interval | 2 500 s |

The thermal timescale of each layer is < 0.1 s; temperature reaches steady state within the first time step.  The species timescale governs the simulation.

---

## How to run

```bash
cd tutorials/case320-shell-tube-hx
./Allrun
```

### Validate

```bash
# Text report with PASS/FAIL on steady-state permeation flux
python3 validation/hx_permeation.py

# Save hx_permeation.png (flux history + temperature + concentration profiles)
python3 validation/hx_permeation.py --plot
```

### Clean

```bash
./Allclean
```

---

## Expected results

| Metric | Value |
|--------|-------|
| Heat flux $J_q$ | $\approx 1.34 \times 10^6$ W/m² |
| T at FLiBe/wall interface | ≈ 951 K |
| T at wall/Hitec interface | ≈ 907 K |
| Steady-state permeation flux $J_C$ | $\approx 9 \times 10^{-8}$ mol/(m²·s) |
| $D$ variation across wall | ≈ 1.56× (hot side to cold side) |

The concentration profile in the wall is visibly non-linear: the bowing toward the Hitec side (the cooler, lower-$D$ region) is a direct signature of the Arrhenius temperature dependence.  The isothermal reference profile (at the hot-side $D$) would be a straight line.

PASS criterion: permeation flux at the last write time within 5% of the analytical steady-state value.

---

## Physical context

**Why this matters:**

In a molten-salt reactor or tritium breeding blanket, ⁶Li + n → ⁴He + T produces tritium (³H) that dissolves in the FLiBe primary coolant.  Tritium behaves chemically like hydrogen and permeates through metal walls — a serious radiological concern.  Understanding the permeation rate from the primary to the secondary loop is essential for:

- Licensing (tritium inventory in the secondary)
- Tritium recovery system design (matching the purge rate to the permeation flux)
- Barrier material selection (e.g., oxide coatings on the inner tube surface)

The temperature dependence of $D$ means that a hotter primary loop permeates tritium significantly faster.  A 100 K increase in the FLiBe temperature roughly doubles $D$ and hence the steady-state flux.

---

## Key references

1. Forcey K.S. *et al.*, "Hydrogen transport and solubility in 316L and 1.4914 steels for fusion reactor applications," *J. Nucl. Mater.* **160**, 153–159 (1988).
2. Romatoski R.R. & Hu L.-W., "Fluoride salt coolant properties for nuclear reactor applications: a review," *Nucl. Technol.* **205**, 1367–1388 (2019).
3. Bradshaw R.W. & Siegel N.P., "Molten nitrate salt development for thermal energy storage in parabolic trough solar power systems," *ASME ESDA2008-59283* (2008).
4. Zheng J. *et al.*, "Tritium permeation in molten salt reactors," *Nucl. Fusion* **61**, 026005 (2021).

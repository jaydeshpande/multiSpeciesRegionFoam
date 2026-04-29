# case313 — McNabb–Foster Trapping

**Solver:** `speciesSolid`  
**Coupling:** Single region, no interface  
**Trapping:** `McNabbFoster`

---

## What this case covers

This case validates the McNabb–Foster trapping model. A constant concentration is maintained at the inlet of a 1 m slab while the outlet is held at zero. Trapping retards the effective propagation of the diffusion front and delays the onset of steady-state permeation (the *breakthrough*). Two physically distinct regimes are provided:

| Regime | Trap energy ε | Behaviour |
|--------|---------------|-----------|
| **Weak trapping** (default) | 0.09 eV | Detrapping is fast; the system reaches a diffusion-modified steady state with effective diffusivity D_eff = D / S, S ≫ 1 |
| **Strong trapping** | 1.5 eV | Detrapping is negligible; traps fill irreversibly, delaying breakthrough by τ_bd ∝ ft·N·L² / (C_in·D) |

---

## Governing equations

**Mobile species (diffusion + trapping sink):**

```
∂C/∂t  =  D · ∂²C/∂x²  −  ∂Ct/∂t
```

**Trapped concentration (McNabb–Foster ODE, per cell, semi-implicit):**

```
∂Ct/∂t  =  (αt / N) · C · (nt − Ct)  −  αd(T) · Ct

αd(T)  =  αd0 · exp(−ε / (kB · T))
```

Discretised as:

```
Ct^{n+1}  =  [Ct^n + Δt·(αt/N)·C·nt]  /  [1 + Δt·((αt/N)·C + αd)]
```

---

## Geometry and parameters

| Parameter | Value |
|-----------|-------|
| Slab length L | 1 m |
| Cells | 100 (uniform) |
| Diffusivity D | 1.0 m² s⁻¹ (toy, to give tractable run time) |
| Inlet concentration C_in | 5.25×10⁻⁶ mol m⁻³ |
| Outlet concentration | 0 |
| Temperature (isothermal) | 1000 K |
| Molar trap-site density N | 0.0525 mol m⁻³ |
| Trap fraction ft | 0.1 → nt = 0.00525 mol m⁻³ |
| Trapping rate αt | 1×10¹⁵ s⁻¹ |
| Detrapping pre-exponential αd0 | 1×10¹³ s⁻¹ |

---

## Analytical steady state

At long times (or for weak trapping) the mobile concentration reaches the linear steady state:

```
C_ss(x)  =  C_in · (1 − x/L)
```

For **weak trapping** the approach to steady state is governed by the effective diffusivity and the trapping factor S:

```
S      =  1 + ft·αt / αd(T)
D_eff  =  D / S
τ_bc   =  L² · S / (2π² · D)       [breakthrough time]
```

With the default parameters at T = 1000 K:  
αd ≈ 3.52×10¹² s⁻¹, S ≈ 29.4, D_eff ≈ 0.034 m² s⁻¹, τ_bc ≈ 1.5 s. Run to t = 20 s.

For **strong trapping** (ε = 1.5 eV) the detrapping rate becomes negligible and traps fill irreversibly. The breakthrough time is instead

```
τ_bd  =  ft · N · L² / (2 · C_in · D)
```

With the strong-trapping parameters: τ_bd ≈ 500 s. Run to t = 3000 s (set `endTime 3000` in `controlDict`).

---

## How to run

**Weak trapping (default):**

```sh
cd tutorials/case313-trapping
./Allrun
```

**Strong trapping:**

Edit `constant/slab/speciesProperties` and change `epsilon 0.09` to `epsilon 1.5`, then set `endTime 3000` in `system/controlDict`, and re-run:

```sh
./Allclean && ./Allrun
```

**Validate:**

```sh
# Weak trapping (default)
python3 validation/trapping_compare.py

# Strong trapping
python3 validation/trapping_compare.py --epsilon 1.5
```

The script checks the steady-state mobile concentration profile against `C_ss(x) = C_in(1 − x/L)` at 2 % relative tolerance. It also prints the computed αd, S, D_eff, and the expected breakthrough time for the chosen regime.

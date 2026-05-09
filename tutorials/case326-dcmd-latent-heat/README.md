# case326 — DCMD with Latent Heat Flux BC

**Solver:** `speciesSolid` (all three regions)  
**Geometry:** Three-region 2-D slab — feed channel | PTFE membrane | permeate channel  
**New BCs:** `antoineEquilibrium` (membrane C_H2O), `latentHeatFlux` (feed/permeate T)  
**Reference:** Khalifa A. et al., *Desalination* **404** (2017) 22–34

---

## Purpose

case326 extends case325 by adding the dominant physical effect that was
missing — **latent heat coupling** at the membrane faces:

| Aspect | case325 | case326 |
|---|---|---|
| Membrane C_H2O BC | `antoineEquilibrium` | `antoineEquilibrium` (same) |
| Feed/permeate T BC | `coupledTemperature` | `latentHeatFlux` |
| Latent heat | absent | q = J · L_vap · M at membrane faces |

---

## New boundary conditions

### `antoineEquilibrium` (for C_H2O in membrane)

A `fixedValueFvPatchScalarField` that computes the equilibrium vapour
concentration from the local face temperature using the Antoine equation:

```
log₁₀(p_sat / mmHg) = A − B / (C + T [°C])
C_face = p_sat [Pa] / (R · T [K])
```

First introduced in case325 (replacing the static `fixedValue` used in
case324 with a dynamic BC that reads the live CFD temperature); case326
inherits it unchanged.

### `latentHeatFlux` (for T in fluid regions)

A `mixedFvPatchScalarField` that replaces `coupledTemperature` at the
membrane-fluid interfaces. It enforces **thermal continuity** (same as
`coupledTemperature`) and additionally applies the **latent heat flux**:

```
q_latent = J · L_vap · M_H2O   [W/m²]
```

where J [mol/(m²·s)] is read from the membrane C_H2O gradient computed
by the species solver — no algebraic duplication.

The interface temperature balance with latent heat:

```
κ_fluid · ∂T_fluid/∂n = κ_mem · ∂T_mem/∂n ± q_latent
```

- Feed side (`side evaporation`): q_latent subtracted — surface cools
- Permeate side (`side condensation`): q_latent added — surface warms

---

## Expected physical effect

Latent heat flux q ≈ J·L_vap·M ≈ 21 kW/m² (for J ≈ 20 L/(m²·h)).

With feed convective coefficient h_plug ≈ 1270 W/(m²·K):

```
ΔT_latent = q / h ≈ 21000 / 1270 ≈ 16 K
```

case326 should show:
- T_mf_326 ≈ T_mf_325 − 16 K ≈ 306 K
- T_mp_326 ≈ T_mp_325 + (q / h_perm) ≈ slightly higher
- TPC_326 < TPC_325 (stronger temperature polarisation)
- J_326 < J_325 (lower flux from lower T_mf − T_mp driving force)

---

## How to run

```sh
cd tutorials/case326-dcmd-latent-heat
./Allrun
```

Requires that the library is built with the new BCs:

```sh
cd ../../..   # → multiSpeciesRegionFoam root
./Allwmake
```

**Validate and compare with cases 324 and 325:**

```sh
python3 validation/dcmd_latent_compare.py
```

---

## Connection to the framework

case326 demonstrates:
1. **Named BCs as first-class library objects** — `antoineEquilibrium` and
   `latentHeatFlux` live in `libspeciesCoupling.so` alongside `sievertsCoupledMixed`
2. **Consistent flux reading** — `latentHeatFlux` reads J from the same
   membrane C_H2O field the species solver updates; no algebraic re-derivation
3. **Progressive model complexity**: case324 (analytical TPC) → case325
   (CFD plug flow, no latent heat) → case326 (CFD plug flow + latent heat)

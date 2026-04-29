# case311 — Slab Diffusion (erfc front)

**Solver:** `solid` (OpenFOAM built-in) — temperature field used as species proxy  
**Coupling:** Single region, no interface  
**Trapping:** None

---

## What this case covers

This is the first end-to-end validation case. It demonstrates that the OpenFOAM `solid` solver correctly reproduces 1-D transient diffusion from a step boundary condition, using the temperature field T [K] as a stand-in for species concentration C [mol m⁻³]. The mapping is exact because both fields satisfy the same parabolic equation with the same boundary conditions.

The `speciesSolid` solver module is not used here; this case exercises the underlying `solid` module so later cases can isolate regressions to the species extensions.

---

## Governing equation

```
∂C/∂t  =  D · ∂²C/∂x²      x ∈ [0, L],  t > 0
```

Boundary and initial conditions:

| Location | Condition |
|----------|-----------|
| x = 0 (inlet) | C = C_in = 1 mol m⁻³ (fixedValue) |
| x = L (outlet) | C = 0 mol m⁻³ (fixedValue) |
| t = 0 | C = 0 everywhere (empty slab) |

---

## Geometry and parameters

| Parameter | Value |
|-----------|-------|
| Slab length L | 1 mm |
| Cells | 200 (uniform) |
| Diffusivity D | 1×10⁻⁸ m² s⁻¹ (= κ / ρCv of the solid) |
| End time | 50 s |
| Write interval | 5 s |

The thermal–species analogy is set up by choosing the `physicalProperties` of the proxy solid so that κ / (ρ Cv) = D = 1×10⁻⁸ m² s⁻¹. The inlet temperature is offset by T_ref = 300 K to satisfy the `heSolidThermo` requirement T > 0; the validation script subtracts the reference before comparing.

---

## Analytical solution

While the diffusion front has not yet reached the far wall (`2√(Dt) ≪ L`), the semi-infinite approximation applies:

```
C(x, t)  =  erfc( x / (2√(Dt)) )
```

This is valid for t ≪ L²/D = 100 s. By t = 30 s the diffusion length is 2√(Dt) ≈ 0.35 mm, comfortably within the 1 mm slab.

At later times (t → ∞) the profile relaxes to the linear steady state `C_ss(x) = 1 − x/L`.

---

## How to run

```sh
cd tutorials/case311-slab-diffusion
./Allrun
```

The `Allrun` script calls `blockMesh` and then `foamMultiRun`.

**Validate against the analytical solution:**

```sh
python3 validation/erfc_compare.py --time 30
# or plot several snapshots:
python3 validation/erfc_compare.py --all
```

The script reads the `slab/T` field at the requested write time, converts T to C, and plots (or checks) against the erfc curve. `matplotlib` is optional; without it the script prints a pass/fail result to stdout.

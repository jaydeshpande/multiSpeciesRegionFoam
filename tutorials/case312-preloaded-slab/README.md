# case312 — Preloaded Slab (Fourier cosine series)

**Solver:** `speciesSolid`  
**Coupling:** Single region, no interface  
**Trapping:** None (`noTrapping`)

---

## What this case covers

This case validates transient diffusion from a non-uniform initial condition in an isolated (zero-flux) slab. The initial field is a step function — the left 10 % of the slab is preloaded with species, the remainder is empty — and the boundaries are sealed (zero-gradient). The concentration then redistributes by diffusion until uniform.

The case tests:

- `speciesSolid` as a standalone single-region solver (first use of the custom module)
- `setFields -region slab` to impose a non-trivial initial condition from the `setFieldsDict`
- The zero-gradient (isolated) boundary condition on both ends
- Accuracy of the transient solution against a known Fourier series

---

## Governing equation

```
∂C/∂t  =  D · ∂²C/∂x²      x ∈ [0, L],  t > 0
```

Boundary and initial conditions:

| Location | Condition |
|----------|-----------|
| x = 0 (left) | ∂C/∂x = 0 (zeroGradient) |
| x = L (right) | ∂C/∂x = 0 (zeroGradient) |
| t = 0 | C = C₀ for x ∈ [0, h], C = 0 for x ∈ (h, L] |

---

## Geometry and parameters

| Parameter | Value |
|-----------|-------|
| Slab length L | 0.1 m |
| Preloaded depth h | 0.01 m (10 % of L) |
| Initial concentration C₀ | 1 mol m⁻³ |
| Cells | 100 (uniform, dx = 1 mm) |
| Diffusivity D | 1×10⁻⁶ m² s⁻¹ |
| End time | 100 s |
| Write interval | 10 s |

---

## Analytical solution

Because the boundaries are sealed, mass is conserved. The exact solution is a Fourier cosine series:

```
C(x, t)  =  C₀·h/L  +  (2C₀/π) · Σ_{n=1}^{∞}  [sin(nπh/L) / n] · cos(nπx/L) · exp(−(nπ/L)²·D·t)
```

The DC term C₀h/L = 0.1 mol m⁻³ is the long-time uniform equilibrium.  
The series converges quickly: 500 terms are used in the validation script. Checked at t = 10, 50, and 100 s with a 2 % relative tolerance.

---

## How to run

```sh
cd tutorials/case312-preloaded-slab
./Allrun
```

The `Allrun` script executes in order:

1. `blockMesh` — creates the 1-D mesh with block label `slab`
2. `mkdir -p constant/slab` then moves `constant/polyMesh` into `constant/slab/polyMesh` — required by `foamMultiRun`'s region layout
3. `setFields -region slab` — imposes C = 1 in the box x ∈ [0, 0.01 m]
4. `foamMultiRun` — runs the `speciesSolid` solver

**Validate against the analytical solution:**

```sh
python3 validation/preloaded_compare.py
```

The script reads the `C_T2` field at t = 10, 50, and 100 s and compares to the Fourier cosine series. It reports the maximum relative error at each time and prints an overall `PASS` / `FAIL`.

# case318 — Membrane Diffusion Benchmark

**Source:** Hattab et al. 2024, *Fusion Eng. Des.* **202**, 114362, Section 3.3.1
**Solver:** `foamMultiRun` with `speciesSolid` region module

---

## What this case covers

A one-dimensional transient diffusion benchmark for hydrogen (H₂) permeating
through a thin metal membrane.  The case simultaneously validates against:

1. The **exact analytical solution** (Fourier sine series)
2. The **Pasler et al. OpenFOAM reference code** reported in the same paper

The geometry is a single 1.2 mm solid domain ("membrane") with constant
concentration boundary conditions — the simplest non-trivial test of the
species transport solver.

---

## Governing equation

$$\frac{\partial C}{\partial t} = D \frac{\partial^2 C}{\partial x^2}$$

| Symbol | Value | Description |
|--------|-------|-------------|
| $L$ | $1.2 \times 10^{-3}$ m | Membrane thickness |
| $D$ | $1.119018 \times 10^{-8}$ m²/s | Constant diffusivity at 300 K |
| $C_1$ | $0.83$ mol/m³ | Inlet (high-pressure side) concentration |
| $C_0$ | $0$ mol/m³ | Outlet (low-pressure side) concentration |
| $\tau_{SS}$ | $L^2/D \approx 128.7$ s | Steady-state diffusion time |

The domain is isothermal at 300 K and the activation energy is set to zero,
so $D$ is strictly constant throughout the simulation.

---

## Analytical solution

With $C(0,t) = C_1$, $C(L,t) = 0$, and $C(x,0) = 0$, the solution is:

$$C(x,t) = C_1\left(1 - \frac{x}{L}\right) - \frac{2C_1}{\pi} \sum_{n=1}^{\infty} \frac{1}{n} \sin\left(\frac{n\pi x}{L}\right) \exp\left[-\left(\frac{n\pi}{L}\right)^2 D\, t\right]$$

**Derivation sketch:**

Let $u = C - C_1(1-x/L)$ so that $u$ satisfies homogeneous Dirichlet BCs.
The eigenfunctions are $\sin(n\pi x/L)$ and the initial condition gives
coefficient $B_n = -2C_1/(n\pi)$ via the standard Fourier sine projection.

The steady state ($t \to \infty$) is the linear profile $C_{SS} = C_1(1-x/L)$,
reached on the time scale $\tau_{SS} = L^2/D$.

The simulation runs to $t = 30\,\text{s} \approx 0.23\,\tau_{SS}$, capturing
the early-to-mid transient when the diffusion front (depth $\sim 2\sqrt{Dt}$)
is approaching the far wall — the most sensitive regime for solver accuracy.

---

## Mesh

| Parameter | Value |
|-----------|-------|
| Geometry | 1-D slab, 1.2 mm × 1 mm × 1 mm |
| Cells | 90 in x, 1×1 in y,z (2-D empty planes) |
| Cell size | $\Delta x = 13.3\,\mu$m |
| Fourier number per step | $Fo = D\,\Delta t/\Delta x^2 \approx 0.32$ |

The 90-cell mesh matches the Pasler et al. reference discretisation.

---

## Setup

| File | Purpose |
|------|---------|
| `system/blockMeshDict` | 90-cell 1-D slab, region `membrane` |
| `system/controlDict` | `speciesSolid` solver, `endTime 30`, `deltaT 5e-3` |
| `system/membrane/fvSchemes` | Euler time, Gauss linear diffusion |
| `system/membrane/fvSolution` | PCG/DIC, tolerance 1e-10 |
| `constant/membrane/physicalProperties` | Metal solid thermo (isothermal proxy) |
| `constant/membrane/speciesProperties` | `C_H2`, $D = 1.119018\times10^{-8}$ m²/s, `noTrapping` |
| `0/membrane/C_H2` | IC: 0; inlet: 0.83 mol/m³; outlet: 0 |
| `0/membrane/T` | Uniform 300 K |

---

## How to run

```bash
cd tutorials/case318-membrane-benchmark
./Allrun
```

This runs `blockMesh`, moves the mesh into `constant/membrane/polyMesh`, and
launches `foamMultiRun`.  On a laptop it completes in under one minute.

### Validate against analytical solution

```bash
python3 validation/membrane_benchmark.py
```

Optional flags:

| Flag | Description |
|------|-------------|
| `--tol 0.02` | Relative tolerance for PASS/FAIL (default 2%) |
| `--times 1 5 10 20 30` | Override which time directories to compare |
| `--plot` | Save `membrane_benchmark.png` with profile overlays |

### Clean

```bash
./Allclean
```

---

## Expected output

At $t = 30\,\text{s}$ the L∞ relative error against the analytical solution
should be well below 0.5% (typically ~0.1% with 90 cells and $\Delta t = 5\,$ms).
A PASS is declared when the error at the last write time is less than 2%.

---

## Key references

- Hattab N., Stroth D., Pasler V. *et al.*, "multiSpeciesRegionFoam —
  A species-transport OpenFOAM solver for fusion applications,"
  *Fusion Eng. Des.* **202**, 114362 (2024).
  DOI: [10.1016/j.fusengdes.2024.114362](https://doi.org/10.1016/j.fusengdes.2024.114362)

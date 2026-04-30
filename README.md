# multiSpeciesRegionFoam

An OpenFOAM-13 addon for multi-region species transport through solid membranes and porous materials. The library models diffusion, Arrhenius-temperature-dependent transport coefficients, McNabb–Foster trapping kinetics, and thermodynamically consistent interface partitioning (Sieverts / Henry law). Typical applications include tritium permeation in fusion materials, dense-film membrane separation (reverse osmosis), and membrane distillation.

The solver plugs into OpenFOAM-13's modular `foamMultiRun` architecture: each solid region runs the `speciesSolid` solver module, which couples the species and energy equations within the same PIMPLE loop. Multi-region coupling (heat and species) is handled by OpenFOAM's existing mapped-patch infrastructure.

---

## Governing equations

### Species diffusion

In each solid region the mobile species concentration C [mol m⁻³] satisfies

```
∂C/∂t  =  ∇·(D(T) ∇C)  −  ∂Ct/∂t  +  S_vol
```

where D(T) is the temperature-dependent diffusivity (Arrhenius form, see below), Ct is the trapped concentration (zero when trapping is disabled), and S_vol [mol m⁻³ s⁻¹] is an optional uniform volumetric source.

### Heat conduction

The energy equation (solved by the inherited `solid` base module) is

```
∂(ρ e)/∂t  +  ∇·q  =  S_e
```

with q computed from `thermophysicalTransport` (Fourier conduction for a solid). The species and energy equations are solved sequentially within each PIMPLE corrector loop; the updated temperature is used to refresh D(T) before assembling the species equation.

### Arrhenius temperature dependence

All material properties (diffusivity D, solubility Ks, dissociation rate Kd, recombination rate Kr) can follow the Arrhenius form

```
X(T)  =  X₀ · exp(−Ea / (R · T))
```

where X₀ is the pre-exponential factor (carries the correct SI units), Ea is the activation energy [J mol⁻¹], and R = 8.314 J mol⁻¹ K⁻¹. Setting Ea = 0 gives a temperature-independent constant.

### McNabb–Foster trapping

When `trappingModel McNabbFoster` is selected, a trapped-species field Ct [mol m⁻³] evolves according to the single-trap ODE (Hattab et al. 2024, Eq. 2):

```
∂Ct/∂t  =  (αt / N) · C · (nt − Ct)  −  αd(T) · Ct
```

| Symbol | Meaning | Units |
|--------|---------|-------|
| N | molar trap-site density | mol m⁻³ |
| nt = ft · N | available trap-site density | mol m⁻³ |
| ft | trap fraction | — |
| αt | trapping rate coefficient | s⁻¹ |
| αd(T) = αd0 · exp(−ε / (kB T)) | detrapping rate (Arrhenius) | s⁻¹ |
| ε | trap binding energy | eV |

The ODE is discretised with a per-cell semi-implicit Euler step, which is unconditionally stable and bounded (0 ≤ Ct ≤ nt):

```
Ct^{n+1}  =  [Ct^n + Δt·(αt/N)·C·nt]  /  [1 + Δt·((αt/N)·C + αd)]
```

The right-hand side `∂Ct/∂t` is subtracted from the mobile species equation as a sink term.

### Interface partition (Sieverts law)

At a coupled interface between two diffusive regions the `sievertsCoupledMixed` boundary condition enforces:

1. **Flux continuity:** `D_s (∂C/∂n)_s = D_n (∂C/∂n)_n`
2. **Thermodynamic equilibrium (linear / Sieverts):** `C_s / Ks_s = C_n / Ks_n`

This is implemented as a mixed (Robin) BC. Denoting the conductance products `selfKD = D_s·δ_s·(Ks_s/Ks_n)` and `nbrKD = D_n·δ_n` (where δ is the inverse cell-to-face distance), the coefficients are

```
refValue  =  (Ks_s / Ks_n) · C_{n,cell}
refGrad   =  0
w         =  nbrKD / (nbrKD + selfKD)
```

so that `C_face = w · refValue + (1−w) · C_{s,cell}`, consistent with the standard conjugate-heat-transfer derivation applied to species transport.

A **quadratic (Picard-frozen)** variant, `partition quadratic`, models the case where one side is not in the Sieverts regime: `C_s = Ks_s · (C_n / Ks_n)²`. The nonlinearity is lagged to the previous PIMPLE outer iteration.

---

## Library architecture

The addon is split into three libraries and one solver module, built in dependency order:

```
src/
├── speciesTransport/    →  libspeciesTransport.so
│   ├── arrheniusProperty/   Temperature-dependent property X(T)=X₀exp(-Ea/RT)
│   ├── speciesModel/        Per-region owner of C, D(T), Ks, source, trapping
│   └── trappingModel/
│       ├── noTrapping           No-op (default for trap-free regions)
│       └── McNabbFoster         Single-trap model with semi-implicit ODE update
│
├── speciesCoupling/     →  libspeciesCoupling.so
│   ├── speciesCoupledMixed/     Abstract base: neighbour-field mapping
│   └── sievertsCoupledMixed/    Linear or quadratic Sieverts partition
│
└── speciesSolid/        →  libspeciesSolid.so  (solver module)
    └── speciesSolid             Extends the solid module with species transport
```

Each `controlDict` that uses this addon must load the required libraries:

```
libs  ("libspeciesTransport.so"  "libspeciesCoupling.so"  "libspeciesSolid.so");
```

---

## Build instructions

**Prerequisites**

- OpenFOAM 13 installed and sourced (e.g. `source /opt/openfoam13/etc/bashrc`)
- Standard C++ build tools (`g++`, `make`)

**Recommended location**

Place the repository under `$WM_PROJECT_USER_DIR`:

```
$WM_PROJECT_USER_DIR/
└── multiSpeciesRegionFoam/
    ├── Allwmake
    ├── Allwclean
    ├── README.md
    ├── src/
    └── tutorials/
```

**Build**

```sh
cd $WM_PROJECT_USER_DIR/multiSpeciesRegionFoam
./Allwmake -j4 2>&1 | tee build.log
```

Compilation takes under a minute. The three libraries and one solver object are installed into `$FOAM_USER_LIBBIN` and `$FOAM_USER_APPBIN` automatically by `wmake`.

**Clean**

```sh
./Allwclean
```

---

## Tutorials

All tutorials are self-contained cases in `tutorials/`. Each directory contains an `Allrun` script, an `Allclean` script, and a `validation/` subdirectory with a Python script that reads the OpenFOAM output and compares it to the analytical solution.

| Case | Physics | Key feature | Solver |
|------|---------|-------------|--------|
| [case311-slab-diffusion](tutorials/case311-slab-diffusion/README.md) | 1-D diffusion, step BC, erfc front | First end-to-end run; temperature field used as species proxy | `solid` |
| [case312-preloaded-slab](tutorials/case312-preloaded-slab/README.md) | 1-D diffusion, step IC, isolated slab | Fourier cosine-series solution; `setFields` for initial condition | `speciesSolid` |
| [case313-trapping](tutorials/case313-trapping/README.md) | Diffusion + McNabb–Foster trapping | Two trapping regimes (weak / strong); breakthrough-time scaling | `speciesSolid` |
| [case314-composite-membrane](tutorials/case314-composite-membrane/README.md) | Two-region diffusion, Sieverts interface | `sievertsCoupledMixed` BC; concentration jump at material interface | `speciesSolid` |
| [case315-membrane-distillation](tutorials/case315-membrane-distillation/README.md) | 1-D isothermal vapour diffusion | Membrane distillation with constant D; linear steady-state profile | `speciesSolid` |
| [case316-membrane-distillation-thermal](tutorials/case316-membrane-distillation-thermal/README.md) | Vapour diffusion with Arrhenius D(T) | Non-linear steady-state profile from temperature gradient; full energy coupling | `speciesSolid` |
| [case317-reverse-osmosis](tutorials/case317-reverse-osmosis/README.md) | Solution-diffusion across two regions | Henry (Sieverts) partition at feed/membrane interface; two-region coupling | `speciesSolid` |
| [case318-membrane-benchmark](tutorials/case318-membrane-benchmark/README.md) | 1-D transient diffusion, constant BCs | Code-to-code benchmark vs Pasler et al. and Fourier sine series; exact analytical solution | `speciesSolid` |
| [case319-permeation-barrier](tutorials/case319-permeation-barrier/README.md) | WC coating + SS316 tube, Sieverts interface | Real Arrhenius D(T) for both materials; PRF ≈ 200; with/without coating comparison | `speciesSolid` |

---

## References

- Hattab, N., Siriano, S., Giannetti, F. "An OpenFOAM multi-region solver for tritium transport modeling in fusion systems." *Fusion Engineering and Design* 202 (2024) 114362.
- McNabb, A., Foster, P. K. "A new analysis of the diffusion of hydrogen in iron and ferritic steels." *Trans. Metall. Soc. AIME* 227 (1963) 618.
- CFD Direct. "Modular Solvers in OpenFOAM." https://cfd.direct/openfoam/free-software/modular-solvers/

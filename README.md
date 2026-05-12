# multiSpeciesRegionFoam

An OpenFOAM-13 addon for multi-region species transport through solid membranes, fluid channels, and composite material systems. The library models advection-diffusion, Arrhenius temperature-dependent transport coefficients, McNabb–Foster trapping kinetics, and thermodynamically consistent interface partitioning (Sieverts / Henry law). Typical applications include tritium permeation in fission and fusion systems, dense-film membrane separation (reverse osmosis, membrane distillation), permeation barriers, and heat-exchanger conjugate transfer with dissolved species.

Three solver modules extend the OpenFOAM-13 `foamMultiRun` framework:

- **`speciesSolid`** — adds species transport to solid regions (or fluid regions with prescribed flow). Inherits the `solid` base module for temperature via `solidThermo`. Suitable when the velocity field is known or negligibly influences the solution.
- **`speciesFluid`** — adds species and temperature transport to fully solved incompressible fluid regions. Inherits the `incompressibleFluid` base module (PIMPLE U/p loop). Use when the velocity field must be computed from the Navier-Stokes equations.
- **`compressibleSpeciesFluid`** — adds species transport to compressible fluid regions. Inherits the OpenFOAM-13 `fluid` base module (compressible PIMPLE + full internal-energy equation). Required for gases, or when the standard `heRhoThermo` property infrastructure is preferred (use `rhoConst` EOS for liquids).

---

## Governing equations

### Species transport — solid regions (`speciesSolid`)

In each solid region the mobile species concentration C [mol m⁻³] satisfies:

```
∂C/∂t  =  ∇·(D(T) ∇C)  −  ∂Ct/∂t  +  S_vol
```

When a velocity field U is registered on the mesh (e.g. prescribed flow in a channel region using `speciesSolid`), an advection term is added automatically:

```
∂C/∂t  +  ∇·(U C)  =  ∇·(D(T) ∇C)  −  ∂Ct/∂t  +  S_vol
```

| Symbol | Meaning | Units |
|--------|---------|-------|
| C | Mobile species concentration | mol m⁻³ |
| D(T) | Temperature-dependent diffusivity (Arrhenius) | m² s⁻¹ |
| Ct | Trapped concentration (zero when trapping disabled) | mol m⁻³ |
| S_vol | Optional uniform volumetric source | mol m⁻³ s⁻¹ |
| U | Velocity field (if prescribed) | m s⁻¹ |

### Species, momentum, and energy transport — incompressible fluid regions (`speciesFluid`)

The `speciesFluid` module solves the coupled system of incompressible flow, scalar temperature, and species concentration. The full equation set solved in sequence within each PIMPLE outer corrector is:

**Continuity (incompressible):**
```
∇·U = 0
```

**Momentum (PIMPLE):**
```
∂U/∂t  +  ∇·(UU)  =  −∇(p/ρ)  +  ∇·(ν∇U)
```

**Energy (scalar transport form):**
```
∂T/∂t  +  ∇·(UT)  =  ∇·(α ∇T)      α = κ / (ρ Cp)
```

**Species (advection-diffusion):**
```
∂C/∂t  +  ∇·(UC)  =  ∇·(D(T) ∇C)  −  ∂Ct/∂t  +  S_vol
```

| Symbol | Meaning | Units |
|--------|---------|-------|
| U | Velocity vector | m s⁻¹ |
| p | Kinematic pressure (p/ρ) | m² s⁻² |
| ν | Kinematic viscosity | m² s⁻¹ |
| T | Temperature | K |
| α = κ/(ρCp) | Thermal diffusivity | m² s⁻¹ |
| κ | Thermal conductivity | W m⁻¹ K⁻¹ |
| ρ | Density | kg m⁻³ |
| Cp | Specific heat capacity | J kg⁻¹ K⁻¹ |
| C, D(T), Ct, S_vol | As above | — |

The PIMPLE loop first converges U and p (inherited from `incompressibleFluid`), then calls `thermophysicalPredictor()` which assembles and solves T followed by C. The updated T is used to refresh D(T) before solving the species equation.

> **Numerical note:** The advection terms `fvm::div(phi, T)` and `fvm::div(phi, C)` create **asymmetric** coefficient matrices. The solver for these fields must be `PBiCGStab` with `DILU` preconditioner. Using `PCG` (symmetric-only) will abort with `Unknown asymmetric matrix solver PCG`. Also, the `"T.*"` wildcard in `fvSolution` is required because the PIMPLE final corrector creates a `TFinal` field.

### Species transport in compressible fluid regions (`compressibleSpeciesFluid`)

The `compressibleSpeciesFluid` module wraps the OpenFOAM-13 `fluid` base class and adds molar concentration transport. The key difference from `speciesFluid` is the **mass flux conversion**: the `fluid` base supplies `phi` [kg/s] (mass flux), but the species equation needs the volumetric flux [m³/s]:

**Continuity (compressible):**
```
∂ρ/∂t  +  ∇·(ρU)  =  0
```

**Internal energy (solved for `e`, T recovered from EOS):**
```
∂(ρe)/∂t  +  ∇·(φe)  +  ∇·(φK)  +  ∇·(φ p/ρ)  =  ∇·(κ∇T)
```

**Species (advection-diffusion with volumetric flux):**
```
∂C/∂t  +  ∇·(U_vol C)  =  ∇·(D(T) ∇C)  −  ∂Ct/∂t  +  S_vol

    where  U_vol = phiU = phi / fvc::interpolate(rho)   [m³/s]
```

| Symbol | Meaning | Units |
|--------|---------|-------|
| φ = ρU·Sf | Face mass flux | kg/s |
| phiU = φ/ρ_face | Volumetric flux | m³/s |
| e | Sensible internal energy | J/kg |
| K = ½\|U\|² | Specific kinetic energy | J/kg |
| p | Thermodynamic pressure | Pa |

**Required `physicalProperties` format (`heRhoThermo`):**

```
thermoType
{
    type            heRhoThermo;  mixture  pureMixture;
    transport       const;        thermo   eConst;
    equationOfState rhoConst;     specie   specie;
    energy          sensibleInternalEnergy;
}
mixture
{
    specie          { molWeight  <M>; }
    equationOfState { rho  <rho>; }
    thermodynamics  { Cv  <Cv>;   Hf  0; }
    transport       { mu  <mu>;   Pr  <Pr>; }
}
```

Also requires `constant/<region>/thermophysicalTransport`:
```
laminar { model Fourier; }
```

**Required `fvSchemes` additions:**
```
divSchemes
{
    div(phi,e)          Gauss linearUpwind grad(e);
    div(phi,K)          Gauss linear;
    div(phi,(p|rho))    Gauss linear;
    div(phiU,C_<name>)  Gauss upwind;
}
```

**Required `fvSolution` settings:**
```
solvers
{
    "e.*"   { solver PBiCGStab; preconditioner DILU; tolerance 1e-10; relTol 0; }
    "C_.*"  { solver PBiCGStab; preconditioner DILU; tolerance 1e-10; relTol 0; }
    "rho.*" { solver diagonal; }
}
```

Pressure field uses thermodynamic dimensions `[1 -1 -2 0 0 0 0]` with a `fixedValue` outlet (e.g. 101325 Pa).

### Heat conduction — solid regions

The energy equation in `speciesSolid` regions is inherited from the OpenFOAM `solid` base module:

```
∂(ρ e)/∂t  +  ∇·q  =  S_e      q = −κ(T) ∇T
```

Properties are read from `physicalProperties` via `solidThermo` (thermoType block with `heSolidThermo`, `constIsoSolid`, etc.).

### Arrhenius temperature dependence

All material properties (D, Ks, Kd, Kr) can follow:

```
X(T)  =  X₀ · exp(−Ea / (R · T))
```

Setting Ea = 0 gives a temperature-independent constant.

```
D { X0  [0 2 -1 0 0 0 0]  5.08e-7;   Ea  71334; }   // Forcey 1988 SS316
```

### McNabb–Foster trapping

When `trappingModel McNabbFoster` is selected, trapped concentration Ct evolves as:

```
∂Ct/∂t  =  (αt / N) · C · (nt − Ct)  −  αd(T) · Ct
```

Discretised with a per-cell semi-implicit Euler step (unconditionally stable, bounded 0 ≤ Ct ≤ nt):

```
Ct^{n+1}  =  [Ct^n + Δt·(αt/N)·C·nt]  /  [1 + Δt·((αt/N)·C + αd)]
```

### Interface partition (Sieverts law)

At a coupled interface, `sievertsCoupledMixed` enforces flux continuity and thermodynamic equilibrium. Three partition modes:

| `partition` keyword | Self regime | Neighbour regime | Condition |
|---------------------|-------------|------------------|-----------|
| `linear` (default) | Sieverts | Sieverts | C_s/Ks_s = C_n/Ks_n |
| `quadratic` | Henry | Sieverts | C_s = Ks_s · (C_n/Ks_n)² |
| `sqrt` | Sieverts | Henry | C_s = Ks_s · √(C_n/Ks_n) |

The mixed (Robin) BC coefficients are:

```
refValue  =  (Ks_s / Ks_n) · C_n,cell           (linear)
w         =  nbrKD / (nbrKD + selfKD)
selfKD    =  D_s · δ_s⁻¹ · (Ks_s/Ks_n)          (linear)
selfKD    =  D_s · δ_s⁻¹ · C_s,face/(2·C_n,cell)  (sqrt — linearised)
selfKD    =  D_s · δ_s⁻¹ · 2·C_n,cell/C_s,face   (quadratic — linearised)
```

The linearised selfKD for nonlinear modes ensures flux continuity at every Picard step and avoids a 6–7% systematic error that occurs when using the constant Ks ratio.

### Surface recombination / dissociation

The `surfaceRecombination` BC models gas-phase atom/molecule exchange at a free surface:

```
D · ∂C/∂n  =  Kd(T) · p_gas  −  Kr(T) · C²
```

Setting p_gas = 0 models a vacuum or purge side.

### Membrane equilibrium BCs

**`antoineEquilibrium`** — sets C at a membrane face from the local temperature using the Antoine equation for vapour pressure:

```
ln(p_vap [Pa])  =  A  −  B / (C_ant + T)
C_face  =  p_vap / (R · T)   [mol/m³ ideal-gas]
```

**`latentHeatFlux`** — mixed BC for T at a fluid-membrane face that adds the latent heat of evaporation/condensation to the thermal balance:

```
q_latent  =  J · Lvap · M        [W/m²]
```

where J [mol/(m²·s)] is the species flux through the membrane, Lvap [J/mol] is the molar latent heat, and M [kg/mol] is the molar mass. The refGrad contribution to the T equation is ±q_latent/κ_fluid.

### Surface molar flux post-processing

The `speciesFlux` function object computes:

```
J  =  −D · ∇C · n̂    [mol/(m²·s)]
```

integrated over user-specified patches.

```
functions
{
    wallPermeation
    {
        type        speciesFlux;
        libs        ("libspeciesPost.so");
        region      wall;
        species     C_H2;
        patches     (wall_to_hitec);
        writeControl writeTime;
    }
}
```

Output: `postProcessing/<name>/<time>/speciesFlux.dat` with columns `time`, `<patch>_J`, `<patch>_total [mol/s]`.

---

## Solver selection guide

### Use `speciesSolid` for fluid regions when

The velocity field is **known in advance** (plug flow, analytically derived Poiseuille profile, or negligible flow effects). `speciesSolid` automatically detects a registered `U` field and adds `∇·(U C)` to the species equation and `∇·(ρU e)` to the energy equation. No U/p solve is performed.

**Prescribe the velocity:** create a `0/<region>/U` file with `fixedValue` at all boundaries and `internalField uniform (Ux Uy 0)`. Set `divSchemes { div(phi,C_H2) Gauss upwind; }` and no `div(phi,U)` entry (there is no momentum equation).

Appropriate when:
- Channel Re is low and flow development can be neglected (plug-flow approximation)
- The primary goal is species or heat transport; flow is a forcing term, not the unknown
- Computational cost must be minimised (no pressure solve or momentum iterations)
- Examples: case323 (FLiBe–SS316–Hitec HX, turbulent effective D absorbs flow effects), case325/326 (DCMD with plug-flow approximation for feed/permeate channels)

### Use `speciesFluid` for fluid regions when

The velocity and pressure fields must be **computed by solving the Navier-Stokes equations** and the fluid is **incompressible**. Required when:

- The velocity profile is unknown (developing flow, arbitrary geometry, pressure-driven flow with unknown parabolic/turbulent profile)
- The interaction between concentration gradients and flow (e.g. concentration polarisation driven by a developing boundary layer) must be captured accurately
- You need the momentum equation to provide a self-consistent pressure field (for multi-region cases where the fluid pressure is relevant)
- Buoyancy or other body forces affect the flow
- Examples: case327 (DCMD with solved laminar channel flow), case328 (Poiseuille + Graetz validation)

### Use `compressibleSpeciesFluid` for fluid regions when

The fluid is a **gas** (density varies with pressure), or when the OpenFOAM `heRhoThermo` property infrastructure is preferred — including for liquids using `rhoConst` EOS. Required when:

- The fluid is a gas at conditions where density variation matters
- You want to use the standard OpenFOAM `heRhoThermo` thermoType hierarchy (e.g. `perfectGas`, `rhoConst`, `Boussinesq`)
- The full internal-energy equation (with pressure-work term) should be solved rather than a scalar temperature equation
- Examples: case329 (FLiBe channel; `rhoConst` validates equivalence to incompressible formulation), case330 (He gas at 700 K, 1 atm)

**Required `physicalProperties` format for `speciesFluid`:**

```
viscosityModel  Newtonian;
nu              [0 2 -1 0 0 0 0] 4.74e-7;   // read by incompressibleFluid base

rho             983;        // kg/m³   — read by speciesFluid
Cp              4183;       // J/(kg·K)
mixture { transport { kappa 0.654; } }      // W/(m·K) — for latentHeatFlux BC
```

**Required `fvSolution` settings for `speciesFluid`:**

```
solvers
{
    "T.*"  { solver PBiCGStab; preconditioner DILU; tolerance 1e-10; relTol 0; }
    "C_.*" { solver PBiCGStab; preconditioner DILU; tolerance 1e-10; relTol 0; }
}
```

---

## Library architecture

The addon is split into five libraries and two solver modules, built in dependency order:

```
src/
├── speciesTransport/              →  libspeciesTransport.so
│   ├── arrheniusProperty/             X(T) = X₀·exp(-Ea/RT)
│   ├── speciesModel/                  Per-region owner of C, D(T), Ks, source, trapping
│   └── trappingModel/
│       ├── noTrapping                 No-op (default)
│       └── McNabbFoster               Single-trap model with semi-implicit ODE
│
├── speciesCoupling/               →  libspeciesCoupling.so
│   ├── sievertsCoupledMixed/          Sieverts/Henry interface partition (3 modes)
│   ├── surfaceRecombination/          Robin BC: D·∂C/∂n = Kd·p − Kr·C²
│   ├── antoineEquilibriumFvPatch/     fixedValue from Antoine equation at membrane face
│   └── latentHeatFluxFvPatch/         Mixed T BC adding q = J·Lvap·M
│
├── speciesSolid/                  →  libspeciesSolid.so  (solver module)
│   └── speciesSolid                   Extends solid: species + optional advection when U registered
│
├── speciesFluid/                  →  libspeciesFluid.so  (solver module)
│   └── speciesFluid                   Extends incompressibleFluid: PIMPLE + scalar T + C
│
├── compressibleSpeciesFluid/      →  libcompressibleSpeciesFluid.so  (solver module)
│   └── compressibleSpeciesFluid       Extends fluid: compressible PIMPLE + energy e + C
│                                      phiU = phi/rho_face  (mass→volumetric flux conversion)
│
└── speciesPost/                   →  libspeciesPost.so  (post-processing)
    └── speciesFlux                    Function object: ∫ -D∇C·n̂ dA over user patches
```

**Library loading in `controlDict`:**

```
libs  ("libspeciesTransport.so"  "libspeciesCoupling.so"
       "libspeciesSolid.so"      "libspeciesPost.so");
```

Add `"libspeciesFluid.so"` when any region uses `speciesFluid`.
Add `"libcompressibleSpeciesFluid.so"` when any region uses `compressibleSpeciesFluid`.

---

## Build instructions

**Prerequisites:** OpenFOAM 13 sourced (`source /opt/openfoam13/etc/bashrc`), standard C++ tools.

**Recommended location:**

```
$WM_PROJECT_USER_DIR/
└── multiSpeciesRegionFoam/
    ├── Allwmake
    ├── src/
    └── tutorials/
```

**Build all libraries:**

```sh
cd $WM_PROJECT_USER_DIR/multiSpeciesRegionFoam
./Allwmake -j4 2>&1 | tee build.log
```

Build order in `Allwmake`: `speciesTransport` → `speciesCoupling` → `speciesSolid` → `speciesFluid` → `compressibleSpeciesFluid` → `speciesPost`. Compilation takes under two minutes. All `.so` objects are installed into `$FOAM_USER_LIBBIN`.

**Clean:**

```sh
./Allwclean
```

---

## Tutorials

All tutorials are self-contained cases in `tutorials/`. Each contains `Allrun`, `Allclean`, and (for validation cases) a Python script in `validation/` that reads OpenFOAM output and compares to an analytical or benchmark solution. The script exits 0 on PASS and 1 on FAIL.

| Case | Physics summary | Key feature | Solver |
|------|----------------|-------------|--------|
| [case311-slab-diffusion](tutorials/case311-slab-diffusion/README.md) | 1-D diffusion, step BC, erfc front | Temperature field used as species proxy; validates pure-diffusion solver | `solid` |
| [case312-preloaded-slab](tutorials/case312-preloaded-slab/README.md) | 1-D diffusion, non-zero IC, isolated slab | Fourier cosine-series analytical solution; `setFields` for IC | `speciesSolid` |
| [case313-trapping](tutorials/case313-trapping/README.md) | Diffusion + McNabb–Foster trapping | Weak and strong trapping regimes; breakthrough-time scaling | `speciesSolid` |
| [case314-composite-membrane](tutorials/case314-composite-membrane/README.md) | Two-region diffusion, Sieverts interface | `sievertsCoupledMixed` linear partition; concentration jump at interface | `speciesSolid` |
| [case315-membrane-distillation](tutorials/case315-membrane-distillation/README.md) | Isothermal vapour diffusion through PTFE | Henry-law solubility; constant D; linear steady-state profile | `speciesSolid` |
| [case316-membrane-distillation-thermal](tutorials/case316-membrane-distillation-thermal/README.md) | Vapour diffusion with Arrhenius D(T) | Non-linear profile from T gradient; full energy–species coupling | `speciesSolid` |
| [case317-reverse-osmosis](tutorials/case317-reverse-osmosis/README.md) | Solution-diffusion, two regions | Henry partition at feed/membrane face (KH=0.5); two-region coupling | `speciesSolid` |
| [case318-membrane-benchmark](tutorials/case318-membrane-benchmark/README.md) | 1-D transient diffusion, constant BCs | Code-to-code benchmark vs Pasler et al.; Fourier sine-series exact solution | `speciesSolid` |
| [case319-permeation-barrier](tutorials/case319-permeation-barrier/README.md) | WC coating on SS316, real Arrhenius D(T) | PRF ≈ 200; with/without coating comparison | `speciesSolid` |
| [case320-shell-tube-hx](tutorials/case320-shell-tube-hx/README.md) | Three-region HX: FLiBe \| SS316 \| Hitec | Conjugate heat + H₂ permeation; Arrhenius D(T) from T gradient; `speciesFlux` FO | `speciesSolid` |
| [case321-henry-law](tutorials/case321-henry-law/README.md) | Henry vs Sieverts partition | `partition sqrt`/`quadratic`; C jumps UP at interface; 23.6% flux difference | `speciesSolid` |
| [case322-surface-recombination](tutorials/case322-surface-recombination/README.md) | Surface recombination kinetics | Robin BC; Da=1 gives Cs=(√5−1)/2; 38.2% flux reduction | `speciesSolid` |
| [case323-counterflow-hx](tutorials/case323-counterflow-hx/README.md) | 2-D counter-flow HX with H₂, prescribed flow | FLiBe–SS316–Hitec; turbulent effective D; advection via registered U in speciesSolid | `speciesSolid` |
| [case324-dcmd-khalifa2017](tutorials/case324-dcmd-khalifa2017/README.md) | DCMD validation vs Khalifa 2017 | Single-region membrane; J_sim=31.1 vs J_exp=35 L/(m²·h) | `speciesSolid` |
| [case325-dcmd-cfd](tutorials/case325-dcmd-cfd/README.md) | DCMD with explicit CFD temperature polarisation | Plug flow in fluid channels; `antoineEquilibrium` BC; three-region | `speciesSolid` |
| [case326-dcmd-latent-heat](tutorials/case326-dcmd-latent-heat/README.md) | DCMD with latent heat coupling | `latentHeatFlux` BC adds q=J·Lvap·M to T balance at membrane face | `speciesSolid` |
| [case327-dcmd-speciesFluid](tutorials/case327-dcmd-speciesFluid/README.md) | DCMD with fully solved incompressible flow | Feed/permeate use `speciesFluid` (PIMPLE U/p + T + C); membrane uses `speciesSolid` | `speciesFluid` + `speciesSolid` |
| [case328-poiseuille-graetz](tutorials/case328-poiseuille-graetz/README.md) | `speciesFluid` solver validation | Poiseuille profile (L2 < 2%) + Graetz Nu∞ = 7.5407 (< 5% error) | `speciesFluid` |
| [case329-flibe-channel](tutorials/case329-flibe-channel/README.md) | FLiBe (873 K) channel with H₂ source → SS316 permeation | heRhoThermo/rhoConst; compressible PIMPLE; Pe=1000; Sieverts interface; validates compressible = incompressible for constant-ρ | `compressibleSpeciesFluid` + `speciesSolid` |
| [case330-compressible-He-channel](tutorials/case330-compressible-He-channel/README.md) | He gas (700 K, 1 atm) channel with H₂ source → SS316 permeation | Same D, Ks, source, U as case329; cross-solver flux comparison verifies compressible solver | `compressibleSpeciesFluid` + `speciesSolid` |

---

## Quick-start example

```sh
# 1. Source OpenFOAM
source /opt/openfoam13/etc/bashrc

# 2. Build libraries
cd $WM_PROJECT_USER_DIR/multiSpeciesRegionFoam
./Allwmake -j4

# 3. Run a tutorial
cd tutorials/case314-composite-membrane
./Allrun

# 4. Validate
python3 validation/validate_composite_membrane.py
```

---

## References

- Hattab, N., Siriano, S., Giannetti, F. "An OpenFOAM multi-region solver for tritium transport modeling in fusion systems." *Fusion Engineering and Design* **202** (2024) 114362.
- McNabb, A., Foster, P. K. "A new analysis of the diffusion of hydrogen in iron and ferritic steels." *Trans. Metall. Soc. AIME* **227** (1963) 618–627.
- Forcey, K. S. et al. "Hydrogen transport and solubility in 316L and 1.4914 steels for fusion reactor applications." *J. Nucl. Mater.* **160** (1988) 153–159.
- Romatoski, R. R., Hu, L. W. "Fluoride salt coolant properties for nuclear reactor applications: A review." *Nucl. Technol.* **205** (2019) 1367–1388.
- Khalifa, A. et al. "Experimental and theoretical investigations on water desalination using direct contact membrane distillation." *Desalination* **404** (2017) 22–34.
- Baker, R. W. *Membrane Technology and Applications*, 3rd ed. Wiley, 2012.
- CFD Direct. "Modular Solvers in OpenFOAM." https://cfd.direct/openfoam/free-software/modular-solvers/

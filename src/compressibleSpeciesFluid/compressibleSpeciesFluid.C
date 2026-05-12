/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | multiSpeciesRegionFoam
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2026
     \\/     M anipulation  |
-------------------------------------------------------------------------------
\*---------------------------------------------------------------------------*/

#include "compressibleSpeciesFluid.H"
#include "surfaceInterpolate.H"
#include "fvmDdt.H"
#include "fvmDiv.H"
#include "fvmLaplacian.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(compressibleSpeciesFluid, 0);
    addToRunTimeSelectionTable(solver, compressibleSpeciesFluid, fvMesh);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::compressibleSpeciesFluid::compressibleSpeciesFluid
(
    fvMesh& mesh
)
:
    fluid(mesh),
    // Temperature is owned by fluidThermo (the fluid base class).
    // Pass a const reference to it; the speciesModel reads T but never writes it.
    species_(mesh, thermo_.T())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solvers::compressibleSpeciesFluid::~compressibleSpeciesFluid()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::compressibleSpeciesFluid::thermophysicalPredictor()
{
    // 1. Solve energy equation and update temperature via thermo_.correct().
    //    After this call, thermo_.T() and rho are fully updated for this
    //    PIMPLE outer corrector.
    fluid::thermophysicalPredictor();

    // 2. Refresh D(T) using the temperature just updated by the energy solver.
    species_.correctProperties();

    // 3. Derive volumetric flux from the compressible mass flux.
    //
    //    phi  [kg/s] = rho * U · Sf   (mass flux, from isothermalFluid base)
    //    phiU [m³/s] = phi / rho_face (volumetric flux = U · Sf)
    //
    //    The species equation governs molar concentration C [mol/m³], so the
    //    advection term is ∇·(U·C), not ∇·(rho·U·C).  Using phi directly
    //    would introduce a spurious factor of rho.
    const surfaceScalarField phiU
    (
        "phiU",
        phi / fvc::interpolate(rho)
    );

    // 4. Species equation:
    //    ∂C/∂t + ∇·(U·C) = ∇·(D(T)·∇C) − ∂Ct/∂t + S_vol
    volScalarField& C = species_.C();

    fvScalarMatrix CEqn
    (
        fvm::ddt(C)
      + fvm::div(phiU, C)
      - fvm::laplacian(species_.D(), C)
     ==
      - species_.trap().sink()
    );

    CEqn.source() +=
        static_cast<const scalarField&>(mesh_.V())
        * species_.source().value();

    CEqn.relax();
    fvConstraints().constrain(CEqn);
    CEqn.solve();
    fvConstraints().constrain(C);

    // 5. Update trapped concentration with the new C and current T.
    species_.trap().update(C, thermo_.T());
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | multiSpeciesRegionFoam
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2026
     \\/     M anipulation  |
-------------------------------------------------------------------------------
\*---------------------------------------------------------------------------*/

#include "speciesSolid.H"
#include "fvmDdt.H"
#include "fvmLaplacian.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(speciesSolid, 0);
    addToRunTimeSelectionTable(solver, speciesSolid, fvMesh);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::speciesSolid::speciesSolid(fvMesh& mesh)
:
    solid(mesh),
    species_(mesh, T_)       // T_ is the protected solidThermo temperature ref
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solvers::speciesSolid::~speciesSolid()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::speciesSolid::thermophysicalPredictor()
{
    // --- Thermal fields (from solid base) ---
    volScalarField& e   = thermo_.he();
    const volScalarField& rho = thermo_.rho();

    // --- Species field ---
    volScalarField& C = species_.C();

    while (pimple.correctNonOrthogonal())
    {
        // 1. Thermal energy equation (same as solid::thermophysicalPredictor)
        fvScalarMatrix eEqn
        (
            fvm::ddt(rho, e)
          + thermophysicalTransport->divq(e)
         ==
            fvModels().source(rho, e)
        );

        eEqn.relax();
        fvConstraints().constrain(eEqn);
        eEqn.solve();
        fvConstraints().constrain(e);
        thermo_.correct();

        // 2. Refresh diffusivity with the updated temperature
        species_.correctProperties();

        // 3. Species diffusion equation
        //    dC/dt = div(D(T)*grad(C)) - trap.sink() + source
        fvScalarMatrix CEqn
        (
            fvm::ddt(C)
          - fvm::laplacian(species_.D(), C)
         ==
          - species_.trap().sink()
        );

        // Add uniform volumetric source [mol/m³/s] * V [m³] → [mol/s per cell].
        // DimensionedField<scalar,volMesh> inherits from Field<scalar>, so the
        // static_cast gives us the underlying scalarField without using the
        // private OldTimeField::field() accessor.
        CEqn.source() +=
            static_cast<const scalarField&>(mesh_.V())
            * species_.source().value();

        CEqn.relax();
        fvConstraints().constrain(CEqn);
        CEqn.solve();
        fvConstraints().constrain(C);
    }

    // 4. Advance trapping ODEs once per PIMPLE outer corrector
    //    using the converged mobile concentration and temperature
    species_.trap().update(C, T_);
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | multiSpeciesRegionFoam
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2026
     \\/     M anipulation  |
-------------------------------------------------------------------------------
\*---------------------------------------------------------------------------*/

#include "speciesFluid.H"
#include "fvmDdt.H"
#include "fvmDiv.H"
#include "fvmLaplacian.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(speciesFluid, 0);
    addToRunTimeSelectionTable(solver, speciesFluid, fvMesh);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::speciesFluid::speciesFluid(fvMesh& mesh)
:
    incompressibleFluid(mesh),
    rho_(0),
    Cp_(0),
    kappa_(0),
    T_
    (
        IOobject
        (
            "T",
            runTime.name(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    species_(mesh, T_)
{
    // physicalProperties is registered as an IOdictionary by the viscosityModel
    // base class during incompressibleFluid construction.
    const IOdictionary& physProps =
        mesh.lookupObject<IOdictionary>("physicalProperties");

    rho_   = physProps.lookup<scalar>("rho");
    Cp_    = physProps.lookup<scalar>("Cp");
    kappa_ = physProps.subDict("mixture").subDict("transport")
                      .lookup<scalar>("kappa");
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solvers::speciesFluid::~speciesFluid()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::speciesFluid::thermophysicalPredictor()
{
    // Thermal diffusivity α = κ/(ρ·Cp)  [m²/s = dimKinematicViscosity]
    const dimensionedScalar alpha
    (
        "alpha",
        dimKinematicViscosity,
        kappa_ / (rho_ * Cp_)
    );

    volScalarField& C = species_.C();

    while (pimple.correctNonOrthogonal())
    {
        // 1. Temperature equation (incompressible constant-property form)
        //    ∂T/∂t + div(phi,T) = laplacian(α,T)
        fvScalarMatrix TEqn
        (
            fvm::ddt(T_)
          + fvm::div(phi_, T_)
          - fvm::laplacian(alpha, T_)
        );

        TEqn.relax();
        fvConstraints().constrain(TEqn);
        TEqn.solve();
        fvConstraints().constrain(T_);

        // 2. Refresh D(T) with updated temperature
        species_.correctProperties();

        // 3. Species equation
        //    ∂C/∂t + div(phi,C) = laplacian(D,C) - trap.sink() + source
        fvScalarMatrix CEqn
        (
            fvm::ddt(C)
          + fvm::div(phi_, C)
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
    }

    species_.trap().update(C, T_);
}


// ************************************************************************* //

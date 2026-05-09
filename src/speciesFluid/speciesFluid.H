/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | multiSpeciesRegionFoam
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2026
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Class
    Foam::solvers::speciesFluid

Description
    Solver module for single-species transport in incompressible fluid regions.

    Inherits the velocity-pressure PIMPLE loop from incompressibleFluid and
    adds:
      - a temperature field T solved via ∂T/∂t + div(phi,T) = laplacian(α,T)
      - a concentration field C_<species> via the speciesModel

    Constant thermal properties (rho, Cp, kappa) are read from
    constant/<region>/physicalProperties alongside the Newtonian viscosity
    entries required by the incompressibleFluid base class:

    \verbatim
        viscosityModel  Newtonian;
        nu              [0 2 -1 0 0 0 0] 4.74e-7;

        rho             983;        // kg/m^3
        Cp              4183;       // J/(kg.K)

        mixture { transport { kappa 0.654; } }   // for latentHeatFlux BC
    \endverbatim

    Species properties are read from constant/<region>/speciesProperties
    by the speciesModel.

SourceFiles
    speciesFluid.C

\*---------------------------------------------------------------------------*/

#ifndef speciesFluid_H
#define speciesFluid_H

#include "incompressibleFluid.H"
#include "speciesModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{

/*---------------------------------------------------------------------------*\
                       Class speciesFluid Declaration
\*---------------------------------------------------------------------------*/

class speciesFluid
:
    public incompressibleFluid
{
    // Constant thermal properties (read from physicalProperties)

        scalar rho_;
        scalar Cp_;
        scalar kappa_;


    // Temperature field (solved directly, not via solidThermo)

        volScalarField T_;


    // Species transport

        speciesModel species_;


public:

    //- Runtime type information
    TypeName("speciesFluid");


    // Constructors

        //- Construct from region mesh
        speciesFluid(fvMesh& mesh);

        //- Disallow default bitwise copy construction
        speciesFluid(const speciesFluid&) = delete;


    //- Destructor
    virtual ~speciesFluid();


    // Member Functions

        //- Solve the temperature and species equations
        virtual void thermophysicalPredictor();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const speciesFluid&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solvers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

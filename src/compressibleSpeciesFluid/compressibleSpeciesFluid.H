/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | multiSpeciesRegionFoam
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2026
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Class
    Foam::solvers::compressibleSpeciesFluid

Description
    Solver module for single-species transport in compressible fluid regions.

    Inherits the compressible PIMPLE loop (momentum, continuity, pressure
    correction) and the full energy equation from the OpenFOAM-13 \c fluid
    base class and adds:
      - a concentration field C_<species> solved via the speciesModel

    The species equation uses the volumetric flux derived from the mass flux
    to correctly advect molar concentration C [mol/m³]:

    \verbatim
        phiU = phi / fvc::interpolate(rho)     // [m³/s]
        ∂C/∂t + ∇·(phiU·C) = ∇·(D(T)·∇C) − ∂Ct/∂t + S_vol
    \endverbatim

    Using \c phi (mass flux [kg/s]) directly in the species equation would
    incorrectly introduce a factor of ρ in the advection term.  The phiU
    conversion is the critical difference from speciesFluid (incompressible).

    Temperature is not stored as an independent field; it is read from
    \c thermo_.T() which is kept up to date by the base class energy solver.

    Thermophysical properties are specified via the standard OpenFOAM
    \c thermoType hierarchy in \c constant/<region>/physicalProperties:

    \verbatim
        thermoType
        {
            type            heRhoThermo;
            mixture         pureMixture;
            transport       const;
            thermo          eConst;
            equationOfState rhoConst;    // or perfectGas for true compressibility
            specie          specie;
            energy          sensibleInternalEnergy;
        }
        mixture
        {
            specie      { molWeight <M>; }
            equationOfState { rho <rho>; }     // for rhoConst
            thermodynamics  { Cv <Cv>; Hf 0; }
            transport       { mu <mu>; Pr <Pr>; }
        }
    \endverbatim

    Required additional dictionary \c constant/<region>/thermophysicalTransport:

    \verbatim
        laminar { model Fourier; }
    \endverbatim

    Species properties are read from constant/<region>/speciesProperties
    by the speciesModel.

    \note fvSchemes must include:
    \verbatim
        div(phi,e)  Gauss linearUpwind grad(e);   // internal energy advection
        div(phi,K)  Gauss linear;                 // kinetic energy
    \endverbatim

    \note fvSolution must use "e.*" (not "T.*") for the energy field:
    \verbatim
        "e.*" { solver PBiCGStab; preconditioner DILU; tolerance 1e-10; relTol 0; }
    \endverbatim

SourceFiles
    compressibleSpeciesFluid.C

\*---------------------------------------------------------------------------*/

#ifndef compressibleSpeciesFluid_H
#define compressibleSpeciesFluid_H

#include "fluid.H"
#include "speciesModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{

/*---------------------------------------------------------------------------*\
                  Class compressibleSpeciesFluid Declaration
\*---------------------------------------------------------------------------*/

class compressibleSpeciesFluid
:
    public fluid
{
    // Species transport

        speciesModel species_;


public:

    //- Runtime type information
    TypeName("compressibleSpeciesFluid");


    // Constructors

        //- Construct from region mesh
        compressibleSpeciesFluid(fvMesh& mesh);

        //- Disallow default bitwise copy construction
        compressibleSpeciesFluid(const compressibleSpeciesFluid&) = delete;


    //- Destructor
    virtual ~compressibleSpeciesFluid();


    // Member Functions

        //- Solve energy (delegated to fluid base) then species equation
        virtual void thermophysicalPredictor();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const compressibleSpeciesFluid&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solvers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

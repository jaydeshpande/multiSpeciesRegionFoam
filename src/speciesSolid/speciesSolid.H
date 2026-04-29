/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | multiSpeciesRegionFoam
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2026
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Class
    Foam::solvers::speciesSolid

Description
    Solver module for tritium/hydrogen species transport in solid regions,
    derived from the OF-13 solid solver module.

    Inherits all thermal-solve behaviour from solid (heSolidThermo, coupled
    temperature BCs) and adds a species concentration equation:

        dC/dt = div(D(T) * grad(C)) - trap.sink() + source

    where D(T) is an Arrhenius diffusivity, trap.sink() is the volumetric
    trapping rate from the selected trappingModel (noTrapping or McNabb-Foster),
    and source is a uniform volumetric source (mol/m³/s).

    Species properties are read from constant/<region>/speciesProperties.
    The concentration field name is taken from the "species" keyword therein
    (e.g. "C_T2").  At region interfaces the sievertsCoupledMixed boundary
    condition enforces the Sieverts-law partition:
        C_self / Ks_self = C_nbr / Ks_nbr

    Usage in system/controlDict:
    \verbatim
        regionSolvers
        {
            pyc     speciesSolid;
            sic     speciesSolid;
        }
    \endverbatim

SourceFiles
    speciesSolid.C

\*---------------------------------------------------------------------------*/

#ifndef speciesSolid_H
#define speciesSolid_H

#include "solid.H"
#include "speciesModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{

/*---------------------------------------------------------------------------*\
                       Class speciesSolid Declaration
\*---------------------------------------------------------------------------*/

class speciesSolid
:
    public solid
{
    // Species transport

        //- Container for the concentration field, diffusivity, and trapping
        speciesModel species_;


public:

    //- Runtime type information
    TypeName("speciesSolid");


    // Constructors

        //- Construct from region mesh
        speciesSolid(fvMesh& mesh);

        //- Disallow default bitwise copy construction
        speciesSolid(const speciesSolid&) = delete;


    //- Destructor
    virtual ~speciesSolid();


    // Member Functions

        //- Construct and solve the energy and species equations,
        //  update thermophysical properties, refresh D(T)
        virtual void thermophysicalPredictor();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const speciesSolid&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solvers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

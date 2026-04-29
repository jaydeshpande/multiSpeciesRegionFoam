/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | multiSpeciesRegionFoam
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2026
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Class
    Foam::trappingModel

Description
    Abstract base class for models of hydrogen trapping in solids.

    A trapping model couples a mobile concentration field C (solved by the
    host solver) to one or more trapped-concentration fields Ct_i via a
    rate equation:

        dCt_i / dt = f_i(C, Ct_i, T)

    The model contributes a sink/source term to the mobile-C equation
    equal to - sum_i dCt_i/dt.

    Runtime selection:
        trappingModel
        {
            type            McNabbFoster;    // or noTrapping
            // ... model-specific entries ...
        }

SourceFiles
    trappingModel.C
    trappingModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef trappingModel_H
#define trappingModel_H

#include "volFields.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "autoPtr.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class trappingModel Declaration
\*---------------------------------------------------------------------------*/

class trappingModel
{
protected:

    // Protected data

        //- Reference to the mesh
        const fvMesh& mesh_;

        //- Model subdictionary
        dictionary dict_;

        //- Name of the mobile concentration field
        word cName_;


public:

    //- Runtime type information
    TypeName("trappingModel");


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            trappingModel,
            dictionary,
            (
                const fvMesh& mesh,
                const dictionary& dict,
                const word& cName
            ),
            (mesh, dict, cName)
        );


    // Constructors

        trappingModel
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const word& cName
        );

        //- Disallow copy
        trappingModel(const trappingModel&) = delete;
        void operator=(const trappingModel&) = delete;


    // Selectors

        //- Return a reference to the selected trapping model
        static autoPtr<trappingModel> New
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const word& cName
        );


    //- Destructor
    virtual ~trappingModel() = default;


    // Member Functions

        //- Return the sink term for the mobile-C equation: sum_i dCt_i/dt.
        //  Sign convention: returned field, when SUBTRACTED from the
        //  RHS of the mobile equation, gives the correct balance.
        //  I.e. mobile eq has (... = ... - sink)
        virtual tmp<volScalarField> sink() const = 0;

        //- Advance the trapped-species ODEs by the current timestep.
        //  Called once per PIMPLE outer iteration, AFTER the mobile-C
        //  solve, using the latest C field.
        virtual void update
        (
            const volScalarField& C,
            const volScalarField& T
        ) = 0;

        //- Return a const-ref to the i-th trapped field (for I/O & probes)
        virtual const volScalarField& Ct(const label i = 0) const = 0;

        //- Number of trap populations
        virtual label nTraps() const = 0;
};


} // End namespace Foam

#endif

// ************************************************************************* //

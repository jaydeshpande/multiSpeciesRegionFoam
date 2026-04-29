/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | multiSpeciesRegionFoam
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Class
    Foam::trappingModels::noTrapping

Description
    No-trapping model. Returns a zero sink term and does not allocate any
    trapped fields. Used in fluid regions and in solids where defects do
    not play a role (Hattab cases 3.1.1, 3.1.2, 3.1.4, 3.2, 3.3.x, 3.4.1).

SourceFiles
    noTrapping.C

\*---------------------------------------------------------------------------*/

#ifndef noTrapping_H
#define noTrapping_H

#include "trappingModel.H"

namespace Foam
{
namespace trappingModels
{

/*---------------------------------------------------------------------------*\
                         Class noTrapping Declaration
\*---------------------------------------------------------------------------*/

class noTrapping
:
    public trappingModel
{
    //- Dummy field returned from Ct(); zero, never written
    mutable autoPtr<volScalarField> dummyPtr_;


public:

    //- Runtime type information
    TypeName("noTrapping");


    // Constructors

        noTrapping
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const word& cName
        );


    //- Destructor
    virtual ~noTrapping() = default;


    // Member Functions

        virtual tmp<volScalarField> sink() const override;

        virtual void update
        (
            const volScalarField& C,
            const volScalarField& T
        ) override
        {}

        virtual const volScalarField& Ct(const label i = 0) const override;

        virtual label nTraps() const override { return 0; }
};


} // End namespace trappingModels
} // End namespace Foam

#endif

// ************************************************************************* //

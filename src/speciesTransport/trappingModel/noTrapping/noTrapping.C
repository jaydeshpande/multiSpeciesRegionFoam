/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | multiSpeciesRegionFoam
\*---------------------------------------------------------------------------*/

#include "noTrapping.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace trappingModels
{
    defineTypeNameAndDebug(noTrapping, 0);
    addToRunTimeSelectionTable(trappingModel, noTrapping, dictionary);
}
}


// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * //

Foam::trappingModels::noTrapping::noTrapping
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& cName
)
:
    trappingModel(mesh, dict, cName),
    dummyPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::trappingModels::noTrapping::sink() const
{
    return volScalarField::New
    (
        "noTrapping::sink",
        mesh_,
        dimensionedScalar
        (
            "zero",
            dimMoles/dimVolume/dimTime,   // sink for atomic/molar conc eqn
            0.0
        )
    );
}


const Foam::volScalarField&
Foam::trappingModels::noTrapping::Ct(const label) const
{
    if (!dummyPtr_.valid())
    {
        dummyPtr_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "Ct_dummy",
                    mesh_.time().name(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar("zero", dimMoles/dimVolume, 0.0)
            )
        );
    }
    return dummyPtr_();
}


// ************************************************************************* //

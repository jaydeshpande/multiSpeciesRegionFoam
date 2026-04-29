/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | multiSpeciesRegionFoam
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2026
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

#include "trappingModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(trappingModel, 0);
    defineRunTimeSelectionTable(trappingModel, dictionary);
}


// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * //

Foam::trappingModel::trappingModel
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& cName
)
:
    mesh_(mesh),
    dict_(dict),
    cName_(cName)
{}


// * * * * * * * * * * * * * * * * * Selector * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::trappingModel>
Foam::trappingModel::New
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& cName
)
{
    const word modelType(dict.lookup("type"));

    Info<< "Selecting trapping model: " << modelType << endl;

    auto cstrIter = dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalIOErrorInFunction(dict)
            << "Unknown trappingModel type " << modelType << nl << nl
            << "Valid trappingModel types are:" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalIOError);
    }

    return autoPtr<trappingModel>(cstrIter()(mesh, dict, cName));
}


// ************************************************************************* //

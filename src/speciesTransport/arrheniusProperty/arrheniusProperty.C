/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | multiSpeciesRegionFoam
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2026
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

#include "arrheniusProperty.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

const Foam::dimensionedScalar Foam::arrheniusProperty::R_
(
    "R",
    dimensionSet(1, 2, -2, -1, -1, 0, 0),  // kg*m^2 / (s^2 * K * mol) = J/(mol*K)
    8.314462618
);


// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * //

Foam::arrheniusProperty::arrheniusProperty
(
    const word& name,
    const dictionary& dict
)
:
    name_(name),
    X0_("X0", dict.lookup("X0")),                  // dimensions inferred from entry
    Ea_("Ea", dimensionSet(1, 2, -2, 0, -1, 0, 0), // J/mol
        dict.lookup<scalar>("Ea"))
{}


// * * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::arrheniusProperty::value(const volScalarField& T) const
{
    // X(T) = X0 * exp(-Ea / (R T))
    tmp<volScalarField> tX
    (
        volScalarField::New
        (
            name_,
            T.mesh(),
            X0_
        )
    );

    volScalarField& X = tX.ref();
    X.primitiveFieldRef() =
        X0_.value() * Foam::exp(-Ea_.value() / (R_.value() * T.primitiveField()));

    // Boundaries: evaluate on every patch to keep the field consistent
    forAll(X.boundaryField(), patchi)
    {
        X.boundaryFieldRef()[patchi] =
            X0_.value()
          * Foam::exp(-Ea_.value() / (R_.value() * T.boundaryField()[patchi]));
    }

    return tX;
}


Foam::tmp<Foam::scalarField>
Foam::arrheniusProperty::patchValue(const scalarField& Tp) const
{
    return tmp<scalarField>
    (
        new scalarField
        (
            X0_.value() * Foam::exp(-Ea_.value() / (R_.value() * Tp))
        )
    );
}


Foam::scalar Foam::arrheniusProperty::pointValue(const scalar T) const
{
    return X0_.value() * std::exp(-Ea_.value() / (R_.value() * T));
}


void Foam::arrheniusProperty::write(const word& keyword, Ostream& os) const
{
    os.writeKeyword(keyword);
    dictionary dict;
    dict.add("X0", X0_);
    dict.add("Ea", Ea_.value());
    os << dict << nl;
}


// ************************************************************************* //

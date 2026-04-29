/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | multiSpeciesRegionFoam
\*---------------------------------------------------------------------------*/

#include "speciesModel.H"

// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * //

Foam::speciesModel::speciesModel
(
    const fvMesh& mesh,
    const volScalarField& T
)
:
    mesh_(mesh),
    T_(T),
    properties_
    (
        IOobject
        (
            "speciesProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    name_(properties_.lookup("species")),
    C_
    (
        IOobject
        (
            name_,
            mesh.time().name(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    D_("D", properties_.subDict("D")),
    KsPtr_
    (
        properties_.found("Ks")
      ? new arrheniusProperty("Ks", properties_.subDict("Ks"))
      : nullptr
    ),
    KdPtr_
    (
        properties_.found("Kd")
      ? new arrheniusProperty("Kd", properties_.subDict("Kd"))
      : nullptr
    ),
    KrPtr_
    (
        properties_.found("Kr")
      ? new arrheniusProperty("Kr", properties_.subDict("Kr"))
      : nullptr
    ),
    source_
    (
        "source",
        dimMoles/dimVolume/dimTime,
        properties_.lookupOrDefault<scalar>("source", 0.0)
    ),
    trapPtr_
    (
        trappingModel::New
        (
            mesh,
            properties_.subDict("trappingModel"),
            name_
        )
    ),
    Dfield_
    (
        IOobject
        (
            "D_" + name_,
            mesh.time().name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimArea/dimTime, 0.0)
    )
{
    Info<< "Species " << name_ << " initialised in region " << mesh.name()
        << nl
        << "    D0 = " << D_.X0().value()
        << ", Ea_D = " << D_.Ea().value() << " J/mol"
        << nl;

    if (KsPtr_.valid())
    {
        Info<< "    Ks0 = " << KsPtr_->X0().value() << endl;
    }
    if (KdPtr_.valid())
    {
        Info<< "    Kd0 = " << KdPtr_->X0().value() << endl;
    }
    if (KrPtr_.valid())
    {
        Info<< "    Kr0 = " << KrPtr_->X0().value() << endl;
    }

    correctProperties();
}


// * * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * //

void Foam::speciesModel::correctProperties()
{
    Dfield_ = D_.value(T_);
}


// ************************************************************************* //

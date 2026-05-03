/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | multiSpeciesRegionFoam
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

#include "surfaceRecombinationFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(surfaceRecombinationFvPatchScalarField, 0);

    addToPatchFieldRunTimeSelection
    (
        fvPatchScalarField,
        surfaceRecombinationFvPatchScalarField
    );
}


// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * //

Foam::surfaceRecombinationFvPatchScalarField::
surfaceRecombinationFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF),
    Kd_(new arrheniusProperty("Kd", dict.subDict("Kd"))),
    Kr_(new arrheniusProperty("Kr", dict.subDict("Kr"))),
    pGas_(dict.lookup<scalar>("pGas"))
{
    if (dict.found("value"))
    {
        fvPatchScalarField::operator=
        (
            scalarField("value", iF.dimensions(), dict, p.size())
        );
    }
    else
    {
        fvPatchScalarField::operator=(patchInternalField());
    }

    gradient() = scalar(0);
}


Foam::surfaceRecombinationFvPatchScalarField::
surfaceRecombinationFvPatchScalarField
(
    const surfaceRecombinationFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
    Kd_(ptf.Kd_.valid() ? new arrheniusProperty(ptf.Kd_()) : nullptr),
    Kr_(ptf.Kr_.valid() ? new arrheniusProperty(ptf.Kr_()) : nullptr),
    pGas_(ptf.pGas_)
{}


Foam::surfaceRecombinationFvPatchScalarField::
surfaceRecombinationFvPatchScalarField
(
    const surfaceRecombinationFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(ptf, iF),
    Kd_(ptf.Kd_.valid() ? new arrheniusProperty(ptf.Kd_()) : nullptr),
    Kr_(ptf.Kr_.valid() ? new arrheniusProperty(ptf.Kr_()) : nullptr),
    pGas_(ptf.pGas_)
{}


// * * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * //

void Foam::surfaceRecombinationFvPatchScalarField::updateCoeffs()
{
    if (this->updated()) return;

    // Patch temperature: try T first, then fall back to 300 K
    const word Tname = "T";
    scalarField Tp(this->patch().size(), scalar(300));
    if (db().foundObject<volScalarField>(Tname))
    {
        Tp = db().lookupObject<volScalarField>(Tname)
                 .boundaryField()[patch().index()];
    }

    // Arrhenius rate constants at patch temperature
    const scalarField Kd = Kd_->patchValue(Tp);
    const scalarField Kr = Kr_->patchValue(Tp);

    // D at the patch face: look up D_<fieldName>
    const word Dname = "D_" + internalField().name();
    scalarField Dp(this->patch().size(), scalar(1));
    if (db().foundObject<volScalarField>(Dname))
    {
        Dp = db().lookupObject<volScalarField>(Dname)
                  .boundaryField()[patch().index()];
    }
    else
    {
        WarningInFunction
            << "Diffusivity field " << Dname << " not found; using 1 m²/s."
            << endl;
    }

    // Current cell-centre concentration (Picard freeze of C^2 nonlinearity)
    const scalarField Ccell = patchInternalField();

    // Robin gradient: D dC/dn = Kd*pGas - Kr*Ccell^2
    // => dC/dn = (Kd*pGas - Kr*Ccell^2) / D
    gradient() = (Kd * pGas_ - Kr * sqr(max(Ccell, scalar(0)))) / (Dp + SMALL);

    fixedGradientFvPatchScalarField::updateCoeffs();
}


void Foam::surfaceRecombinationFvPatchScalarField::write(Ostream& os) const
{
    fixedGradientFvPatchScalarField::write(os);

    if (Kd_.valid()) Kd_->write("Kd", os);
    if (Kr_.valid()) Kr_->write("Kr", os);
    writeEntry(os, "pGas", pGas_);
    writeEntry(os, "value", *this);
}


// ************************************************************************* //

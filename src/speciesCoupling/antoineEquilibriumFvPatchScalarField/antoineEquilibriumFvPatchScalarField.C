/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | multiSpeciesRegionFoam
\*---------------------------------------------------------------------------*/

#include "antoineEquilibriumFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(antoineEquilibriumFvPatchScalarField, 0);

    addToPatchFieldRunTimeSelection
    (
        fvPatchScalarField,
        antoineEquilibriumFvPatchScalarField
    );
}


// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * //

Foam::antoineEquilibriumFvPatchScalarField::
antoineEquilibriumFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict, false),
    A_(dict.lookup<scalar>("A")),
    B_(dict.lookup<scalar>("B")),
    C_(dict.lookup<scalar>("C")),
    TName_(dict.lookupOrDefault<word>("TName", "T")),
    mmHg_(dict.lookupOrDefault<scalar>("mmHg", scalar(133.322))),
    Rgas_(scalar(8.314))
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
        fvPatchScalarField::operator=(scalar(0));
    }
}


Foam::antoineEquilibriumFvPatchScalarField::
antoineEquilibriumFvPatchScalarField
(
    const antoineEquilibriumFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    A_(ptf.A_),
    B_(ptf.B_),
    C_(ptf.C_),
    TName_(ptf.TName_),
    mmHg_(ptf.mmHg_),
    Rgas_(ptf.Rgas_)
{}


Foam::antoineEquilibriumFvPatchScalarField::
antoineEquilibriumFvPatchScalarField
(
    const antoineEquilibriumFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(ptf, iF),
    A_(ptf.A_),
    B_(ptf.B_),
    C_(ptf.C_),
    TName_(ptf.TName_),
    mmHg_(ptf.mmHg_),
    Rgas_(ptf.Rgas_)
{}


// * * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * //

void Foam::antoineEquilibriumFvPatchScalarField::updateCoeffs()
{
    if (updated()) return;

    // Read the T face field from the same patch.
    // For membrane-fluid interfaces, coupledTemperature sets this T value
    // from the adjacent fluid region before the species equation is solved.
    const fvPatchScalarField& Tp =
        patch().lookupPatchField<volScalarField, scalar>(TName_);

    scalarField& Cface = *this;
    forAll(Cface, fi)
    {
        const scalar T_K    = Tp[fi];
        const scalar T_C    = T_K - scalar(273.15);
        const scalar log10p = A_ - B_ / (C_ + T_C);
        const scalar p_Pa   = mmHg_ * pow(scalar(10), log10p);
        Cface[fi]           = p_Pa / (Rgas_ * T_K);
    }

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::antoineEquilibriumFvPatchScalarField::write(Ostream& os) const
{
    fixedValueFvPatchScalarField::write(os);
    writeEntry(os, "A",     A_);
    writeEntry(os, "B",     B_);
    writeEntry(os, "C",     C_);
    writeEntryIfDifferent<word>(os, "TName", "T",      TName_);
    writeEntryIfDifferent<scalar>(os, "mmHg", 133.322, mmHg_);
}


// ************************************************************************* //

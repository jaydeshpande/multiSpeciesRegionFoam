/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | multiSpeciesRegionFoam
\*---------------------------------------------------------------------------*/

#include "speciesCoupledMixedFvPatchScalarField.H"
#include "mappedFvPatchBaseBase.H"
#include "mappedPatchBaseBase.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "fvPatchFieldMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(speciesCoupledMixedFvPatchScalarField, 0);
    addToPatchFieldRunTimeSelection
    (
        fvPatchScalarField,
        speciesCoupledMixedFvPatchScalarField
    );
}


// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * //

Foam::speciesCoupledMixedFvPatchScalarField::
speciesCoupledMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF, dict, false),
    nbrFieldName_(dict.lookupOrDefault<word>("nbrField", iF.name())),
    TName_(dict.lookupOrDefault<word>("T", "T")),
    nbrTName_(dict.lookupOrDefault<word>("Tnbr", "T")),
    propertiesName_
    (
        dict.lookupOrDefault<word>("properties", "speciesProperties")
    )
{
    // Validate that the underlying patch is a mapped type pointing to a
    // different patch (the coupled interface in another region).
    mappedPatchBaseBase::validateMapForField
    (
        *this,
        iF,
        dict,
        mappedPatchBaseBase::from::differentPatch
    );

    // Read the patch value field.
    fvPatchScalarField::operator=
    (
        scalarField("value", iF.dimensions(), dict, p.size())
    );

    // Initialise mixed-BC parts from restart data if present, otherwise
    // assume the field starts from the given value (zeroGradient-like).
    if (dict.found("refValue"))
    {
        refValue() =
            scalarField("refValue", iF.dimensions(), dict, p.size());
        refGrad() =
            scalarField
            (
                "refGradient",
                iF.dimensions()/dimLength,
                dict,
                p.size()
            );
        valueFraction() =
            scalarField("valueFraction", unitFraction, dict, p.size());
    }
    else
    {
        refValue() = *this;
        refGrad()  = 0.0;
        valueFraction() = 0.0;
    }
}


Foam::speciesCoupledMixedFvPatchScalarField::
speciesCoupledMixedFvPatchScalarField
(
    const speciesCoupledMixedFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    nbrFieldName_(ptf.nbrFieldName_),
    TName_(ptf.TName_),
    nbrTName_(ptf.nbrTName_),
    propertiesName_(ptf.propertiesName_)
{}


Foam::speciesCoupledMixedFvPatchScalarField::
speciesCoupledMixedFvPatchScalarField
(
    const speciesCoupledMixedFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    nbrFieldName_(ptf.nbrFieldName_),
    TName_(ptf.TName_),
    nbrTName_(ptf.nbrTName_),
    propertiesName_(ptf.propertiesName_)
{}


// * * * * * * * * * * * * * Protected Member Functions * * * * * * * * * * //

const Foam::fvPatch&
Foam::speciesCoupledMixedFvPatchScalarField::nbrPatch() const
{
    return mappedFvPatchBaseBase::getMap(patch()).nbrFvPatch();
}


Foam::tmp<Foam::scalarField>
Foam::speciesCoupledMixedFvPatchScalarField::nbrCellC() const
{
    const mappedFvPatchBaseBase& mapper =
        mappedFvPatchBaseBase::getMap(patch());
    const fvPatchScalarField& nbrFld =
        mapper.nbrFvPatch()
           .lookupPatchField<volScalarField, scalar>(nbrFieldName_);
    return mapper.fromNeighbour(nbrFld.patchInternalField());
}


Foam::tmp<Foam::scalarField>
Foam::speciesCoupledMixedFvPatchScalarField::nbrT() const
{
    const mappedFvPatchBaseBase& mapper =
        mappedFvPatchBaseBase::getMap(patch());
    const fvPatchScalarField& nbrTp =
        mapper.nbrFvPatch()
           .lookupPatchField<volScalarField, scalar>(nbrTName_);
    return mapper.fromNeighbour(nbrTp.patchInternalField());
}


Foam::tmp<Foam::scalarField>
Foam::speciesCoupledMixedFvPatchScalarField::nbrPatchT() const
{
    const mappedFvPatchBaseBase& mapper =
        mappedFvPatchBaseBase::getMap(patch());
    const fvPatchScalarField& nbrTp =
        mapper.nbrFvPatch()
           .lookupPatchField<volScalarField, scalar>(nbrTName_);
    return mapper.fromNeighbour(nbrTp);
}


Foam::tmp<Foam::scalarField>
Foam::speciesCoupledMixedFvPatchScalarField::thisPatchT() const
{
    const fvPatchScalarField& Tp =
        patch().lookupPatchField<volScalarField, scalar>(TName_);
    return tmp<scalarField>(new scalarField(Tp));
}


Foam::tmp<Foam::scalarField>
Foam::speciesCoupledMixedFvPatchScalarField::thisD() const
{
    const word DName = "D_" + this->internalField().name();
    const fvPatchScalarField& Dp =
        patch().lookupPatchField<volScalarField, scalar>(DName);
    return tmp<scalarField>(new scalarField(Dp));
}


Foam::tmp<Foam::scalarField>
Foam::speciesCoupledMixedFvPatchScalarField::nbrD() const
{
    const mappedFvPatchBaseBase& mapper =
        mappedFvPatchBaseBase::getMap(patch());
    const word DName = "D_" + nbrFieldName_;
    const fvPatchScalarField& nbrDp =
        mapper.nbrFvPatch()
           .lookupPatchField<volScalarField, scalar>(DName);
    return mapper.fromNeighbour(nbrDp);
}


Foam::tmp<Foam::scalarField>
Foam::speciesCoupledMixedFvPatchScalarField::nbrDeltaCoeffs() const
{
    const mappedFvPatchBaseBase& mapper =
        mappedFvPatchBaseBase::getMap(patch());
    return mapper.fromNeighbour(mapper.nbrFvPatch().deltaCoeffs());
}


// * * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * //

void Foam::speciesCoupledMixedFvPatchScalarField::write(Ostream& os) const
{
    mixedFvPatchScalarField::write(os);
    writeEntryIfDifferent<word>(os, "nbrField", this->internalField().name(),
                                nbrFieldName_);
    writeEntryIfDifferent<word>(os, "T", word("T"), TName_);
    writeEntryIfDifferent<word>(os, "Tnbr", word("T"), nbrTName_);
    writeEntryIfDifferent<word>(os, "properties", word("speciesProperties"),
                                propertiesName_);
}


// ************************************************************************* //

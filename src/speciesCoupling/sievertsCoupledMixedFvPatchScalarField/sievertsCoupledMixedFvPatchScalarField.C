/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | multiSpeciesRegionFoam
\*---------------------------------------------------------------------------*/

#include "sievertsCoupledMixedFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "UPstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(sievertsCoupledMixedFvPatchScalarField, 0);

    addToPatchFieldRunTimeSelection
    (
        fvPatchScalarField,
        sievertsCoupledMixedFvPatchScalarField
    );
}


// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * //

Foam::sievertsCoupledMixedFvPatchScalarField::
sievertsCoupledMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    speciesCoupledMixedFvPatchScalarField(p, iF, dict),
    partition_
    (
        dict.lookupOrDefault<word>("partition", "linear") == "quadratic"
      ? quadratic : linear
    ),
    Ks_(new arrheniusProperty("Ks", dict.subDict("Ks"))),
    KsNbr_(new arrheniusProperty("KsNbr", dict.subDict("KsNbr")))
{}


Foam::sievertsCoupledMixedFvPatchScalarField::
sievertsCoupledMixedFvPatchScalarField
(
    const sievertsCoupledMixedFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fieldMapper& mapper
)
:
    speciesCoupledMixedFvPatchScalarField(ptf, p, iF, mapper),
    partition_(ptf.partition_),
    Ks_(ptf.Ks_.valid() ? new arrheniusProperty(ptf.Ks_()) : nullptr),
    KsNbr_(ptf.KsNbr_.valid() ? new arrheniusProperty(ptf.KsNbr_()) : nullptr)
{}


Foam::sievertsCoupledMixedFvPatchScalarField::
sievertsCoupledMixedFvPatchScalarField
(
    const sievertsCoupledMixedFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    speciesCoupledMixedFvPatchScalarField(ptf, iF),
    partition_(ptf.partition_),
    Ks_(ptf.Ks_.valid() ? new arrheniusProperty(ptf.Ks_()) : nullptr),
    KsNbr_(ptf.KsNbr_.valid() ? new arrheniusProperty(ptf.KsNbr_()) : nullptr)
{}


// * * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * //

void Foam::sievertsCoupledMixedFvPatchScalarField::updateCoeffs()
{
    if (this->updated()) return;

    // Change comm tag to avoid collisions during parallel BC evaluation
    // (mirrors the pattern in coupledTemperatureFvPatchScalarField)
    int oldTag = UPstream::msgType();
    UPstream::msgType() = oldTag + 1;

    // ----- neighbour-side samples -----
    const scalarField nbrC     = nbrCellC();       // C_n,c mapped to this
    const scalarField nbrD_    = nbrD();            // D on neighbour patch
    const scalarField nbrDelta = nbrDeltaCoeffs();  // 1/dist on nbr patch
    const scalarField TnbrP    = nbrPatchT();       // neighbour patch T

    // ----- this-side samples -----
    const scalarField& selfDelta = patch().deltaCoeffs();
    const scalarField selfD_     = thisD();
    const scalarField TselfP     = thisPatchT();

    // ----- solubilities at the interface -----
    const scalarField Ks_self = Ks_->patchValue(TselfP);
    const scalarField Ks_nbr  = KsNbr_->patchValue(TnbrP);

    // ----- mixed-BC coefficients -----
    scalarField& refV = this->refValue();
    scalarField& refG = this->refGrad();
    scalarField& w    = this->valueFraction();

    refG = 0.0;

    // Conductance-weighted by the Sieverts partition ratio.
    // The partition ratio r = Ks_self/Ks_nbr enters on the self side because
    // flux continuity + Sieverts (C_s = r*C_n) gives:
    //   w = D_n*delta_n / (D_n*delta_n + r*D_s*delta_s)
    const scalarField selfKD = selfD_ * selfDelta * (Ks_self / Ks_nbr);
    const scalarField nbrKD  = nbrD_  * nbrDelta;

    w = nbrKD / (nbrKD + selfKD + SMALL);

    if (partition_ == linear)
    {
        // C_w,self = (Ks_self / Ks_nbr) * C_n,c
        refV = (Ks_self / Ks_nbr) * nbrC;
    }
    else
    {
        // Picard-frozen quadratic: C_s = Ks_s * (C_n / Ks_n)^2
        // The quadratic factor is evaluated at the current iterate;
        // PIMPLE outer iterations close the fixed-point.
        const scalarField ratio = nbrC / (Ks_nbr + SMALL);
        refV = Ks_self * ratio * ratio;
    }

    mixedFvPatchScalarField::updateCoeffs();

    UPstream::msgType() = oldTag;
}


void Foam::sievertsCoupledMixedFvPatchScalarField::write(Ostream& os) const
{
    speciesCoupledMixedFvPatchScalarField::write(os);

    writeEntryIfDifferent<word>
    (
        os,
        "partition",
        word("linear"),
        partition_ == quadratic ? word("quadratic") : word("linear")
    );

    if (Ks_.valid())    Ks_->write("Ks", os);
    if (KsNbr_.valid()) KsNbr_->write("KsNbr", os);
}


// ************************************************************************* //

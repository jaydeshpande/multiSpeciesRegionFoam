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
        [&]() -> partitionLaw
        {
            const word p = dict.lookupOrDefault<word>("partition", "linear");
            if (p == "quadratic") return quadratic;
            if (p == "sqrt")      return sqrt_;
            return linear;
        }()
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

    const scalarField nbrKD = nbrD_ * nbrDelta;

    if (partition_ == linear)
    {
        // C_s = (Ks_s/Ks_n)*C_n  →  dC_s/dC_n = Ks_s/Ks_n
        refV = (Ks_self / Ks_nbr) * nbrC;
        const scalarField selfKD = selfD_ * selfDelta * (Ks_self / Ks_nbr);
        w = nbrKD / (nbrKD + selfKD + SMALL);
    }
    else if (partition_ == quadratic)
    {
        // Picard-frozen quadratic: C_s = Ks_s*(C_n/Ks_n)^2
        // Self is Henry, neighbour is Sieverts.
        // dC_s/dC_n = 2*C_s/C_n  →  use linearised derivative for selfKD
        // so that the mixed-BC weight correctly enforces flux continuity.
        const scalarField ratio = nbrC / (Ks_nbr + SMALL);
        refV = Ks_self * ratio * ratio;
        const scalarField selfKD =
            selfD_ * selfDelta * 2.0 * refV / (nbrC + SMALL);
        w = nbrKD / (nbrKD + selfKD + SMALL);
    }
    else  // sqrt_: self is Sieverts, neighbour is Henry
    {
        // Picard-frozen sqrt: C_s = Ks_s*sqrt(C_n/Ks_n)
        // dC_s/dC_n = C_s/(2*C_n)  →  use linearised derivative for selfKD.
        refV = Ks_self * sqrt(max(nbrC, scalar(0)) / (Ks_nbr + SMALL));
        const scalarField selfKD =
            selfD_ * selfDelta * refV / (2.0 * max(nbrC, SMALL));
        w = nbrKD / (nbrKD + selfKD + SMALL);
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
        partition_ == quadratic ? word("quadratic")
      : partition_ == sqrt_    ? word("sqrt")
      : word("linear")
    );

    if (Ks_.valid())    Ks_->write("Ks", os);
    if (KsNbr_.valid()) KsNbr_->write("KsNbr", os);
}


// ************************************************************************* //

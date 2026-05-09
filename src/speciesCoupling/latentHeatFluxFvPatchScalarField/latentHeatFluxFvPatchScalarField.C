/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | multiSpeciesRegionFoam
\*---------------------------------------------------------------------------*/

#include "latentHeatFluxFvPatchScalarField.H"
#include "mappedFvPatchBaseBase.H"
#include "mappedPatchBaseBase.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "UPstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(latentHeatFluxFvPatchScalarField, 0);

    addToPatchFieldRunTimeSelection
    (
        fvPatchScalarField,
        latentHeatFluxFvPatchScalarField
    );
}


// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * //

Foam::latentHeatFluxFvPatchScalarField::latentHeatFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF, dict, false),
    Lvap_(dict.lookup<scalar>("Lvap")),
    M_(dict.lookup<scalar>("M")),
    CName_(dict.lookupOrDefault<word>("CName", "C_H2O")),
    evaporation_
    (
        [&]() -> bool
        {
            const word side = dict.lookup<word>("side");
            if (side == "evaporation") return true;
            if (side == "condensation") return false;
            FatalIOErrorInFunction(dict)
                << "side must be 'evaporation' or 'condensation', got: "
                << side << exit(FatalIOError);
            return true;
        }()
    ),
    selfMode_
    (
        [&]() -> bool
        {
            const word src = dict.lookupOrDefault<word>("source", "neighbour");
            if (src == "neighbour") return false;
            if (src == "self")      return true;
            FatalIOErrorInFunction(dict)
                << "source must be 'neighbour' or 'self', got: "
                << src << exit(FatalIOError);
            return false;
        }()
    )
{
    mappedPatchBaseBase::validateMapForField
    (
        *this,
        iF,
        dict,
        mappedPatchBaseBase::from::differentPatch
    );

    fvPatchScalarField::operator=
    (
        scalarField("value", iF.dimensions(), dict, p.size())
    );

    if (dict.found("refValue"))
    {
        refValue()       = scalarField("refValue", iF.dimensions(), dict, p.size());
        refGrad()        = scalarField("refGradient", iF.dimensions()/dimLength, dict, p.size());
        valueFraction()  = scalarField("valueFraction", unitFraction, dict, p.size());
    }
    else
    {
        refValue()      = *this;
        refGrad()       = scalar(0);
        valueFraction() = scalar(0);
    }
}


Foam::latentHeatFluxFvPatchScalarField::latentHeatFluxFvPatchScalarField
(
    const latentHeatFluxFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    Lvap_(ptf.Lvap_),
    M_(ptf.M_),
    CName_(ptf.CName_),
    evaporation_(ptf.evaporation_),
    selfMode_(ptf.selfMode_)
{}


Foam::latentHeatFluxFvPatchScalarField::latentHeatFluxFvPatchScalarField
(
    const latentHeatFluxFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    Lvap_(ptf.Lvap_),
    M_(ptf.M_),
    CName_(ptf.CName_),
    evaporation_(ptf.evaporation_),
    selfMode_(ptf.selfMode_)
{}


// * * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * //

void Foam::latentHeatFluxFvPatchScalarField::updateCoeffs()
{
    if (updated()) return;

    int oldTag = UPstream::msgType();
    UPstream::msgType() = oldTag + 1;

    // ── Access neighbour via mapped patch ────────────────────────────────────

    const mappedFvPatchBaseBase& mapper =
        mappedFvPatchBaseBase::getMap(patch());

    const fvPatch& nbrPatch = mapper.nbrFvPatch();
    const fvMesh&  nbrMesh  = nbrPatch.boundaryMesh().mesh();

    // ── Thermal conductivities ────────────────────────────────────────────

    const scalar kappaSelf =
        this->db().lookupObject<IOdictionary>("physicalProperties")
                  .subDict("mixture").subDict("transport")
                  .lookup<scalar>("kappa");

    const scalar kappaNbr =
        nbrMesh.lookupObject<IOdictionary>("physicalProperties")
               .subDict("mixture").subDict("transport")
               .lookup<scalar>("kappa");

    // ── Thermal coupling ─────────────────────────────────────────────────────

    // Neighbour T cell-centre values mapped to this patch
    const fvPatchScalarField& nbrTp =
        nbrPatch.lookupPatchField<volScalarField, scalar>("T");
    const scalarField nbrTcell = mapper.fromNeighbour(nbrTp.patchInternalField());

    // Mapped neighbour deltaCoeffs
    const scalarField nbrDelta = mapper.fromNeighbour(nbrPatch.deltaCoeffs());

    // Thermal conductance per unit area [W/m²/K]
    const scalarField selfKD = kappaSelf * patch().deltaCoeffs();
    const scalarField nbrKD  = kappaNbr  * nbrDelta;

    // ── Latent heat flux ─────────────────────────────────────────────────────

    const word DName = "D_" + CName_;
    scalarField J(patch().size(), scalar(0));

    if (selfMode_)
    {
        // Membrane side: read C and D from own patch (Antoine BC already applied)
        const fvPatchScalarField& selfCp =
            patch().lookupPatchField<volScalarField, scalar>(CName_);
        const scalarField selfCface = selfCp;                          // face (Antoine)
        const scalarField selfCcell = selfCp.patchInternalField();     // membrane cell
        const fvPatchScalarField& selfDp =
            patch().lookupPatchField<volScalarField, scalar>(DName);
        const scalarField selfDface = selfDp;

        if (evaporation_)
            J = max(selfDface * (selfCface - selfCcell) * patch().deltaCoeffs(), scalar(0));
        else
            J = max(selfDface * (selfCcell - selfCface) * patch().deltaCoeffs(), scalar(0));
    }
    else
    {
        // Fluid side: read C and D from neighbouring membrane patch
        const fvPatchScalarField& nbrCp =
            nbrPatch.lookupPatchField<volScalarField, scalar>(CName_);
        const scalarField nbrCface = mapper.fromNeighbour(nbrCp);
        const scalarField nbrCcell =
            mapper.fromNeighbour(nbrCp.patchInternalField());
        const fvPatchScalarField& nbrDp =
            nbrPatch.lookupPatchField<volScalarField, scalar>(DName);
        const scalarField nbrDface = mapper.fromNeighbour(nbrDp);

        if (evaporation_)
            J = max(nbrDface * (nbrCface - nbrCcell) * nbrDelta, scalar(0));
        else
            J = max(nbrDface * (nbrCcell - nbrCface) * nbrDelta, scalar(0));
    }

    const scalarField q = J * Lvap_ * M_;

    // ── Set mixed BC coefficients ─────────────────────────────────────────
    //
    // T_face = w·T_nbr_cell + (1-w)·T_self_cell + (1-w)·refGrad/selfDelta
    //
    //   w        = nbrKD / (selfKD + nbrKD)
    //   refValue = T_nbr_cell
    //   refGrad  = -q/kappa_self  (evaporation — heat sink on self)
    //   refGrad  = +q/kappa_self  (condensation — heat source on self)

    valueFraction() = nbrKD / (selfKD + nbrKD + SMALL);
    refValue()      = nbrTcell;

    if (evaporation_)
        refGrad() = -q / (kappaSelf + SMALL);
    else
        refGrad() = q / (kappaSelf + SMALL);

    mixedFvPatchScalarField::updateCoeffs();

    UPstream::msgType() = oldTag;
}


void Foam::latentHeatFluxFvPatchScalarField::write(Ostream& os) const
{
    mixedFvPatchScalarField::write(os);
    writeEntry(os, "Lvap", Lvap_);
    writeEntry(os, "M",    M_);
    writeEntryIfDifferent<word>(os, "CName", "C_H2O", CName_);
    writeEntry(os, "side", evaporation_ ? word("evaporation") : word("condensation"));
    if (selfMode_) writeEntry(os, "source", word("self"));
    writeEntry(os, "value", *this);
}


// ************************************************************************* //

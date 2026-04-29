/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | multiSpeciesRegionFoam
-------------------------------------------------------------------------------
Class
    Foam::speciesCoupledMixedFvPatchScalarField

Description
    Abstract base class for coupling a species concentration field across a
    fluid-solid or solid-solid region interface in foamMultiRun.

    The BC is implemented on top of OpenFOAM's mixedFvPatchField, which
    represents a wall value as

        phi_wall = w * refValue + (1-w) * (phi_cellCentre + delta * refGrad)

    where w = valueFraction, in [0,1].

    The subclass's responsibility is to supply:
      - refValue_: the "target" wall value implied by the interface jump
      - refGrad_:  the "target" wall normal gradient (usually 0 for
                   diffusion-limited cases; nonzero for surface-kinetics)
      - valueFraction_: how strongly the target value is enforced

    At construction the BC locates the conjugate patch in the neighbour
    region via OF-13's mappedFvPatchBaseBase (the same mechanism used by
    coupledTemperatureFvPatchScalarField for CHT).

    Nonlinear jump conditions (e.g. Sieverts' law involves Cs^2) are
    linearised across PIMPLE outer iterations using a Picard split: the
    previous outer-iterate value of the quadratic factor is held fixed
    when computing the current iterate. See Hattab et al. (2024) §2.2.

SourceFiles
    speciesCoupledMixedFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef speciesCoupledMixedFvPatchScalarField_H
#define speciesCoupledMixedFvPatchScalarField_H

#include "mixedFvPatchFields.H"
#include "mappedFvPatchBaseBase.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class speciesCoupledMixedFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class speciesCoupledMixedFvPatchScalarField
:
    public mixedFvPatchScalarField
{
protected:

    // Protected data

        //- Name of the species field in the neighbour region
        word nbrFieldName_;

        //- Name of the temperature field (this region), needed for
        //  evaluating Arrhenius properties on the patch.
        word TName_;

        //- Name of the temperature field in the neighbour region
        word nbrTName_;

        //- Name of the species-properties dictionary in this region
        //  (defaults to "speciesProperties")
        word propertiesName_;


    // Protected member functions

        //- Return the neighbour patch via mappedFvPatchBaseBase
        const fvPatch& nbrPatch() const;

        //- Get the neighbour's concentration cell-centre values, mapped to
        //  this patch
        tmp<scalarField> nbrCellC() const;

        //- Get the neighbour's temperature cell-centres, mapped to this patch
        tmp<scalarField> nbrT() const;

        //- Get the neighbour's patch-face temperature, mapped to this patch
        tmp<scalarField> nbrPatchT() const;

        //- Get THIS side's temperature on the patch
        tmp<scalarField> thisPatchT() const;

        //- Get THIS side's diffusivity on the patch (from "D_<speciesName>")
        tmp<scalarField> thisD() const;

        //- Get the neighbour's diffusivity on its patch, mapped to this patch
        //  (from "D_<nbrFieldName>" in the neighbour region)
        tmp<scalarField> nbrD() const;

        //- Get the neighbour's deltaCoeffs (1/dist), mapped to this patch
        tmp<scalarField> nbrDeltaCoeffs() const;


public:

    //- Runtime type information
    TypeName("speciesCoupledMixed");


    // Constructors

        //- Construct from patch, internal field and dictionary
        speciesCoupledMixedFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given field onto a new patch
        speciesCoupledMixedFvPatchScalarField
        (
            const speciesCoupledMixedFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fieldMapper&
        );

        //- Copy constructor setting internal field reference
        speciesCoupledMixedFvPatchScalarField
        (
            const speciesCoupledMixedFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new speciesCoupledMixedFvPatchScalarField(*this, iF)
            );
        }


    // Member Functions

        //- Write
        virtual void write(Ostream&) const;
};


} // End namespace Foam

#endif

// ************************************************************************* //

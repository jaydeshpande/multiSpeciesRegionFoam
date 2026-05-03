/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | multiSpeciesRegionFoam
-------------------------------------------------------------------------------
Class
    Foam::sievertsCoupledMixedFvPatchScalarField

Description
    Diffusion-limited interface coupling with a Sieverts-law partition.

    At a coupled interface between two diffusive regions (solid-solid or
    solid-fluid, as in Hattab et al. 2024 Eq. 8-9), two conditions hold:

      1. Flux continuity:    D_s (dC/dn)_s  =  D_n (dC/dn)_n
      2. Partition law (one of):
           linear    C_s/Ks_s = C_n/Ks_n          (Sieverts-Sieverts, default)
           quadratic C_s = Ks_s * (C_n/Ks_n)^2    (Henry self, Sieverts nbr)
           sqrt      C_s = Ks_s * sqrt(C_n/Ks_n)  (Sieverts self, Henry nbr)

    Cast as a mixed BC, the standard CHT derivation gives:

        refGrad  = 0
        refValue = partition-dependent (see below)
        w        = nbrKD / (nbrKD + selfKD)
        selfKD   = D_s * selfDeltaCoeff * (Ks_s / Ks_n)
        nbrKD    = D_n * nbrDeltaCoeff

    The quadratic and sqrt modes use Picard (fixed-point) linearisation:
    refValue is evaluated at the previous-iterate neighbour concentration
    and PIMPLE outer iterations drive the solution to the correct fixed point.

SourceFiles
    sievertsCoupledMixedFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef sievertsCoupledMixedFvPatchScalarField_H
#define sievertsCoupledMixedFvPatchScalarField_H

#include "speciesCoupledMixedFvPatchScalarField.H"
#include "arrheniusProperty.H"
#include "autoPtr.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
         Class sievertsCoupledMixedFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class sievertsCoupledMixedFvPatchScalarField
:
    public speciesCoupledMixedFvPatchScalarField
{
public:

    enum partitionLaw { linear, quadratic, sqrt_ };


private:

    // Private data

        //- Partition law (linear or quadratic)
        partitionLaw partition_;

        //- This side's solubility (Arrhenius)
        autoPtr<arrheniusProperty> Ks_;

        //- Neighbour side's solubility (Arrhenius)
        autoPtr<arrheniusProperty> KsNbr_;


public:

    //- Runtime type information
    TypeName("sievertsCoupledMixed");


    // Constructors

        //- Construct from patch, internal field and dictionary
        sievertsCoupledMixedFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given field onto a new patch
        sievertsCoupledMixedFvPatchScalarField
        (
            const sievertsCoupledMixedFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fieldMapper&
        );

        //- Copy constructor setting internal field reference
        sievertsCoupledMixedFvPatchScalarField
        (
            const sievertsCoupledMixedFvPatchScalarField&,
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
                new sievertsCoupledMixedFvPatchScalarField(*this, iF)
            );
        }


    // Member Functions

        virtual void updateCoeffs();
        virtual void write(Ostream&) const;
};


} // End namespace Foam

#endif

// ************************************************************************* //

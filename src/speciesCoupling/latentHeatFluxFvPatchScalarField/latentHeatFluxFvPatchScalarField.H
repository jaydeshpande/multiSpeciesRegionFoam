/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | multiSpeciesRegionFoam
-------------------------------------------------------------------------------
Class
    Foam::latentHeatFluxFvPatchScalarField

Description
    Mixed BC for the temperature field T at a membrane-fluid interface.
    Combines the thermal coupling of coupledTemperature with the latent heat
    flux J·L_vap·M due to phase change across the membrane.

    Can be applied on either side of the interface:

      source neighbour  (default — fluid side):
        C_H2O and D are read from the neighbouring membrane patch.
        Used on feed_to_membrane and permeate_to_membrane.

      source self  (membrane side):
        C_H2O and D are read from the patch's own region (the membrane).
        Used on membrane_to_feed and membrane_to_permeate in place of
        coupledTemperature, so that the latent heat correction is applied
        consistently on both sides of the interface.

    Physics
    -------
    Heat flux balance at the fluid-membrane interface:

        kappa_self · ∂T_self/∂n  =  kappa_nbr · ∂T_nbr/∂n  ±  q_latent

    where q_latent = J · L_vap · M  [W/m²].

    Mixed BC form (valueFraction w, refValue Tv, refGrad g):

        T_face = w · T_nbr_cell + (1-w) · T_self_cell + (1-w) · g / selfDelta

    where
        selfKD = kappa_self · selfDelta
        nbrKD  = kappa_nbr  · nbrDelta
        w      = nbrKD / (selfKD + nbrKD)
        g      = -q / kappa_self   (evaporation — heat sink)
        g      = +q / kappa_self   (condensation — heat source)

    Flux computation
    ----------------
    neighbour mode:
        Reads C_face (Antoine BC), C_cell, D from the neighbouring patch.
        J = D_nbr · (C_nbr_face − C_nbr_cell) · nbrDelta   (evaporation)
        J = D_nbr · (C_nbr_cell − C_nbr_face) · nbrDelta   (condensation)

    self mode:
        Reads C_face (Antoine BC), C_cell, D from own patch.
        J = D_self · (C_self_face − C_self_cell) · selfDelta  (evaporation)
        J = D_self · (C_self_cell − C_self_face) · selfDelta  (condensation)

    Dictionary format
    -----------------
    Fluid side (default):
    feed_to_membrane
    {
        type        latentHeatFlux;
        Lvap        2.45e6;         // latent heat of vaporisation [J/kg]
        M           0.018;          // molar mass of vapour species [kg/mol]
        CName       C_H2O;          // vapour concentration field name
        side        evaporation;    // "evaporation" or "condensation"
        value       uniform 333;
    }

    Membrane side:
    membrane_to_feed
    {
        type        latentHeatFlux;
        Lvap        2.45e6;
        M           0.018;
        CName       C_H2O;
        side        evaporation;
        source      self;           // read C/D from own patch
        value       uniform 333;
    }

SourceFiles
    latentHeatFluxFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef latentHeatFluxFvPatchScalarField_H
#define latentHeatFluxFvPatchScalarField_H

#include "mixedFvPatchFields.H"
#include "mappedPatchBaseBase.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class latentHeatFluxFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class latentHeatFluxFvPatchScalarField
:
    public mixedFvPatchScalarField
{
    // Private data

        scalar Lvap_;          // latent heat of vaporisation [J/kg]
        scalar M_;             // molar mass of vapour [kg/mol]
        word   CName_;         // concentration field name
        bool   evaporation_;   // true = heat sink (feed), false = heat source (permeate)
        bool   selfMode_;      // true = read C/D from own patch (membrane side)


public:

    //- Runtime type information
    TypeName("latentHeatFlux");


    // Constructors

        latentHeatFluxFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        latentHeatFluxFvPatchScalarField
        (
            const latentHeatFluxFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fieldMapper&
        );

        latentHeatFluxFvPatchScalarField
        (
            const latentHeatFluxFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new latentHeatFluxFvPatchScalarField(*this, iF)
            );
        }


    // Member Functions

        virtual void updateCoeffs();
        virtual void write(Ostream&) const;
};


} // End namespace Foam

#endif

// ************************************************************************* //

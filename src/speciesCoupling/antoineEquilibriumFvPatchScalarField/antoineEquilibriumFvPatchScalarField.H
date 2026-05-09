/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | multiSpeciesRegionFoam
-------------------------------------------------------------------------------
Class
    Foam::antoineEquilibriumFvPatchScalarField

Description
    Fixed-value BC for a vapour concentration field at a membrane face.
    The face value is set from the Antoine vapour-pressure equation evaluated
    at the local temperature:

        log10(p_sat / mmHg) = A − B / (C + T_Celsius)
        C_face = p_sat [Pa] / (R_gas · T [K])

    Intended for use on the membrane region at fluid-membrane interfaces
    where the adjacent fluid T is supplied via coupledTemperature.  The BC
    reads T from the same patch (lookupPatchField) so it picks up whatever
    value the T BC has set — no mapped-patch infrastructure required.

    Dictionary format
    -----------------
    patchName
    {
        type        antoineEquilibrium;
        A           8.10765;      // Antoine A coefficient (log10, mmHg scale)
        B           1750.286;     // Antoine B coefficient [°C]
        C           235.0;        // Antoine C offset [°C]
        TName       T;            // temperature field name (default "T")
        mmHg        133.322;      // Pa per mmHg (default 133.322)
        value       uniform 4.0;  // initial face value [mol/m³]
    }

SourceFiles
    antoineEquilibriumFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef antoineEquilibriumFvPatchScalarField_H
#define antoineEquilibriumFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class antoineEquilibriumFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class antoineEquilibriumFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private data

        scalar A_;
        scalar B_;
        scalar C_;          // Antoine C offset [°C]
        word   TName_;      // temperature field name
        scalar mmHg_;       // Pa per mmHg
        scalar Rgas_;       // universal gas constant [J/(mol·K)]


public:

    //- Runtime type information
    TypeName("antoineEquilibrium");


    // Constructors

        antoineEquilibriumFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        antoineEquilibriumFvPatchScalarField
        (
            const antoineEquilibriumFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fieldMapper&
        );

        antoineEquilibriumFvPatchScalarField
        (
            const antoineEquilibriumFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new antoineEquilibriumFvPatchScalarField(*this, iF)
            );
        }


    // Member Functions

        virtual void updateCoeffs();
        virtual void write(Ostream&) const;
};


} // End namespace Foam

#endif

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | multiSpeciesRegionFoam
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Class
    Foam::surfaceRecombinationFvPatchScalarField

Description
    Robin boundary condition for surface recombination/dissociation kinetics.

    Models the balance between gas-phase dissociation (supplying atomic species
    to the surface) and recombinative desorption (removing species as molecules):

        D * dC/dn  =  Kd(T) * pGas  -  Kr(T) * C_surface^2

    where
        D      = diffusivity   [m^2/s]      (looked up from D_<species>)
        Kd(T)  = dissociation rate          [mol/(m^2 s Pa)]  (Arrhenius)
        Kr(T)  = recombination rate         [m^4/(mol s)]     (Arrhenius)
        pGas   = imposed gas partial pressure [Pa]            (constant)
        C      = mobile species concentration [mol/m^3]

    The nonlinear C^2 term is Picard-frozen at the previous-iterate cell
    concentration; PIMPLE outer iterations close the fixed point.

    Dictionary entries:
    \verbatim
        type    surfaceRecombination;
        Kd
        {
            X0   <scalar>;   // pre-exponential [mol/(m^2 s Pa)]
            Ea   <scalar>;   // activation energy [J/mol]  (0 = constant)
        }
        Kr
        {
            X0   <scalar>;   // pre-exponential [m^4/(mol s)]
            Ea   <scalar>;
        }
        pGas   <scalar>;     // gas-phase partial pressure [Pa]
        value  uniform 0;    // required placeholder
    \endverbatim

SourceFiles
    surfaceRecombinationFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceRecombinationFvPatchScalarField_H
#define surfaceRecombinationFvPatchScalarField_H

#include "fixedGradientFvPatchFields.H"
#include "arrheniusProperty.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class surfaceRecombinationFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class surfaceRecombinationFvPatchScalarField
:
    public fixedGradientFvPatchScalarField
{
    // Private data

        //- Dissociation rate coefficient Kd(T)  [mol/(m^2 s Pa)]
        autoPtr<arrheniusProperty> Kd_;

        //- Recombination rate coefficient Kr(T)  [m^4/(mol s)]
        autoPtr<arrheniusProperty> Kr_;

        //- Gas-phase partial pressure [Pa]
        scalar pGas_;


public:

    //- Runtime type information
    TypeName("surfaceRecombination");


    // Constructors

        //- Construct from patch, internal field and dictionary
        surfaceRecombinationFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given field onto a new patch
        surfaceRecombinationFvPatchScalarField
        (
            const surfaceRecombinationFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fieldMapper&
        );

        //- Copy constructor setting internal field reference
        surfaceRecombinationFvPatchScalarField
        (
            const surfaceRecombinationFvPatchScalarField&,
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
                new surfaceRecombinationFvPatchScalarField(*this, iF)
            );
        }


    // Member Functions

        //- Update the gradient field from Kd, Kr, pGas and T
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

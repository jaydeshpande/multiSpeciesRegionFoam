/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | multiSpeciesRegionFoam
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2026
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    Free to use, modify and distribute.

Class
    Foam::arrheniusProperty

Description
    A temperature-dependent scalar property of Arrhenius form:

    \f[
        X(T) = X_0 \, \exp\!\left(\frac{-E_a}{R\,T}\right)
    \f]

    Used for diffusivity, solubility, dissociation and recombination
    coefficients in species transport. The pre-exponential factor X_0
    is a dimensioned scalar (its units carry the dimensions of the
    property, e.g. m^2/s for diffusivity), and E_a is in J/mol.

    Evaluation is deferred to a value(T) method, which takes a volume
    scalar temperature field and returns a volume scalar field of the
    property. A patchValue(T) overload is also provided for use inside
    boundary conditions, where only the patch temperature is available.

    Dictionary entry:
    \verbatim
        D
        {
            D0      1e-7;
            Ea      40000;   // J/mol
        }
    \endverbatim

SourceFiles
    arrheniusProperty.C

\*---------------------------------------------------------------------------*/

#ifndef arrheniusProperty_H
#define arrheniusProperty_H

#include "dimensionedScalar.H"
#include "volFields.H"
#include "dictionary.H"
#include "Ostream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class arrheniusProperty Declaration
\*---------------------------------------------------------------------------*/

class arrheniusProperty
{
    // Private data

        //- Name of the property (e.g. "D", "Ks")
        word name_;

        //- Pre-exponential factor (dimensioned)
        dimensionedScalar X0_;

        //- Activation energy [J/mol]
        dimensionedScalar Ea_;

        //- Universal gas constant [J/(mol*K)]
        static const dimensionedScalar R_;


public:

    // Constructors

        //- Construct from dictionary, picking up units from the X0 entry
        arrheniusProperty
        (
            const word& name,
            const dictionary& dict
        );


    //- Destructor
    ~arrheniusProperty() = default;


    // Member Functions

        //- Return the property name
        const word& name() const { return name_; }

        //- Return the pre-exponential factor
        const dimensionedScalar& X0() const { return X0_; }

        //- Return the activation energy
        const dimensionedScalar& Ea() const { return Ea_; }

        //- Evaluate over a volume temperature field
        tmp<volScalarField> value(const volScalarField& T) const;

        //- Evaluate on a patch (for use in boundary conditions)
        tmp<scalarField> patchValue(const scalarField& Tp) const;

        //- Evaluate at a single point (scalar T)
        scalar pointValue(const scalar T) const;

        //- Write as a named sub-dictionary entry (e.g. "D { X0 ...; Ea ...; }")
        void write(const word& keyword, Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | multiSpeciesRegionFoam
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

#ifndef speciesFlux_H
#define speciesFlux_H

#include "fvMeshFunctionObject.H"
#include "wordList.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                         Class speciesFlux Declaration
\*---------------------------------------------------------------------------*/

class speciesFlux
:
    public fvMeshFunctionObject
{
    // Private Data

        //- Species field name (e.g. "C_H2")
        word species_;

        //- Patch names to integrate over
        wordList patches_;

        //- Output file stream
        autoPtr<OFstream> filePtr_;

    // Private Member Functions

        //- Open (or re-open) the output file
        void openFile();

        //- No copy construct
        speciesFlux(const speciesFlux&) = delete;

        //- No copy assignment
        void operator=(const speciesFlux&) = delete;


public:

    //- Runtime type information
    TypeName("speciesFlux");


    // Constructors

        //- Construct from Time and dictionary
        speciesFlux
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~speciesFlux();


    // Member Functions

        //- Read the speciesFlux data
        virtual bool read(const dictionary&);

        //- Required by functionObject base class; this FO has no field deps
        virtual wordList fields() const { return wordList(); }

        //- Execute — nothing to do each timestep
        virtual bool execute();

        //- Compute and write fluxes at write time
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

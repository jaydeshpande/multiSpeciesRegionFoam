/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | multiSpeciesRegionFoam
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

#include "speciesFlux.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvPatch.H"
#include "OFstream.H"
#include "OSspecific.H"
#include "Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(speciesFlux, 0);
    addToRunTimeSelectionTable(functionObject, speciesFlux, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * //

void Foam::functionObjects::speciesFlux::openFile()
{
    if (!filePtr_.valid())
    {
        // postProcessing/<name>/<startTime>/speciesFlux.dat
        fileName outDir =
            mesh_.time().globalPath()
          / "postProcessing"
          / name()
          / mesh_.time().name();

        mkDir(outDir);

        filePtr_.reset(new OFstream(outDir / "speciesFlux.dat"));

        OFstream& os = filePtr_();
        os  << "# time";
        for (const word& pName : patches_)
        {
            os  << tab << pName + "_J[mol/(m2.s)]"
                << tab << pName + "_total[mol/s]";
        }
        os  << nl;
    }
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::functionObjects::speciesFlux::speciesFlux
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    species_(dict.lookup("species")),
    patches_(dict.lookup("patches"))
{
    read(dict);
}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

Foam::functionObjects::speciesFlux::~speciesFlux()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::speciesFlux::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    dict.lookup("species") >> species_;
    dict.lookup("patches") >> patches_;
    return true;
}


bool Foam::functionObjects::speciesFlux::execute()
{
    return true;
}


bool Foam::functionObjects::speciesFlux::write()
{
    // Look up concentration field and diffusivity field
    const word Dname = "D_" + species_;

    if (!mesh_.foundObject<volScalarField>(species_))
    {
        WarningInFunction
            << "Field " << species_ << " not found in mesh " << mesh_.name()
            << endl;
        return false;
    }
    if (!mesh_.foundObject<volScalarField>(Dname))
    {
        WarningInFunction
            << "Field " << Dname << " not found in mesh " << mesh_.name()
            << ". Ensure speciesModel AUTO_WRITE is enabled." << endl;
        return false;
    }

    const volScalarField& C = mesh_.lookupObject<volScalarField>(species_);
    const volScalarField& D = mesh_.lookupObject<volScalarField>(Dname);

    openFile();
    OFstream& os = filePtr_();
    os  << mesh_.time().value();

    const surfaceScalarField& magSf = mesh_.magSf();

    for (const word& pName : patches_)
    {
        label patchId = mesh_.boundaryMesh().findIndex(pName);
        if (patchId == -1)
        {
            WarningInFunction
                << "Patch " << pName << " not found in mesh " << mesh_.name()
                << endl;
            os  << tab << 0 << tab << 0;
            continue;
        }

        const scalarField& Dp = D.boundaryField()[patchId];
        // snGrad: (C_face - C_cell) / delta  — at fixedValue face this is
        // (C_bc - C_cell)*delta.  Use C.boundaryField().snGrad() which is
        // the outward normal gradient.
        const scalarField snGradC = C.boundaryField()[patchId].snGrad();
        const scalarField& Sfp = magSf.boundaryField()[patchId];

        // J = -D * snGrad(C)  [mol/(m²·s)]  (positive = into domain)
        scalarField Jface = -Dp * snGradC;

        scalar Jmean = gSum(Jface * Sfp) / (gSum(Sfp) + SMALL);
        scalar Jtotal = gSum(Jface * Sfp);

        os  << tab << Jmean << tab << Jtotal;

        Info<< name() << " [" << mesh_.name() << "] patch " << pName
            << "  J = " << Jmean << " mol/(m²·s)"
            << "  total = " << Jtotal << " mol/s"
            << nl;
    }

    os  << nl;
    return true;
}


// ************************************************************************* //

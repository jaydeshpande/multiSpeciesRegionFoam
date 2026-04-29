/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | multiSpeciesRegionFoam
\*---------------------------------------------------------------------------*/

#include "McNabbFoster.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace trappingModels
{
    defineTypeNameAndDebug(McNabbFoster, 0);
    addToRunTimeSelectionTable(trappingModel, McNabbFoster, dictionary);
}
}


// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * //

Foam::trappingModels::McNabbFoster::McNabbFoster
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& cName
)
:
    trappingModel(mesh, dict, cName),
    N_      (dict.lookup<scalar>("N")),
    ft_     (dict.lookup<scalar>("ft")),
    nt_     (ft_ * N_),
    alphat_ (dict.lookup<scalar>("alphat")),
    alphad0_(dict.lookup<scalar>("alphad0")),
    epsilon_eV_(dict.lookup<scalar>("epsilon")),
    Ct_
    (
        IOobject
        (
            "Ct_" + cName,
            mesh.time().name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimMoles/dimVolume, 0.0)
    )
{
    Info<< "    McNabbFoster: N = " << N_
        << ", ft = " << ft_
        << ", nt = " << nt_
        << ", alphat = " << alphat_
        << ", alphad0 = " << alphad0_
        << ", epsilon = " << epsilon_eV_ << " eV"
        << endl;
}


// * * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::trappingModels::McNabbFoster::sink() const
{
    // Report dCt/dt from the most recent update.
    // The mobile equation subtracts this, per Hattab Eq. 1.
    // We reconstruct it from (Ct - Ct.oldTime()) / dt, which is the
    // update's effective rate (consistent with semi-implicit Euler).
    const scalar dt = mesh_.time().deltaTValue();

    tmp<volScalarField> tS
    (
        volScalarField::New
        (
            "trapSink_" + cName_,
            (Ct_ - Ct_.oldTime()) / dimensionedScalar("dt", dimTime, dt)
        )
    );

    return tS;
}


void Foam::trappingModels::McNabbFoster::update
(
    const volScalarField& C,
    const volScalarField& T
)
{
    const scalar dt = mesh_.time().deltaTValue();

    // Update oldTime BEFORE modifying Ct so sink() reports this step's rate
    Ct_.oldTime();

    // Semi-implicit per-cell update:
    //   Ct_new = (Ct_old + dt*(alphat/N)*C*nt)
    //          / (1 + dt*((alphat/N)*C + alphad(T)))
    const scalar ainv = alphat_ / N_;

    scalarField& CtI = Ct_.primitiveFieldRef();
    const scalarField& CtOldI = Ct_.oldTime().primitiveField();
    const scalarField& CI = C.primitiveField();
    const scalarField& TI = T.primitiveField();

    forAll(CtI, celli)
    {
        const scalar ad = alphad(TI[celli]);
        const scalar num = CtOldI[celli] + dt * ainv * CI[celli] * nt_;
        const scalar den = 1.0 + dt * (ainv * CI[celli] + ad);
        CtI[celli] = num / den;
    }

    // Boundary patches: do the same (boundary values matter for visualization
    // and for any BC that might reference Ct; the physics is cell-based).
    forAll(Ct_.boundaryField(), patchi)
    {
        fvPatchScalarField& Ctp = Ct_.boundaryFieldRef()[patchi];
        const fvPatchScalarField& Cp = C.boundaryField()[patchi];
        const fvPatchScalarField& Tp = T.boundaryField()[patchi];
        const scalarField& CtpOld = Ct_.oldTime().boundaryField()[patchi];

        forAll(Ctp, facei)
        {
            const scalar ad = alphad(Tp[facei]);
            const scalar num = CtpOld[facei] + dt * ainv * Cp[facei] * nt_;
            const scalar den = 1.0 + dt * (ainv * Cp[facei] + ad);
            Ctp[facei] = num / den;
        }
    }

    Ct_.correctBoundaryConditions();
}


// ************************************************************************* //

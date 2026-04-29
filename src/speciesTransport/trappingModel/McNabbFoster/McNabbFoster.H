/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | multiSpeciesRegionFoam
-------------------------------------------------------------------------------
Class
    Foam::trappingModels::McNabbFoster

Description
    Single-trap McNabb-Foster trapping model (Hattab et al. 2024, Eq. 2):

    \f[
        \frac{\partial C_t}{\partial t}
        = \frac{\alpha_t}{N}\,C\,(n_t - C_t) - \alpha_d\,C_t
    \f]

    where
      - \f$ C   \f$ is the mobile concentration solved by the host
      - \f$ C_t \f$ is the trapped concentration (this field)
      - \f$ N   \f$ is the trap site molar density [mol/m^3]
      - \f$ n_t = f_t \, N \f$ is the available trap site density [mol/m^3]
      - \f$ \alpha_t \f$ is the trapping rate [1/s] (constant, from dict)
      - \f$ \alpha_d \f$ is the detrapping rate [1/s]
                      (Arrhenius form: alpha_d0 * exp(-eps / (k_B T)))

    Solved per-cell using a semi-implicit Euler step, implicit in Ct:

    \f[
      C_t^{n+1} = \frac{C_t^n + \Delta t\,(\alpha_t/N)\,C\,n_t}
                       {1 + \Delta t\,[(\alpha_t/N)\,C + \alpha_d]}
    \f]

    This is unconditionally stable for positive coefficients and bounded
    (0 <= Ct <= nt) in exact arithmetic.

Dictionary entries:
\verbatim
    trappingModel
    {
        type            McNabbFoster;
        N               1.053e5;      // trap site molar density [mol/m^3]  (= N_atomic / N_A)
        ft              0.1;          // trap fraction
        alphat          1e15;         // trapping rate [1/s]
        alphad0         1e13;         // detrapping pre-exp [1/s]
        epsilon         1.0;          // trap energy [eV]
    }
\endverbatim

SourceFiles
    McNabbFoster.C

\*---------------------------------------------------------------------------*/

#ifndef McNabbFoster_H
#define McNabbFoster_H

#include "trappingModel.H"

namespace Foam
{
namespace trappingModels
{

/*---------------------------------------------------------------------------*\
                       Class McNabbFoster Declaration
\*---------------------------------------------------------------------------*/

class McNabbFoster
:
    public trappingModel
{
    // Private data

        //- Total trap site molar density N [mol/m^3]
        scalar N_;

        //- Trap fraction [-]
        scalar ft_;

        //- Trap site density nt = ft * N [mol/m^3]
        scalar nt_;

        //- Trapping rate coefficient alpha_t [1/s]
        scalar alphat_;

        //- Detrapping pre-exponential [1/s]
        scalar alphad0_;

        //- Trap energy epsilon [eV]
        scalar epsilon_eV_;

        //- Boltzmann constant [eV/K]
        static constexpr scalar kB_eV_ = 8.617333262e-5;

        //- The trapped-species field (atomic volumetric concentration)
        volScalarField Ct_;


    // Private member functions

        //- Evaluate detrapping rate at a temperature [1/s]
        inline scalar alphad(const scalar T) const
        {
            return alphad0_ * std::exp(-epsilon_eV_ / (kB_eV_ * T));
        }


public:

    //- Runtime type information
    TypeName("McNabbFoster");


    // Constructors

        McNabbFoster
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const word& cName
        );


    //- Destructor
    virtual ~McNabbFoster() = default;


    // Member Functions

        virtual tmp<volScalarField> sink() const override;

        virtual void update
        (
            const volScalarField& C,
            const volScalarField& T
        ) override;

        virtual const volScalarField& Ct(const label i = 0) const override
        {
            return Ct_;
        }

        virtual label nTraps() const override { return 1; }
};


} // End namespace trappingModels
} // End namespace Foam

#endif

// ************************************************************************* //

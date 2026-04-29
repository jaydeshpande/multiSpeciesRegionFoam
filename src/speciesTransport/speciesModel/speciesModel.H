/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | multiSpeciesRegionFoam
-------------------------------------------------------------------------------
Class
    Foam::speciesModel

Description
    Per-region container for a single transported chemical species.

    Owns:
      - the volumetric concentration field C (mol/m^3)
      - the diffusivity (Arrhenius or constant)
      - optional solubility, dissociation, recombination rate (Arrhenius)
        needed by interface coupling BCs
      - a trapping model (may be noTrapping)
      - a uniform volumetric source term (optional)

    Reads its configuration from constant/<region>/speciesProperties:
    \verbatim
        species         T2;                         // name of species

        D
        {
            X0          [0 2 -1 0 0 0 0] 1e-7;
            Ea          40000;                      // J/mol
        }

        // optional -- only needed if an interface BC uses them
        Ks
        {
            X0          [0 -3 0 0 0 0 0] 4.8e-2;    // mol/m^3/Pa^0.5 for Sieverts
            Ea          8900;
        }

        Kd  { X0 ... ; Ea ... ; }
        Kr  { X0 ... ; Ea ... ; }

        source          0;                          // mol/m^3/s

        trappingModel
        {
            type        noTrapping;
        }
    \endverbatim

SourceFiles
    speciesModel.C

\*---------------------------------------------------------------------------*/

#ifndef speciesModel_H
#define speciesModel_H

#include "volFields.H"
#include "IOdictionary.H"
#include "arrheniusProperty.H"
#include "trappingModel.H"
#include "autoPtr.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class speciesModel Declaration
\*---------------------------------------------------------------------------*/

class speciesModel
{
    // Private data

        //- Reference to the mesh
        const fvMesh& mesh_;

        //- Reference to the region's temperature field
        //  (for Arrhenius evaluation)
        const volScalarField& T_;

        //- Properties dictionary
        IOdictionary properties_;

        //- Name of the species (e.g. "T2")
        word name_;

        //- Concentration field (the primary solved variable)
        volScalarField C_;

        //- Diffusivity property
        arrheniusProperty D_;

        //- Optional: solubility (present only if entry "Ks" exists)
        autoPtr<arrheniusProperty> KsPtr_;

        //- Optional: dissociation rate
        autoPtr<arrheniusProperty> KdPtr_;

        //- Optional: recombination rate
        autoPtr<arrheniusProperty> KrPtr_;

        //- Uniform volumetric source term [mol/m^3/s]
        dimensionedScalar source_;

        //- Trapping model
        autoPtr<trappingModel> trapPtr_;

        //- Cached diffusivity field, updated each timestep
        mutable volScalarField Dfield_;


public:

    // Constructors

        //- Construct from a mesh and a reference temperature field
        speciesModel
        (
            const fvMesh& mesh,
            const volScalarField& T
        );

        //- Disallow copy
        speciesModel(const speciesModel&) = delete;
        void operator=(const speciesModel&) = delete;


    //- Destructor
    ~speciesModel() = default;


    // Access

        const word& name() const { return name_; }
        const fvMesh& mesh() const { return mesh_; }
        const volScalarField& T() const { return T_; }

        volScalarField& C() { return C_; }
        const volScalarField& C() const { return C_; }

        //- Return the (cached) diffusivity field for the current T
        const volScalarField& D() const { return Dfield_; }

        //- Return the source term
        const dimensionedScalar& source() const { return source_; }

        //- Access to the trapping model
        const trappingModel& trap() const { return trapPtr_(); }
        trappingModel& trap() { return trapPtr_(); }

        //- Has each optional surface property?
        bool hasKs() const { return KsPtr_.valid(); }
        bool hasKd() const { return KdPtr_.valid(); }
        bool hasKr() const { return KrPtr_.valid(); }

        //- Access the optional surface properties (throws if absent)
        const arrheniusProperty& Ks() const { return KsPtr_(); }
        const arrheniusProperty& Kd() const { return KdPtr_(); }
        const arrheniusProperty& Kr() const { return KrPtr_(); }


    // Member Functions

        //- Refresh temperature-dependent fields (D, cached quantities).
        //  Call at the start of each PIMPLE outer iteration.
        void correctProperties();
};


} // End namespace Foam

#endif

// ************************************************************************* //

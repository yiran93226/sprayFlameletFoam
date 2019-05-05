#include "readGravitationalAcceleration.H"

Info<< "Reading thermophysical properties\n" << endl;
autoPtr<psiReactionThermo> pThermo(psiReactionThermo::New(mesh));
psiReactionThermo& thermo = pThermo();
thermo.validate(args.executable(), "h", "e");

SLGThermo slgThermo(mesh, thermo);

basicSpecieMixture& composition = thermo.composition();
PtrList<volScalarField>& Y = composition.Y();

const word defaultOxid("O2");
const word oxidizerSpecie(thermo.lookupOrDefault("oxidizerSpecie",defaultOxid));

volScalarField& p = thermo.p();

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    thermo.rho()
);

Info<< "\nReading field Z" << endl;
volScalarField Z
(
    IOobject
    (
        "Z",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field varZ" << endl;
volScalarField varZ
(
    IOobject
    (
        "varZ",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField chi
(
    IOobject
    (
        "chi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "\nReading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "compressibleCreatePhi.H"

mesh.setFluxRequired(p.name());

dimensionedScalar rhoMax
(
    dimensionedScalar::lookupOrDefault
    (
        "rhoMax",
        pimple.dict(),
        dimDensity,
        great
    )
);

dimensionedScalar rhoMin
(
    dimensionedScalar::lookupOrDefault
    (
        "rhoMin",
        pimple.dict(),
        dimDensity,
        0
    )
);

Info<< "Creating turbulence model\n" << endl;
autoPtr<compressible::turbulenceModel> turbulence
(
    compressible::turbulenceModel::New
    (
        rho,
        U,
        phi,
        thermo
    )
);

// Info<< "Creating combustion model\n" << endl;
// autoPtr<CombustionModel<psiReactionThermo>> combustion
// (
//     CombustionModel<psiReactionThermo>::New(thermo, turbulence())
// );

Info<< "Creating field dpdt\n" << endl;
volScalarField dpdt
(
    IOobject
    (
        "dpdt",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("dpdt", p.dimensions()/dimTime, 0)
);

Info<< "Creating field kinetic energy K\n" << endl;
volScalarField K("K", 0.5*magSqr(U));
multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

forAll(Y, i)
{
    fields.add(Y[i]);
}
fields.add(thermo.he());

// volScalarField Qdot
// (
//     IOobject
//     (
//         "Qdot",
//         runTime.timeName(),
//         mesh,
//         IOobject::NO_READ,
//         IOobject::AUTO_WRITE
//     ),
//     mesh,
//     dimensionedScalar("Qdot", dimEnergy/dimVolume/dimTime, 0.0)
// );

Info<< "Reading flamelet tables" << endl;
selector flameletTable("laminar");

#include "createMRF.H"
#include "createClouds.H"
#include "createRadiationModel.H"
#include "createFvOptions.H"

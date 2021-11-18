#include "createRDeltaT.H"

Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
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

#include "createPhi.H"


Info<< "Reading transportProperties\n" << endl;

immiscibleIncompressibleTwoPhaseMixture mixture(U, phi);

volScalarField& alpha1(mixture.alpha1());
volScalarField& alpha2(mixture.alpha2());

const dimensionedScalar& rho1 = mixture.rho1();
const dimensionedScalar& rho2 = mixture.rho2();


// Need to store rho for ddt(rho, U)
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT
    ),
    alpha1*rho1 + alpha2*rho2
);
rho.oldTime();


// Mass flux
surfaceScalarField rhoPhi
(
    IOobject
    (
        "rhoPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(rho)*phi
);


// Construct incompressible turbulence model
autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, mixture)
);


#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"


volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p_rgh + rho*gh
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    p_rgh,
    pimple.dict(),
    pRefCell,
    pRefValue
);

if (p_rgh.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
    p_rgh = p - rho*gh;
}

mesh.setFluxRequired(p_rgh.name());
mesh.setFluxRequired(alpha1.name());

#include "createMRF.H"
#include "createFvOptions.H"




// MULES correction
tmp<surfaceScalarField> talphaPhiCorr0;

Info<< "Reading concentration field A,B and C\n" << endl;

volScalarField A
(
	IOobject
	(
		"A",
		runTime.timeName(),
		mesh,
		IOobject::MUST_READ,
		IOobject::AUTO_WRITE	
	),
	mesh
);

volScalarField B
(
	IOobject
	(
		"B",
		runTime.timeName(),
		mesh,
		IOobject::MUST_READ,
		IOobject::AUTO_WRITE	
	),
	mesh
);

volScalarField C
(
	IOobject
	(
		"C",
		runTime.timeName(),
		mesh,
		IOobject::MUST_READ,
		IOobject::AUTO_WRITE	
	),
	mesh
);

/* volScalarField T
(
	IOobject
	(
		"T",
		runTime.timeName(),
		mesh,
		IOobject::MUST_READ,
		IOobject::AUTO_WRITE	
	),
	mesh
); */

volScalarField k_re
(
	IOobject
	(
		"k_re",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE	
	),
    C.mesh(),
	dimensionedScalar("k_re", dimensionSet(0, 3, -1, 0, -1, 0, 0), 0.0) // initial value will be overwritten in CEqn.H
);

IOdictionary transportProperties
(
	IOobject
	(
		"transportProperties",
		runTime.constant(),
		mesh,
		IOobject::MUST_READ_IF_MODIFIED,
		IOobject::NO_WRITE
	)
);

dimensionedScalar DT
(
	transportProperties.lookup ("DT")
);

dimensionedScalar C_thr
(
	transportProperties.lookup ("C_thr")
);

dimensionedScalar C_0
(
	transportProperties.lookup ("C_0")
);

dimensionedScalar He
(
	transportProperties.lookup("He")
);

dimensionedScalar Diff1
(
	transportProperties.lookup("Diff1")
);

dimensionedScalar Diff2
(
	transportProperties.lookup("Diff2")
);

dimensionedScalar k_re0
(
	transportProperties.lookup("k_re0")
);

dimensionedScalar beta
(
	transportProperties.lookup("beta")
);

dimensionedScalar k_pc
(
	transportProperties.lookup("k_pc")
);

dimensionedScalar k_T
(
	transportProperties.lookup("k_T")
);

volScalarField Diff ("Diff", Diff1* alpha1 + Diff2 *(1. - alpha1 ));




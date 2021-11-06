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

multiphaseMixture mixture(U, phi);

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
    mixture.rho()
);
rho.oldTime();

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
}

mesh.setFluxRequired(p_rgh.name());

#include "createMRF.H"
#include "createFvOptions.H"


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

dimensionedScalar cflux
(
	transportProperties.lookup ("cflux")
);

dimensionedScalar H
(
	transportProperties.lookup("H")
);

dimensionedScalar Diff1
(
	transportProperties.lookup("Diff1")
);

dimensionedScalar Diff2
(
	transportProperties.lookup("Diff2")
);

dimensionedScalar kappa
(
	transportProperties.lookup("kappa")
);

dimensionedScalar kC
(
	transportProperties.lookup("kC")
);

volScalarField alpha1= mixture.phases()["liquid1"];
//volScalarField& alpha1(mixture.alpha1());
//volScalarField alpha2= mixture.phases()["liquid2"];
volScalarField alpha3= mixture.phases()["oil1"];
volScalarField alpha4= mixture.phases()["gas1"];

//volScalarField alpha3= mixture.phases()["gas1"];

volScalarField Diff ("Diff", Diff1* alpha1);
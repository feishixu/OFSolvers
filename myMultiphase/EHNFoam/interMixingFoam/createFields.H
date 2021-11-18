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

immiscibleIncompressibleThreePhaseMixture mixture(U, phi);

volScalarField& alpha1(mixture.alpha1());
volScalarField& alpha2(mixture.alpha2());
volScalarField& alpha3(mixture.alpha3());

const dimensionedScalar& rho1 = mixture.rho1();
const dimensionedScalar& rho2 = mixture.rho2();
const dimensionedScalar& rho3 = mixture.rho3();

dimensionedScalar D23("D23", dimViscosity, mixture);

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
    alpha1*rho1 + alpha2*rho2 + alpha3*rho3
);
rho.oldTime();


// Mass flux
// Initialisation does not matter because rhoPhi is reset after the
// alpha solution before it is used in the U equation.
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
    rho1*phi
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
mesh.setFluxRequired(alpha2.name());

#include "createMRF.H"
#include "createFvOptions.H"

// Variable for species transport


/* Info<< "Reading concentration field A,B and C\n" << endl;

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
); */

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

dimensionedScalar He
(
	transportProperties.lookup("He")
);

dimensionedScalar Diff2
(
	transportProperties.lookup("Diff2")
);

dimensionedScalar Diff3
(
	transportProperties.lookup("Diff3")
);

dimensionedScalar k_re
(
	transportProperties.lookup("k_re")
);

dimensionedScalar k_pc
(
	transportProperties.lookup("k_pc")
);

dimensionedScalar spvalue
(
	transportProperties.lookup("spvalue")
);

volScalarField Diff ("Diff", Diff2* alpha2 + Diff3 *alpha3);

// ==========================================

/*     volScalarField Su1
    (
        IOobject
        (
            "Su1",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        //mesh
        (0.0 / rho1 *kC*C)
    );

    volScalarField Su2
    (
        IOobject
        (
            "Su2",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        //mesh
        (0.0 / rho2 *kC*C)
    ); */

/*     volScalarField Sp3
    (
        IOobject
        (
            "Sp3",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        (1.0 / rho1 *kC*C)/((1.0-alpha1) + SMALL)
    );

    volScalarField Sp4
    (
        IOobject
        (
            "Sp4",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        (1.0 / rho2 *kC*C)/(alpha2 + SMALL)
    ); */



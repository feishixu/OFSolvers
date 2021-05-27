// - C field
Diff = Diff_l* alpha1 + Diff_g *(1. - alpha1);

/*
volScalarField phiH
(
	"phiH",
	alpha1 * Diff_l * Diff_g / (alpha1 * Diff_g + (1 - alpha1) * Diff_l) * (1 - H) / (alpha1 * H + (1 - alpha1))
);
*/

volScalarField phiH
(
	"phiH",
	Diff*(H-1.)/(alpha1 + H*(1-alpha1))
	//rho*Diff*(H-1.)/((1.-alpha2) + H*alpha2)
);

surfaceScalarField phiHf
(
	"phiHf",
	fvc::interpolate(phiH)*fvc::snGrad(alpha2)*mesh.magSf()
);

fvScalarMatrix CEqn
(
			  fvm::ddt(C)
			+ fvm::div(phi, C)
			- fvm::laplacian(Diff, C)
			==
			- fvm::div(phiHf, C)
/*

			  fvm::ddt(rho, C)
			+ fvm::div(rhoPhi, C)
			- fvm::laplacian(rho*Diff, C)
			==
			- fvm::div(phiHf, C)
*/
);

CEqn.solve();

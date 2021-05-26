// - C field
Diff = Diff_l* alpha1 + Diff_g *(1. - alpha1);

volScalarField phiH
(
	"phiH",
	alpha1 * Diff_l * Diff_g / (alpha1 * Diff_g + (1 - alpha1) * Diff_l) * (1 - H) / (alpha1 * H + (1 - alpha1))
);

fvScalarMatrix AEqn
(
			  fvm::ddt(A)
			+ fvm::div(phi, A)
			- fvm::laplacian(Diff, A)
			==
			  fvc::laplacian(phiH, A)
			- kappa*A*B*alpha1
);

fvScalarMatrix BEqn
(
			  fvm::ddt(B)
			+ fvm::div(phi, B)
			- fvm::laplacian(Diff, B)
			==
			  fvc::laplacian(phiH, B)
			- kappa*A*B*alpha1
);

fvScalarMatrix CEqn
(
			  fvm::ddt(C)
			+ fvm::div(phi, C)
			- fvm::laplacian(Diff, C)
			==
			  fvc::laplacian(phiH, C)
			+ kappa*A*B*alpha1
);

AEqn.solve();
BEqn.solve();
CEqn.solve();




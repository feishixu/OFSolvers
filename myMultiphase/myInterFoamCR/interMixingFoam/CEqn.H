// - C field
Diff = Diff2* alpha2 + Diff3 *alpha3;

volScalarField phiH
(
	"phiH",
	alpha1 * Diff * (1 - H) / (alpha1 * H + (1 - alpha1))
);

fvScalarMatrix AEqn
(
			  fvm::ddt(A)
			+ fvm::div(phi, A)
			- fvm::laplacian(Diff, A)
			==
			  //fvc::laplacian(phiH, A)
			- kappa*A*B*(1-alpha1)
);

fvScalarMatrix BEqn
(
			  fvm::ddt(B)
			+ fvm::div(phi, B)
			- fvm::laplacian(Diff, B)
			==
			  //fvc::laplacian(phiH, B)
			- kappa*A*B*(1-alpha1)
);

fvScalarMatrix CEqn
(
			  fvm::ddt(C)
			+ fvm::div(phi, C)
			- fvm::laplacian(Diff, C)
			==
			  //fvc::laplacian(phiH, C)
			  kappa*A*B*(1-alpha1)
);

AEqn.solve();
BEqn.solve();
CEqn.solve();




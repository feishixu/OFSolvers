// - C field
D = D1* alpha1 + D2 *(1. - alpha1);

//D = D1 * D2 / (alpha1 * D1+ (1 - alpha1) * D2);

volScalarField phiH
(
	"phiH",
	 C * D * (1 - H) / (alpha1 + H*(1 - alpha1))
);

fvScalarMatrix CEqn
(
			  fvm::ddt(C)
			+ fvm::div(phi, C)
			- fvm::laplacian(D, C)
			==
			  fvc::laplacian(phiH, alpha1)
);

CEqn.solve();

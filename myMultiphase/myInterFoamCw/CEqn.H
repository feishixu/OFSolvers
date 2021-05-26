// - C field
Diff = Diff_l* alpha1 + Diff_g *(1. - alpha1);

volScalarField phiH
(
	"phiH",
	alpha1 * Diff_l * Diff_g / (alpha1 * Diff_g + (1 - alpha1) * Diff_l) * (1 - H) / (alpha1 * H + (1 - alpha1))
);

fvScalarMatrix CEqn
(
			  fvm::ddt(C)
			+ fvm::div(phi, C)
			- fvm::laplacian(Diff, C)
			==
			  fvc::laplacian(phiH, C)
);

CEqn.solve();

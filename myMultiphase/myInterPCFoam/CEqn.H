// - C field
Diff = Diff1* alpha1 + Diff2 *(1. - alpha1);

volScalarField phiH
(
	"phiH",
	alpha1 * Diff * (1 - He) / (alpha1 * He + (1 - alpha1))
);

fvScalarMatrix AEqn
(
			  fvm::ddt(A)
			+ fvm::div(phi, A)
			- fvm::laplacian(Diff, A)
			==
			  fvc::laplacian(phiH, A)
			- k_re*A*B*alpha1
);

fvScalarMatrix BEqn
(
			  fvm::ddt(B)
			+ fvm::div(phi, B)
			- fvm::laplacian(Diff, B)
			==
			  fvc::laplacian(phiH, B)
			- k_re*A*B*alpha1
);

fvScalarMatrix CEqn
(
			  fvm::ddt(C)
			+ fvm::div(phi, C)
			- fvm::laplacian(Diff, C)
			==
			  fvc::laplacian(phiH, C)
			+ k_re*A*B*alpha1
);

AEqn.solve();
BEqn.solve();
CEqn.solve();




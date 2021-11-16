// - C field
// Diff = Diff1 * alpha1 + Diff2 *(1. - alpha1);
Diff = Diff1 * Diff2 / (alpha1*Diff2 + (1 - alpha1)*Diff1);

/* volScalarField phiH
(
	"phiH",
	 Diff * (1 - He) / (alpha1 * He + (1 - alpha1))
	//-Diff * (1 - He) / (alpha1 + He * (1 - alpha1))

); */

scalar sign = 1.;
scalar H = He.value();

// if(He.value()>1.)
// {
// 	sign = -1.;
// 	H = 1./He.value();
// }

volScalarField phiHA
(
	"phiHA",
	-sign * A * Diff * (1.-H) / (alpha1 + H*(1.-alpha1))
);

volScalarField phiHB
(
	"phiHB",
	-sign * B * Diff * (1.-H) / (alpha1 + H*(1.-alpha1))
);

volScalarField phiHC
(
	"phiHC",
	-sign * C * Diff * (1.-H) / (alpha1 + H*(1.-alpha1))
);

fvScalarMatrix AEqn
(
			  fvm::ddt(A)
			+ fvm::div(phi, A)
			- fvm::laplacian(Diff, A)
			==
			  fvc::laplacian(phiHA, alpha1)
			- k_re*A*B*alpha1
);

fvScalarMatrix BEqn
(
			  fvm::ddt(B)
			+ fvm::div(phi, B)
			- fvm::laplacian(Diff, B)
			==
			  fvc::laplacian(phiHB, alpha1)
			- k_re*A*B*alpha1
);

fvScalarMatrix CEqn
(
			  fvm::ddt(C)
			+ fvm::div(phi, C)
			- fvm::laplacian(Diff, C)
			==
			  fvc::laplacian(phiHC, alpha1)
			+ k_re*A*B*alpha1
);

AEqn.solve();
BEqn.solve();
CEqn.solve();




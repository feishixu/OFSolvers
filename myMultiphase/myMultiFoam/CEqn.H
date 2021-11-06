// - C field
//Diff = Diff1* alpha1 + Diff2 *alpha2;
//Diff = Diff1 * Diff2 / (alpha1 * Diff2 + alpha2 * Diff1);
Diff = Diff1* alpha1;
/*
volScalarField phiH
(
	"phiH",
	alpha1 * Diff * (1 - H) / (alpha1 * H + alpha2)
	
);
*/

alpha1= mixture.phases()["liquid1"];

fvScalarMatrix AEqn
(
			  fvm::ddt(A)
			+ fvm::div(phi, A)
			//- fvm::laplacian(Diff, A)
			- fvm::laplacian(Diff*alpha1, A)
			==
			  //fvc::laplacian(phiH, A)
			- kappa*A*B*alpha1
);

fvScalarMatrix BEqn
(
			  fvm::ddt(B)
			+ fvm::div(phi, B)
			//- fvm::laplacian(Diff, B)
			- fvm::laplacian(Diff*alpha1, B)			
			==
			  //fvc::laplacian(phiH, B)
			- kappa*A*B*alpha1
);

fvScalarMatrix CEqn
(
			  fvm::ddt(C)
			+ fvm::div(phi, C)
			//- fvm::laplacian(Diff, C)
			- fvm::laplacian(Diff*alpha1, C)			
			==
			  //fvc::laplacian(phiH, C)
			  kappa*A*B*alpha1
);

AEqn.solve();
BEqn.solve();
CEqn.solve();




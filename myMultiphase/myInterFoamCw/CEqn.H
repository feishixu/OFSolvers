// - C field
//D = D1*alpha1 + D2*(1.-alpha1);

D = D1 * D2 / (alpha1 * D2+ (1 - alpha1) * D1);
scalar sign = 1.;
scalar He = H.value();

if(H.value()>1.)
{
	sign = -1.;
	He = 1./H.value();
}


volScalarField phiH
(
	"phiH",
	-sign * C * D * (1.-He) / (alpha1*He+(1.-alpha1))
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

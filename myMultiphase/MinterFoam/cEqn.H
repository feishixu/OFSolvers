//const surfaceScalarField& fvc::interpolate(alpha1)(fvc::interpolate(alpha1));
surfaceScalarField phiCi
(
 	"phiCi",
 	((fvc::interpolate(D)*
	  (1-2*H-fvc::interpolate(alpha1)+H*fvc::interpolate(alpha1))
	  /(fvc::interpolate(alpha1)+(1-fvc::interpolate(alpha1))*H + (1-fvc::interpolate(alpha1))*fvc::interpolate(alpha1)+(1-fvc::interpolate(alpha1))*(1-fvc::interpolate(alpha1))*H))
	  *fvc::snGrad(alpha1))*mesh.magSf()
);



fvScalarMatrix CEqn
(
 	fvm::ddt(C)+fvm::div(phi,C)-fvm::laplacian(fvc::interpolate(D),C) + fvc::div(phiCi,C)
);

CEqn.solve();




//		  + (1-fvc::interpolate(alpha1))*fvc::interpolate(alpha1)*(1-H)))





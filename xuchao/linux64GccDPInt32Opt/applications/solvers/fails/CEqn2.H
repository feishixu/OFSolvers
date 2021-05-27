surfaceScalarField phiCi
(
 	"phiCi",
 	((fvc::interpolate(D)*(1-H)/(fvc::interpolate(alpha1)+(1-fvc::interpolate(alpha1))*H))*fvc::snGrad(alpha1))*mesh.magSf()
);



fvScalarMatrix CEqn
(
 	fvm::ddt(C)+fvm::div(phi,C)-fvm::laplacian(fvc::interpolate(D),C)+fvm::div(phiCi, C)
);

CEqn.solve();









surfaceScalarField phipsngF
(
 "phipsngF",
 fvc::interpolate(D)*(1-H)
 *(sqr(max(fvc::interpolate(alpha1), scalar(0))) + H - H*fvc::interpolate(alpha1)*fvc::interpolate(alpha1))
 /((fvc::interpolate(alpha1)+(1-fvc::interpolate(alpha1))*H)*(fvc::interpolate(alpha1)+(1-fvc::interpolate(alpha1))*H))
 *fvc::snGrad(alpha1)*mesh.magSf()
);

surfaceScalarField DpsF
(
 "DpsF",
 fvc::interpolate(D)
 *(sqr(max(fvc::interpolate(alpha1), scalar(0))) + H - H*fvc::interpolate(alpha1)*fvc::interpolate(alpha1))
 /(fvc::interpolate(alpha1)+(1-fvc::interpolate(alpha1))*H)
);
/*
surfaceScalarField phinsngF
(
 "phinsngF",
 fvc::interpolate(D)*(1-fvc::interpolate(alpha1)*fvc::interpolate(alpha1))
 *fvc::snGrad(neg(alpha1-1+1e-03)*H/(alpha1+(1-alpha1)*H))*mesh.magSf()
);

surfaceScalarField DngF
(
 "DngF",
 fvc::interpolate(D)*(1-fvc::interpolate(alpha1))*H
 *neg(fvc::interpolate(alpha1)-1+1e-03)/(fvc::interpolate(alpha1)+(1-fvc::interpolate(alpha1))*H)
);
*/
/*
surfaceScalarField phiCi
(
 "phiCi",
 fvc::interpolate(D)*(pos(fvc::interpolate(alpha1)-1e-06)-H*neg(fvc::interpolate(alpha1)-1+1e-06))
 /(fvc::interpolate(alpha1)+(1-fvc::interpolate(alpha1))*H)
 *fvc::snGrad(alpha1)*mesh.magSf()
);
*/

fvScalarMatrix CEqn
(
	fvm::ddt(C) + fvm::div(phi, C) + fvm::div(phipsngF, C) - fvm::laplacian(DpsF, C)
	
);

/*
fvScalarMatrix CEqn
(
 	fvm::ddt(C) + fvm::div(phi,C) - fvm::laplacian(fvc::interpolate(D), C) + fvm::div(phiCi, C)
);

*/
CEqn.solve();



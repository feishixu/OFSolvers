
dimensionedScalar deltasng("deltasng", dimLength/dimArea, small);
dimensionedScalar dLength("dLength", dimLength, 1);
surfaceScalarField xI(mag(fvc::snGrad(alpha1))*dLength);
forAll(xI,i)
{
	if(xI[i]<500)
		xI[i]=1;
	else
		xI[i]=0;
}


surfaceScalarField phiDi
(
 "phiDi",
 fvc::interpolate(D)*xI
);
/*
surfaceScalarField phiCi
(
 "phiCi",
 ((fvc::interpolate(D)*(1-H)/(fvc::interpolate(alpha1)+(1-fvc::interpolate(alpha1))*H))*fvc::snGrad(alpha1))*mesh.magSf()
);
*/
fvScalarMatrix CEqn
(
	fvm::ddt(C) +fvm::div(phi,C)-fvm::laplacian(phiDi, C) 
);
CEqn.solve();

/*
	volVectorField dI(fvc::grad(alpha1));
	surfaceVectorField dIf(fvc::interpolate(dI));
	vectorField ndIf(dIf/(mag(dIf)+deltasng));

	scalar dx=5e-5;
	Balance(C, mesh,ndIf,alpha1,dx);
	ishBalance(C, mesh,ndIf,alpha1,dx);


*/

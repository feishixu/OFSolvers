scalar Pe(readScalar(transportProperties.lookup("Pe")));

scalar N_n=4*Pe;

D == alpha1*D1+(1-alpha1)*D2;

surfaceScalarField phipsngF
(
 "phipsngF",
 fvc::interpolate(D)*(1-H)*N_n
 /(fvc::interpolate(alpha1)+(1-fvc::interpolate(alpha1))*H)
 *fvc::snGrad(alpha1)*mesh.magSf()
);

forAll(phipsngF.boundaryField(), patchI)
{
	if (phipsngF.boundaryField()[patchI].patch().name() == "walli" || phipsngF.boundaryField()[patchI].patch().name() == "wallo")
	{
	        fvsPatchField<scalar>& pbf = phipsngF.boundaryFieldRef()[patchI];
		forAll(pbf, faceI)
		{
			pbf[faceI]=0;
		}
	}
}

surfaceScalarField DpsF
(
 "DpsF",
 fvc::interpolate(D)*
 (fvc::interpolate(alpha1)+neg(fvc::interpolate(alpha1)-1)*(N_n-1)*fvc::interpolate(alpha1)
  + (1-fvc::interpolate(alpha1)+neg(-fvc::interpolate(alpha1))*(N_n-1)*(1-fvc::interpolate(alpha1)))*H
  )
 /(fvc::interpolate(alpha1)+(1-fvc::interpolate(alpha1))*H)
);

forAll(DpsF.boundaryField(), patchI)
{
	if (DpsF.boundaryField()[patchI].patch().name() == "walli" || DpsF.boundaryField()[patchI].patch().name() == "wallo" )
	{
	        fvsPatchField<scalar>& pbf = DpsF.boundaryFieldRef()[patchI];
		forAll(pbf, faceI)
		{
			pbf[faceI]=0;
		}
	}
}

fvScalarMatrix CEqn
(

	fvm::ddt(C) + fvm::div(phi, C) + fvm::div(phipsngF, C) - fvm::laplacian(DpsF, C)
);	
CEqn.solve();



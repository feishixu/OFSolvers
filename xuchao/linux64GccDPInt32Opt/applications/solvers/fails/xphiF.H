#ifndef xphiF_H
#define xphiF_H

namespace Foam
{
	surfaceScalarField xphiF(const volScalarField& C,  const volScalarField& F, const fvMesh& mesh)
	{
		volScalarField CC(C);
		const labelUList& owner=mesh.owner();
		const labelUList& neighbour=mesh.neighbour();
		forAll(CC, i)
		{
			if(F[i]>0.005)
			{
				CC[i]=C[i]/11;
			}
		}
	
		surfaceScalarField HC(fvc::snGrad(CC)*mesh.magSf());
		forAll(CC, i)
		{
			if(F[i]>0.005)
			{
				CC[i]=10*C[i]/11;
			}
		}
		surfaceScalarField LC(fvc::snGrad(CC)*mesh.magSf());
		surfaceScalarField rC(LC);
		forAll(rC, i)
		{
			rC[i]=0;
		}

		forAll(LC, facei)
		{	
			label oc=owner[facei];
			label nc=neighbour[facei];
			if(mag(F[oc]-F[nc])>0.005)
			{
				label hc=((F[oc]-F[nc])>0)? oc: nc;
				label lc=((F[oc]-F[nc])<0)? oc: nc;
				if(mag(F[hc]-1)<1e-05)
					rC[facei]=HC[facei];
				if(mag(F[lc])<1e-05)
					rC[facei]=LC[facei];

			}
		}
		return rC;
	}
}



#endif

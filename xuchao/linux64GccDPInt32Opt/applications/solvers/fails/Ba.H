#ifndef Ba_H
#define Ba_H
namespace Foam
{
	void Balance(volScalarField& C,  const volScalarField& F, const fvMesh& mesh)
	{
		
		const surfaceScalarField faceFlux(fvc::snGrad(F));
		surfaceScalarField Cf(fvc::interpolate(C));
		surfaceScalarField Ff(fvc::interpolate(F));
		/********************************************************/

		const labelUList& owner=mesh.owner();
		const labelUList& neighbour=mesh.neighbour();
		/********************************************************/

		forAll(faceFlux, facei)
		{
		
			if(mag(F[owner[facei]]-F[neighbour[facei]])>0.995)
			{
				label hC=(F[owner[facei]]>F[neighbour[facei]])? owner[facei]:neighbour[facei];
				label lC=(F[owner[facei]]<F[neighbour[facei]])? owner[facei]:neighbour[facei];
				scalar CJ=C[hC]*Ff[facei]+(1-Ff[facei])*C[lC];
			     
				C[hC]=CJ/(Ff[facei]+(1-Ff[facei])*10);
				C[lC]=10*C[hC];

			}	
			


		}


		Info<<"Balance--------------------------"<<endl;
		

		

	}
}
			

#endif

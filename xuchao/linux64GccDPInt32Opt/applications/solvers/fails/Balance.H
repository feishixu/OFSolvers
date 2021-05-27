#ifndef Balance_H
#define Balance_H
namespace Foam
{
	void Balance(volScalarField& C,const fvMesh& mesh, const vectorField& dI, const volScalarField& F, scalar dx)
	{

		const vectorField Mf=mesh.Cf();

		vectorField Ll=dI*dx+Mf;
		vectorField Lr=Mf-dI*dx;

		const surfaceScalarField faceFlux(fvc::snGrad(F));
		surfaceScalarField Cf(fvc::interpolate(C));
		surfaceScalarField Ff(fvc::interpolate(F));


		forAll(faceFlux, facei)
		{
		
			if (mag(faceFlux[facei])>500)
			{
			
			     label  cFl=mesh.findCell(Ll[facei]);
			     label  cFr=mesh.findCell(Lr[facei]);
			     if(mag(F[cFl]-F[cFr])>0.5 && (cFl != -1) && (cFr != -1))
			     {
				     label cC=((F[cFl]-F[cFr])>0)? cFl : cFr;
				     label dC=((F[cFl]-F[cFr])>0)? cFr : cFl;
				     scalar CJ=C[cC]*Ff[facei]+(1-Ff[facei])*C[dC];
				    // Info<<"CjCjCJCJCJCJ--------------------"<<CJ
				    //     <<"********************"<<Ff[facei]<<endl;
				     C[cC]=CJ/(Ff[facei]+(1-Ff[facei])*10);
				     C[dC]=10*C[cC];
				  //  C[cC]=111;
				   // C[dC]=111;
			     }

			}

		}


		Info<<"Balance--------------------------"<<endl;
		

		

	}
}
			

#endif

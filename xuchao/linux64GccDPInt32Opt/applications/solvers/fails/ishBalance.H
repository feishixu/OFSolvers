#ifndef ishBalance_H
#define ishBalance_H
namespace Foam
{
	bool ishBalance(const volScalarField& C,const fvMesh& mesh, const vectorField& dI, const volScalarField& F, scalar dx)
	{
		bool result = false;	

		const vectorField Cf=mesh.Cf();

		vectorField Ll=dI*dx+Cf;
		vectorField Lr=Cf-dI*dx;

		const surfaceScalarField faceFlux(fvc::snGrad(F));

		scalar sum=0;
		label index=0;

		forAll(faceFlux, facei)
		{
		
			if (mag(faceFlux[facei])>1)
			{
			
			     label  cFl=mesh.findCell(Ll[facei]);
			     label  cFr=mesh.findCell(Lr[facei]);
			     if(mag(F[cFl]-F[cFr])>0.5 && (cFl != -1) && (cFr != -1))
			     {
				     label cC=((F[cFl]-F[cFr])>0)? cFl : cFr;
				     label dC=((F[cFl]-F[cFr])>0)? cFr : cFl;
				     sum += (10*C[cC]-C[dC]);
				     ++index;
			     }

			}

		}
		sum = sum / (index+small);


		if (sum < 1e-2 && index && sum >=0)
		{
			result=true;
		
		}	
		Info<<"ishBalance--------------------------"<<result
			<<"------------------  "<<sum<<endl;
		

		

		return result;
	}

	bool ishOut(const volScalarField& C,const fvMesh& mesh, const vectorField& dI, const volScalarField& F, scalar dx)
	{
		bool result = false;	
		const vectorField& Cf=mesh.C();
		vectorField Ll=dI*dx+Cf;
		vectorField Lr=Cf-dI*dx;


		const surfaceScalarField& faceFlux(fvc::snGrad(F));

		scalar sum=0;
		label index=1;

		forAll(faceFlux, facei)
		{
		
			if (mag(faceFlux[facei])>1)
			{
			
			     label  cFl=mesh.findCell(Ll[facei]);
			     label  cFr=mesh.findCell(Lr[facei]);
			     if(mag(F[cFl]-F[cFr])>0.5 && (cFl != -1) && (cFr != -1))
			     {
				     label cC=((F[cFl]-F[cFr])>0)? cFl : cFr;
				     label dC=((F[cFl]-F[cFr])>0)? cFr : cFl;
				     sum += (10*C[cC]-C[dC]);
				     ++index;
			     }

			}

		}
		if(index>1)
			--index;
		sum = sum / (index);

		if (sum < 0 && index)
		{
			result=true;
		}

		
		Info<<"ishOut--------------------------"<<result
			<<"------------------  "<<sum<<endl;
		

		return result;
	}
		
		
}

#endif

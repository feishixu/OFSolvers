#ifndef smoother_H
#define smoother_H
#include "extrapolatedCalculatedFvPatchField.H"

namespace Foam
{
	volScalarField vofSmoother(const volScalarField& alphal)
	{
		const fvMesh& mesh=alphal.mesh();
		tmp<GeometricField<scalar, fvPatchField, volMesh>> tvf
		(
			new GeometricField<scalar, fvPatchField, volMesh>
			(
				IOobject
				(
				 	"su",
					alphal.instance(),
					mesh,
					IOobject::NO_READ,
					IOobject::NO_WRITE
				),
				mesh,
				dimensioned<scalar>
				(
				 	"0",
					dimless,
					Zero
				),
				extrapolatedCalculatedFvPatchField<scalar>::typeName
			)
		);
		volScalarField& alphabt=tvf.ref();
		
		tmp<GeometricField<scalar, fvPatchField, volMesh>> tvs
		(
			new GeometricField<scalar, fvPatchField, volMesh>
			(
				IOobject
				(
				 	"ss",
					alphal.instance(),
					mesh,
					IOobject::NO_READ,
					IOobject::NO_WRITE
				),
				mesh,
				dimensioned<scalar>
				(
				 	"0",
					dimless,
					Zero
				),
				extrapolatedCalculatedFvPatchField<scalar>::typeName
			)
		);
		volScalarField& ss=tvs.ref();

		const surfaceScalarField alphaf(fvc::interpolate(alphal));
		const labelUList& owner=mesh.owner();
		const labelUList& neighbour=mesh.neighbour();

		forAll(owner, facei)
		{
			alphabt[owner[facei]] += (alphaf[facei]*mesh.magSf()[facei]);
			alphabt[neighbour[facei]] += (alphaf[facei]*mesh.magSf()[facei]);
			ss[owner[facei]] += mesh.magSf()[facei];
			ss[neighbour[facei]] += mesh.magSf()[facei];
		}

		forAll(mesh.boundary(), patchi)
		{
			const labelUList& pFaceCells = mesh.boundary()[patchi].faceCells();
			const fvsPatchField<scalar>& paf = alphaf.boundaryField()[patchi];
			const fvPatch& patch = mesh.boundary()[patchi];

			forAll(mesh.boundary()[patchi], facei)
			{
				alphabt[pFaceCells[facei]] += (paf[facei]*patch.magSf()[facei]);
				ss[pFaceCells[facei]] += patch.magSf()[facei];
			}
		}

			ss.correctBoundaryConditions();
			alphabt /= ss;
			alphabt.correctBoundaryConditions();

			Info<<"smoother.................."<<endl;

			return alphabt;
	}
}

#endif

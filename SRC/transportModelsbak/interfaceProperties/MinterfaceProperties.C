/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "MinterfaceProperties.H"
#include "MalphaContactAngleFvPatchScalarField.H"
#include "mathematicalConstants.H"
#include "surfaceInterpolate.H"
#include "fvcDiv.H"
#include "fvcGrad.H"
#include "fvcSnGrad.H"
#include "smoother.H"

// * * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * //

const Foam::scalar Foam::MinterfaceProperties::convertToRad =
    Foam::constant::mathematical::pi/180.0;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Correction for the boundary condition on the unit normal nHat on
// walls to produce the correct contact angle.

// The dynamic contact angle is calculated from the component of the
// velocity on the direction of the interface, parallel to the wall.

void Foam::MinterfaceProperties::correctContactAngle
(
    surfaceVectorField::Boundary& nHatb,
    const surfaceVectorField::Boundary& gradAlphaf,
    const volVectorField::Boundary& gradAlphab,
    const surfaceVectorField& nHatfv
) const
{
    const fvMesh& mesh = alpha1_.mesh();
    const volScalarField::Boundary& abf = alpha1_.boundaryField();
    const surfaceVectorField& Sf=mesh.Sf();
    const scalarField& Vc=mesh.V();
    

    const fvBoundaryMesh& boundary = mesh.boundary();

    vectorField nHatfvvf(nHatfv);
    vectorField Sfvf(Sf);



    forAll(boundary, patchi)
    {
	const vectorField& pnHat=nHatfv.boundaryField()[patchi];
	const vectorField& pSf=Sf.boundaryField()[patchi];
	nHatfvvf.append(pnHat);
	Sfvf.append(pSf);
    }


    forAll(boundary, patchi)
    {
        if (isA<MalphaContactAngleFvPatchScalarField>(abf[patchi]))
        {
            MalphaContactAngleFvPatchScalarField& acap =
                const_cast<MalphaContactAngleFvPatchScalarField&>
                (
                    refCast<const MalphaContactAngleFvPatchScalarField>
                    (
                        abf[patchi]
                    )
                );
	/***********************************************************************/
    		const  labelUList& faceCells=acap.patch().faceCells();
   		const  label startFace=acap.patch().start();
    		scalarField K (acap.size(), 0);
    		const labelUList& owner=mesh.owner();
    		const labelUList& neighbour=mesh.neighbour();
    		forAll(faceCells, pFacei)
    		{
	    		label celli=faceCells[pFacei];//cell global index
	    		const cell& cFaces=mesh.cells()[celli];//cell
	    		forAll(cFaces, i)// for every face
	    		{
		    		label facei=cFaces[i];//one global face index
		    		if(facei < mesh.nInternalFaces()) // for interanl face
		    		{
			    		if(owner[facei]==celli)//owner face
			    		{
				    		K[pFacei] += nHatfv[facei] & Sf[facei]; 
			    		}
			    		else if(neighbour[facei]==celli)//neighbour face
			    		{
				    		K[pFacei] -= nHatfv[facei] & Sf[facei]; 
			    		}
		    		}
		   		 else if(facei < startFace || (facei< Sfvf.size() && facei >= (startFace + acap.size())))//empty faces?
		    		{
					Info<<"facei "<<facei<<endl;
				    K[pFacei] += nHatfvvf[facei] & Sfvf[facei]; //other patch 

		    		}
		    

	    		}
	    
    		}

   		 forAll(K, i)
    		{
	    		K[i] /= Vc[i];
    		}
    







	/***********************************************************************/


            fvsPatchVectorField& nHatp = nHatb[patchi];
            const scalarField theta
            (
                convertToRad*acap.theta(U_.boundaryField()[patchi], nHatp, gradAlphab[patchi], K)
                //convertToRad*acap.theta(U_.boundaryField()[patchi], nHatp)
            );

            const vectorField nf
            (
                boundary[patchi].nf()
            );

            // Reset nHatp to correspond to the contact angle

            const scalarField a12(nHatp & nf);
            const scalarField b1(cos(theta));

            scalarField b2(nHatp.size());
            forAll(b2, facei)
            {
                b2[facei] = cos(acos(a12[facei]) - theta[facei]);
            }

            const scalarField det(1.0 - a12*a12);

            scalarField a((b1 - a12*b2)/det);
            scalarField b((b2 - a12*b1)/det);

            nHatp = a*nf + b*nHatp;
            nHatp /= (mag(nHatp) + deltaN_.value());

            acap.gradient() = (nf & nHatp)*mag(gradAlphaf[patchi]);
            acap.evaluate();
        }
    }
}


void Foam::MinterfaceProperties::calculateK()
{
    const fvMesh& mesh = alpha1_.mesh();
    const surfaceVectorField& Sf = mesh.Sf();
/********************************************************************************/

    volScalarField alpha1b(vofSmoother(alpha1_));
    int m=2;
    for(int i=0; i<m-1; i++)
    {
	    alpha1b=vofSmoother(alpha1b);
    }

/*******************************************************************************/
    // Cell gradient of alpha
    const volVectorField gradAlpha(fvc::grad(alpha1b, "nHat"));

    // Interpolated face-gradient of alpha
    surfaceVectorField gradAlphaf(fvc::interpolate(gradAlpha));

    // gradAlphaf -=
    //    (mesh.Sf()/mesh.magSf())
    //   *(fvc::snGrad(alpha1_) - (mesh.Sf() & gradAlphaf)/mesh.magSf());

    // Face unit interface normal
    surfaceVectorField nHatfv(gradAlphaf/(mag(gradAlphaf) + deltaN_));
    // surfaceVectorField nHatfv
    // (
    //     (gradAlphaf + deltaN_*vector(0, 0, 1)
    //    *sign(gradAlphaf.component(vector::Z)))/(mag(gradAlphaf) + deltaN_)
    // );
    correctContactAngle(nHatfv.boundaryFieldRef(), gradAlphaf.boundaryField(), gradAlpha.boundaryField(), nHatfv);

    // Face unit interface normal flux
    nHatf_ = nHatfv & Sf;

    // Simple expression for curvature
    K_ = -fvc::div(nHatf_);

    // Complex expression for curvature.
    // Correction is formally zero but numerically non-zero.
    /*
    volVectorField nHat(gradAlpha/(mag(gradAlpha) + deltaN_));
    forAll(nHat.boundaryField(), patchi)
    {
        nHat.boundaryField()[patchi] = nHatfv.boundaryField()[patchi];
    }

    K_ = -fvc::div(nHatf_) + (nHat & fvc::grad(nHatfv) & nHat);
    */
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::MinterfaceProperties::MinterfaceProperties
(
    const volScalarField& alpha1,
    const volVectorField& U,
    const IOdictionary& dict
)
:
    transportPropertiesDict_(dict),
    cAlpha_
    (
        readScalar
        (
            alpha1.mesh().solverDict(alpha1.name()).lookup("cAlpha")
        )
    ),

    sigmaPtr_(surfaceTensionModel::New(dict, alpha1.mesh())),

    deltaN_
    (
        "deltaN",
        1e-8/pow(average(alpha1.mesh().V()), 1.0/3.0)
    ),

    alpha1_(alpha1),
    U_(U),

    nHatf_
    (
        IOobject
        (
            "nHatf",
            alpha1_.time().timeName(),
            alpha1_.mesh()
        ),
        alpha1_.mesh(),
        dimensionedScalar(dimArea, 0)
    ),

    K_
    (
        IOobject
        (
            "MinterfaceProperties:K",
            alpha1_.time().timeName(),
            alpha1_.mesh()
        ),
        alpha1_.mesh(),
        dimensionedScalar(dimless/dimLength, 0)
    )
{
    calculateK();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::MinterfaceProperties::sigmaK() const
{
    return sigmaPtr_->sigma()*K_;
}


Foam::tmp<Foam::surfaceScalarField>
Foam::MinterfaceProperties::surfaceTensionForce() const
{
    return fvc::interpolate(sigmaK())*fvc::snGrad(alpha1_);
}


Foam::tmp<Foam::volScalarField>
Foam::MinterfaceProperties::nearInterface() const
{
    return pos0(alpha1_ - 0.01)*pos0(0.99 - alpha1_);
}


void Foam::MinterfaceProperties::correct()
{
    calculateK();
}


bool Foam::MinterfaceProperties::read()
{
    alpha1_.mesh().solverDict(alpha1_.name()).lookup("cAlpha") >> cAlpha_;
    sigmaPtr_->readDict(transportPropertiesDict_);

    return true;
}


// ************************************************************************* //

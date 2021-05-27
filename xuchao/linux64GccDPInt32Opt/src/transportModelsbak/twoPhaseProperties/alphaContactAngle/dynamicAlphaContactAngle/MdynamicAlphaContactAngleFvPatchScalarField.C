/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

#include "MdynamicAlphaContactAngleFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volMesh.H"
#include "mathematicalConstants.H"
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::MdynamicAlphaContactAngleFvPatchScalarField::
MdynamicAlphaContactAngleFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    MalphaContactAngleFvPatchScalarField(p, iF),
    theta0_(0.0),
    uTheta_(0.0),
    thetaA_(0.0),
    thetaR_(0.0)
{}


Foam::MdynamicAlphaContactAngleFvPatchScalarField::
MdynamicAlphaContactAngleFvPatchScalarField
(
    const MdynamicAlphaContactAngleFvPatchScalarField& gcpsf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    MalphaContactAngleFvPatchScalarField(gcpsf, p, iF, mapper),
    theta0_(gcpsf.theta0_),
    uTheta_(gcpsf.uTheta_),
    thetaA_(gcpsf.thetaA_),
    thetaR_(gcpsf.thetaR_)
{}


Foam::MdynamicAlphaContactAngleFvPatchScalarField::
MdynamicAlphaContactAngleFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    MalphaContactAngleFvPatchScalarField(p, iF, dict),
    theta0_(readScalar(dict.lookup("theta0"))),
    uTheta_(readScalar(dict.lookup("uTheta"))),
    thetaA_
    (
        dict.found("thetaA")
      ? readScalar(dict.lookup("thetaA"))
      : readScalar(dict.lookup("thetaRec"))
    ),
    thetaR_
    (
        dict.found("thetaR")
      ? readScalar(dict.lookup("thetaR"))
      : readScalar(dict.lookup("thetaAdv"))
    )
{
    evaluate();
}


Foam::MdynamicAlphaContactAngleFvPatchScalarField::
MdynamicAlphaContactAngleFvPatchScalarField
(
    const MdynamicAlphaContactAngleFvPatchScalarField& gcpsf
)
:
    MalphaContactAngleFvPatchScalarField(gcpsf),
    theta0_(gcpsf.theta0_),
    uTheta_(gcpsf.uTheta_),
    thetaA_(gcpsf.thetaA_),
    thetaR_(gcpsf.thetaR_)
{}


Foam::MdynamicAlphaContactAngleFvPatchScalarField::
MdynamicAlphaContactAngleFvPatchScalarField
(
    const MdynamicAlphaContactAngleFvPatchScalarField& gcpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    MalphaContactAngleFvPatchScalarField(gcpsf, iF),
    theta0_(gcpsf.theta0_),
    uTheta_(gcpsf.uTheta_),
    thetaA_(gcpsf.thetaA_),
    thetaR_(gcpsf.thetaR_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::MdynamicAlphaContactAngleFvPatchScalarField::theta
(
    const fvPatchVectorField& Up,
    const fvsPatchVectorField& nHat,
    const fvPatchVectorField& gradAlphap,
    const scalarField& K
) const
{

	Info<<"theta---------------"<<endl;
    const  fvMesh& mesh=(this->internalField()).mesh();
    const scalar toRad=Foam::constant::mathematical::pi/180.0;


    IOdictionary transportProperties
    (
     	IOobject
    	(
 	"transportProperties",
     	"constant",
     	mesh,
     	IOobject::MUST_READ_IF_MODIFIED,
     	IOobject::NO_WRITE
    	)
    );
	Istream& is=transportProperties.lookup("phases");
	token ptoken(is);
	word name1(is);
	word name2(is);

	dictionary phase1 (transportProperties.subDict(name1));
	dictionary phase2 (transportProperties.subDict(name2));

	scalar xrho1 (readScalar(phase1.lookup("rho")));
	scalar xrho2 (readScalar(phase2.lookup("rho")));
	scalar xmu1 (readScalar(phase1.lookup("nu")));
	scalar xmu2 (readScalar(phase2.lookup("nu")));
	
	xmu1 *= xrho1;
	xmu2 *= xrho2;

	scalar sigma (readScalar(transportProperties.lookup("sigma")));

	const scalarField alphac(this->patchInternalField());

	scalarField xrho (xrho1*alphac+xrho2*(1-alphac)); 
	scalarField xmu (xmu1*alphac+xmu2*(1- alphac)); 

    IOdictionary controlDict
    (
     	IOobject
    	(
 	"controlDict",
     	"system",
     	mesh,
     	IOobject::MUST_READ_IF_MODIFIED,
     	IOobject::NO_WRITE
    	)
    );

	scalar dt (readScalar(controlDict.lookup("deltaT")));

    vectorField Uc(Up.patchInternalField());
    vectorField Un(Uc/(mag(Uc)+small));
    vectorField gradAlphaC(gradAlphap.patchInternalField());

    vectorField Force(xrho*Uc/dt);

    scalarField Kfull(mag(Force)/((gradAlphaC & Un)*sigma +small));
    

    const  labelUList& faceCells=this->patch().faceCells();
    const  label startFace=this->patch().start();

    scalarField pSf(this->patch().magSf());
    const scalarField& Vc=mesh.V();
    scalarField pV(this->patch().size(), 0);
    forAll(pV, pFacei)
    {
	    pV[pFacei]=Vc[faceCells[pFacei]];
    }

    scalarField Mtheta(this->patch().size(),0);

    
    forAll(Mtheta, pFacei)
    {
	    scalar vtheta=(Kfull[pFacei]-K[pFacei])*pV[pFacei]/(pSf[pFacei]+small);
   
	    if(vtheta>=-1 && vtheta<=1)
	    {
	    	Mtheta[pFacei]=acos(vtheta)/toRad;
	    }
	    else if( vtheta < -1)
	    {
		    Mtheta[pFacei]=180;
	    }
    
    }
    


    

/************************************************************/		
/*
    const scalarField a12(nHat & nf);
    
    vectorField Upwall(Up-(nf & Up)*nf);
    scalarField AorR(Upwall & nHat);

    forAll(AorR, facei)
    {
	  if(AorR[facei]<0) //Adv  
	  {
		  if(acos(a12[facei])<=thetaA_) //static angle blance

		  {
			  uTheata_=0;

*/

/*******************************************/
    const vectorField nf(patch().nf());

    if (uTheta_ < small)
    {
        return tmp<scalarField>(new scalarField(size(), theta0_));
    }


    vectorField Uwall(Up.patchInternalField() - Up);
    Uwall -= (nf & Uwall)*nf;

    // Find the direction of the interface parallel to the wall
    vectorField nWall(nHat - (nf & nHat)*nf);

    // Normalise nWall
    nWall /= (mag(nWall) + small);

    // Calculate Uwall resolved normal to the interface parallel to
    // the interface
    scalarField uwall(nWall & Uwall);


    /*****************************************************/
    vectorField Ucl=


    /*****************************************************/

    return theta0_ + (thetaA_ - thetaR_)*tanh(uwall/uTheta_);
//      return tmp<scalarField>(new scalarField(size(),theta0_));
}


void Foam::MdynamicAlphaContactAngleFvPatchScalarField::write(Ostream& os) const
{
    MalphaContactAngleFvPatchScalarField::write(os);
    writeEntry(os, "theta0", theta0_);
    writeEntry(os, "uTheta", uTheta_);
    writeEntry(os, "thetaA", thetaA_);
    writeEntry(os, "thetaR", thetaR_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        MdynamicAlphaContactAngleFvPatchScalarField
    );
}


// ************************************************************************* //

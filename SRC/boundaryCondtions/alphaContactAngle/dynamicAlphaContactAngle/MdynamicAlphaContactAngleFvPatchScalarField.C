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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::MdynamicAlphaContactAngleFvPatchScalarField::
MdynamicAlphaContactAngleFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    alphaContactAngleFvPatchScalarField(p, iF),
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
    alphaContactAngleFvPatchScalarField(gcpsf, p, iF, mapper),
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
    alphaContactAngleFvPatchScalarField(p, iF, dict),
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
    alphaContactAngleFvPatchScalarField(gcpsf),
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
    alphaContactAngleFvPatchScalarField(gcpsf, iF),
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
    const fvsPatchVectorField& nHat
) const
{
    const  fvMesh& mesh=(this->internalField()).mesh();

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

	Info<<"name1 "<<name1<<endl;
	Info<<"name2 "<<name2<<endl;
	dictionary phase1 (transportProperties.subDict(name1));
	dictionary phase2 (transportProperties.subDict(name2));

	scalar xrho1 (readScalar(phase1.lookup("rho")));
	scalar xrho2 (readScalar(phase2.lookup("rho")));

	Info<<"xrho1 "<<xrho1<<endl;
	Info<<"xrho2 "<<xrho2<<endl;

	scalarField xrho (xrho1*(*this)+xrho2*(1-(*this))); 
	Info<<"xrho "<<xrho<<endl;

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
	Info<<"dT "<<dt<<endl;

    vectorField Uc(Up.patchInternalField());
    Info<<"Uc "<<Uc<<endl;

    vectorField Force(xrho*Uc/dt);
    Info<<"Force "<<Force<<endl;
/*		
    const vectorField nf(patch().nf());
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
/*
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
*/
//    return theta0_ + (thetaA_ - thetaR_)*tanh(uwall/uTheta_);
      return tmp<scalarField>(new scalarField(size(),theta0_));
}


void Foam::MdynamicAlphaContactAngleFvPatchScalarField::write(Ostream& os) const
{
    alphaContactAngleFvPatchScalarField::write(os);
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

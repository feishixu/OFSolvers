/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021 OpenFOAM Foundation
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

#include "parabolicVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

//Foam::scalar Foam::parabolicVelocityFvPatchVectorField::t() const
//{
//    return db().time().timeOutputValue();
//}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// 1
Foam::parabolicVelocityFvPatchVectorField::
parabolicVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
	maxvalue_(0),
	n_(1,0,0),
	y_(0,1,0)
{
}

// 2
Foam::parabolicVelocityFvPatchVectorField::
parabolicVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    maxvalue_(readScalar(dict.lookup("maxvalue"))),
	n_(dict.lookup("n")),
	y_(dict.lookup("y"))
{

	Info << "Using the parabolicVelocity boundary condition" << endl;
	if (mag(n_) < SMALL || mag(y_) < SMALL)
	{
		FatalErrorIn("parabolicVelocityFvPatchVectorField(dict)")
		<< "n or y given with zero size not correct"
		<< abort(FatalError);
	}
	n_ /= mag(n_);
	y_ /= mag(y_);

    fixedValueFvPatchVectorField::evaluate();

}

// 3
Foam::parabolicVelocityFvPatchVectorField::
parabolicVelocityFvPatchVectorField
(
    const parabolicVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    maxvalue_(ptf.maxvalue_),
	n_(ptf.n_),
	y_(ptf.y_)
{}

// 4
Foam::parabolicVelocityFvPatchVectorField::
parabolicVelocityFvPatchVectorField
(
    const parabolicVelocityFvPatchVectorField& ptf
)
:
    fixedValueFvPatchVectorField(ptf),
    maxvalue_(ptf.maxvalue_),
	n_(ptf.n_),
	y_(ptf.y_)
{}

// 5
Foam::parabolicVelocityFvPatchVectorField::
parabolicVelocityFvPatchVectorField
(
    const parabolicVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(ptf, iF),
    maxvalue_(ptf.maxvalue_),
	n_(ptf.n_),
	y_(ptf.y_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//void Foam::parabolicVelocityFvPatchVectorField::autoMap
//(
//    const fvPatchFieldMapper& m
//)
//{
//    fixedValueFvPatchVectorField::autoMap(m);
//    m(fieldData_, fieldData_);
//}


//void Foam::parabolicVelocityFvPatchVectorField::rmap
//(
//    const fvPatchVectorField& ptf,
//    const labelList& addr
//)
//{
//    fixedValueFvPatchVectorField::rmap(ptf, addr);

//    const parabolicVelocityFvPatchVectorField& tiptf =
//        refCast<const parabolicVelocityFvPatchVectorField>(ptf);

//    fieldData_.rmap(tiptf.fieldData_, addr);
//}


void Foam::parabolicVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

	boundBox bb(patch().patch().localPoints(),true);
	vector ctr = 0.5*(bb.max()+bb.min());
	const vectorField &c = patch().Cf();

	scalarField coord = 2*((c-ctr)&y_)/((bb.max()-bb.min())&y_);
	vectorField :: operator =(n_*maxvalue_ * (1.0-sqr(coord)));

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::parabolicVelocityFvPatchVectorField::write
(
    Ostream& os
) const
{
    fvPatchVectorField::write(os);
    writeEntry(os, "n", n_);
    writeEntry(os, "y", y_);
    writeEntry(os, "maxvalue", maxvalue_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * Build Macro Function  * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        parabolicVelocityFvPatchVectorField
    );
}

// ************************************************************************* //

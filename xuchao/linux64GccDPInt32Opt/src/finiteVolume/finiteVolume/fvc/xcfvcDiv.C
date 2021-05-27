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

#include "xfvcDiv.H"
#include "fvMesh.H"
#include "fvcSurfaceIntegrate.H"
#include "divScheme.H"
#include "convectionScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fvc
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
xdiv
(
    const GeometricField<Type, fvPatchField, volMesh>& F,
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const word& name
)
{
    
    return fv::convectionScheme<Type>::New
    (
        vf.mesh(),
        ,
        vf.mesh().divScheme(name)
    ).ref().fvcDiv(flux, vf);
}




template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
xdiv
(
    const GeometricField<Type, fvPatchField, volMesh>& F,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return fvc::xdiv
    (
        F, vf, "xdiv("+F.name()+','+vf.name()+')'
    );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvc

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

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

Class
    Foam::fv::xcgaussConvectionScheme

Description
    Basic second-order convection using face-gradients and Gauss' theorem.

SourceFiles
    xcgaussConvectionScheme.C

\*---------------------------------------------------------------------------*/

#ifndef xcgaussConvectionScheme_H
#define xcgaussConvectionScheme_H

#include "convectionScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                       Class xcgaussConvectionScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class xcgaussConvectionScheme
:
    public fv::convectionScheme<Type>
{
    // Private Data

        tmp<surfaceInterpolationScheme<Type>> tinterpScheme_;


public:

    //- Runtime type information
    TypeName("xcGauss");


    // Constructors

        //- Construct from flux and interpolation scheme
        xcgaussConvectionScheme
        (
            const fvMesh& mesh,
            const volScalarField& F,
            const tmp<surfaceInterpolationScheme<Type>>& scheme
        )
        :
            convectionScheme<Type>(mesh, fvc::interpolate(F)),
            tinterpScheme_(scheme)
        {}

        //- Construct from flux and Istream
        xcgaussConvectionScheme
        (
            const fvMesh& mesh,
            const volScalarField& F,
            Istream& is
        )
        :
            convectionScheme<Type>(mesh,  fvc::interpolate(F)),
            tinterpScheme_
            (
                surfaceInterpolationScheme<Type>::New(mesh, (fvc::interpolate(F)-0.5)*fvc::snGrad(F), is)
            )
        {}

        //- Disallow default bitwise copy construction
        xcgaussConvectionScheme(const xcgaussConvectionScheme&) = delete;


    // Member Functions

        const surfaceInterpolationScheme<Type>& interpScheme() const;

        tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> interpolate
        (
            const surfaceScalarField&,
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;

        tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> flux
        (
            const surfaceScalarField&,
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;

        tmp<fvMatrix<Type>> fvmDiv
        (
            const surfaceScalarField&,
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;

        tmp<GeometricField<Type, fvPatchField, volMesh>> fvcDiv
        (
            const surfaceScalarField&,
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const xcgaussConvectionScheme&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "xcgaussConvectionScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

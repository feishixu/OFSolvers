namespace Foam
{

namespace fvm
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<fvMatrix<Type>>
xdiv
(
    const surfaceScalarField& flux,
    const surfaceScalarField& F,
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const word& name
)
{
    return fv::convectionScheme<Type>::New
    (
        vf.mesh(),
        F,
        vf.mesh().divScheme(name)
    )().fvmDiv(flux, vf);
}



template<class Type>
tmp<fvMatrix<Type>>
xdiv
(
    const surfaceScalarField& flux,
    const surfaceScalarField& F,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return fvm::xdiv(flux, F, vf, "xdiv("+flux.name()+','+F.name()+ ',' +vf.name()+')');
}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvm

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

The smoother solver is based on code:
git@github.com:feishixu/OF-kva_interfaceProperties

The solver kva_interFoam is compiled for `of5` in $MY_FOAM_APPBIN
To use the solver:

1. add the following to your constant/transportProperties dictionary:

```C++
curvatureModel      vofsmooth; // normal;
vofsmoothCoeffs
{
    numSmoothingIterations 2; // default: 2
}
```

2. then run the solver
```C++
kva_interFoam
```

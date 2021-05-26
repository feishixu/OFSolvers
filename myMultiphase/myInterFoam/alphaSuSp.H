zeroField Su;
//zeroField Sp;
zeroField divU;

volScalarField Sp ("Sp", (1.0/rho1 * cflux)/(alpha1 + SMALL));
//volScalarField Su ("Su", 0*divU);
/*

            MULES::explicitSolve
            (
                geometricOneField(),
                alpha1,
                phi,
		Sp,
		Su,
                talphaPhiCorr.ref(),
                vDotvmcAlphal,
                (divU*alpha1 + vDotcAlphal)(),
                oneField(),
                zeroField()
            );
*/

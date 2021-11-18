zeroField Su;
//zeroField Sp;
zeroField divU;

// volScalarField Sp ("Sp", (1.0/rho1 * k_pc * C)/(alpha1 + SMALL));

volScalarField C1 ("C1", alpha1*max(C_0,C-C_thr)+(1-alpha1)*min(C_0,C-C_thr));

volScalarField Sp ("Sp", (1.0/rho1 * k_pc*C1)/(alpha1 + SMALL));


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

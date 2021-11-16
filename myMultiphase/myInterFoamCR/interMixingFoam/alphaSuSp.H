zeroField Su;
zeroField Sp;
zeroField divU;

//** Creating Source Terms

volScalarField Sp1 ("Sp1", (1.0 / rho1 *k_pc*alpha2*alpha3)/(alpha1 + SMALL));
volScalarField Sp2 ("Sp2", (1.0 / rho2 *k_pc*alpha2*alpha3)/(alpha2 + SMALL));

/* volScalarField Su1 ("Su1", 1.0 / rho1 *kC*C);
volScalarField Su2 ("Su2", 1.0 / rho2 *kC*C);

    //Evaluating Source Terms
    forAll (alpha1, celli)
    {
    if(C[celli] >= 1e-4)
    {
         //Su1[celli] = 1.0 / rho1.value()[celli] *0.1;
         Su1[celli] = 1.0 / rho1.value()*kC.value()*C[celli];
    }
    else
    {
         Su1[celli] = 0.0;
    }
    }

    forAll (alpha2, celli)
    {
    if(C[celli] >= 1e-4)
    {
         Su2[celli] = 1.0 / rho2.value()*kC.value()*C[celli];
         //Su2[celli] = 1.0;
    }
    else
    {
         Su2[celli] = 0.0;
    }
    } */





// =======================================================

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

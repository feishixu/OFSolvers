// - T field

volScalarField C1 ("C1", alpha1*max(C_0,C-C_thr)+(1-alpha1)*min(C_0,C-C_thr));

//volScalarField Sp ("Sp", (1.0/rho1 * k_pc*C1)/(alpha1 + SMALL));

{
	fvScalarMatrix TEqn
(
			  fvm::ddt(T)
			+ fvm::div(phi, T)
			- fvm::laplacian(DT, T)
			==
			  k_T*alpha1*(1-alpha1)
);

	TEqn.solve();
}





{
    fvScalarMatrix TEqn
    (
        fvm::ddt(rho, T) + fvm::div(rhoPhi, T) - fvm::Sp(contErr, T)
      - fvm::laplacian(turbulence.alphaEff(), T)
      + (
             divUp()// - contErr/rho*p
          +  (fvc::ddt(rho, K) + fvc::div(rhoPhi, K))() - contErr*K
        )
       *(
           alpha1()/mixture.thermo1().Cv()()
         + alpha2()/mixture.thermo2().Cv()()
        )
     ==
        fvOptions(rho, T)
    );

    TEqn.relax();

    fvOptions.constrain(TEqn);

    TEqn.solve();

    fvOptions.correct(T);

    Info<< "Min/max T:" << min(T).value() << ' '
            << max(T).value() << endl;
    T.max(TMin);
    T.min(TMax);

    mixture.correctThermo();
    mixture.correct();
}

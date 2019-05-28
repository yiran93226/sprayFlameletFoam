{
    fvScalarMatrix dhEqn
    (
        fvm::ddt(rho, dh) + mvConvection->fvmDiv(phi, dh)
      - fvm::laplacian(turbulence->alphaEff(), dh)
     ==
        fvOptions(rho, dh)
      // + parcels.hsTrans() / (mesh.V()*runTime.deltaT())
      + parcels.Sh(dh)
    ); 

    dhEqn.relax();

    fvOptions.constrain(dhEqn);

    dhEqn.solve(mesh.solver("h"));

    fvOptions.correct(dh);
}
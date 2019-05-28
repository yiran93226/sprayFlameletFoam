{
    fvScalarMatrix dhEqn
    (
        fvm::ddt(rho, dh) + mvConvection->fvmDiv(phi, dh)
      - fvm::laplacian(turbulence->muEff(), dh)
     ==
        fvOptions(rho, dh)
      + parcels.hsTrans() / (mesh.V()*runTime.deltaT())
    ); 

    dhEqn.relax();

    fvOptions.constrain(dhEqn);

    dhEqn.solve(mesh.solver("dh"));

    fvOptions.correct(dh);
}
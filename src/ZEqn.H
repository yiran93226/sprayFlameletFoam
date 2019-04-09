tmp<fv::convectionScheme<scalar>> mvConvection
(
    fv::convectionScheme<scalar>::New
    (
       mesh,
       fields,
       phi,
       mesh.divScheme("div(phi,Z)")
    )
);

{
    fvScalarMatrix ZEqn
    (
        fvm::ddt(rho, Z) + mvConvection->fvmDiv(phi, Z)
      - fvm::laplacian(turbulence->muEff(), Z)
     ==
        fvOptions(rho, Z)
      + parcels.Srho()
    );

    ZEqn.relax();
    fvOptions.constrain(ZEqn);
    ZEqn.solve(mesh.solver("Z"));
    fvOptions.correct(Z);
    Z.max(0.0);


    if (modelType == "RAS")
    {
        fvScalarMatrix varZEqn
        (
            fvm::ddt(rho, varZ) + fvm::div(phi, varZ)
          - fvm::laplacian(turbulence->muEff(), varZ)
          - 2.0*turbulence->muEff()*magSqr(fvc::grad(Z))
          + 2.0*rho*turbulence->epsilon()/turbulence->k()*varZ
        );
        varZEqn.relax();
        varZEqn.solve(mesh.solver("varZ"));
        varZ.correctBoundaryConditions();
        varZ.max(0.0);
        chi = 2.0*turbulence->epsilon()/turbulence->k()*varZ;
    }
    else
        chi = 2*turbulence->muEff()/rho*mag(fvc::grad(Z))*mag(fvc::grad(Z));

    Info<< "Z min/max : " << min(Z).value() << ", "
        << max(Z).value() << endl;
    Info<< "varZ min/max : " << min(varZ).value() << ", "
        << max(varZ).value() << endl;
    Info<< "chi min/max : " << min(chi).value() << ", "
        << max(chi).value() << endl;
}

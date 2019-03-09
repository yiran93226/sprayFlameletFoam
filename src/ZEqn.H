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

    Info<< "Z min/max : " << min(Z).value() << ", "
        << max(Z).value() << endl;

    chi = 2*turbulence->muEff()/rho*mag(fvc::grad(Z))*mag(fvc::grad(Z));
    Info<< "chi min/max : " << min(chi).value() << ", "
        << max(chi).value() << endl;
}

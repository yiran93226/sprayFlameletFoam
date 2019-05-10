{
    volScalarField& he = thermo.he();

    const scalarField& ZCells = Z.internalField();
    const scalarField& varZCells = varZ.internalField();
    const scalarField& chiCells = chi.internalField();
    const scalarField& pCells = p.internalField();

    scalarField& heCells = he.primitiveFieldRef();
    scalarField& TCells = T.primitiveFieldRef();

    scalar Zeta = 0.0;
    scalar tT = 0.0;
    scalar tY = 0.0;
    double y = 0.0;

    volScalarField Yt(0.0*Y[0]);
    forAll(ZCells, cellI)
    {
        Zeta = Foam::sqrt(varZCells[cellI]/Foam::max(ZCells[cellI]*(1 - ZCells[cellI]), SMALL));
        flameletTable.find(ZCells[cellI], Zeta, chiCells[cellI]);
        heCells[cellI] = 0.0;
        tT = flameletTable.lookupT();
        
        // T is a function of y coordinate
        // y = mesh.C()[cellI][1];
        // tT = 300*(Foam::exp(-(1000*(y-0.01))*(1000*(y-0.01))/(2*6))/(Foam::sqrt(0.005 * Foam::constant::mathematical::pi)*2) + 1);
        TCells[cellI] = tT;

        forAll(Y, i)
        {
            scalarField& YCells = Y[i].primitiveFieldRef();
            // YCells[cellI] = flameletTable.lookupY(i);
            // tY = flameletTable.lookupY(i);
            heCells[cellI] += YCells[cellI]*thermo.composition().Hs(i, 1e5, tT);
            // heCells[cellI] += tY*thermo.composition().Hs(i, pCells[cellI], tT);
            // heCells[cellI] += YCells[cellI]*thermo.composition().Hs(i, 1e5, tT);
            Yt[cellI] += Y[i][cellI];
        }
    }

    forAll(Z.boundaryField(), patchI)
    {
        const fvPatchScalarField& pchi = chi.boundaryField()[patchI];
        const fvPatchScalarField& pvarZ = varZ.boundaryField()[patchI];
        const fvPatchScalarField& pZ = Z.boundaryField()[patchI];
        const fvPatchScalarField& pp = p.boundaryField()[patchI];

        fvPatchScalarField& phe = he.boundaryFieldRef()[patchI];
        fvPatchScalarField& pT = T.boundaryFieldRef()[patchI];
        fvPatchScalarField& pYt = Yt.boundaryFieldRef()[patchI];

        forAll(pZ, faceI)
        {
            Zeta = Foam::sqrt(pvarZ[faceI]/Foam::max(pZ[faceI]*(1 - pZ[faceI]), SMALL));
            flameletTable.find(pZ[faceI], Zeta, pchi[faceI]);
            phe[faceI] = 0.0;
            tT = flameletTable.lookupT();

            // T(y)
            // y = pZ.patch().Cf()[faceI][1];
            // tT = 300*(Foam::exp(-(1000*(y-0.01))*(1000*(y-0.01))/(2*6))/(Foam::sqrt(0.005 * Foam::constant::mathematical::pi)*2) + 1);
            pT[faceI] = tT;
            
            forAll(Y, i)
            {
                fvPatchScalarField& pY = Y[i].boundaryFieldRef()[patchI];
                // pY[faceI] = flameletTable.lookupY(i);
                // tY = flameletTable.lookupY(i);
                // phe[faceI] += pY[faceI]*thermo.composition().Hs(i, pp[faceI], tT);
                // phe[faceI] += tY*thermo.composition().Hs(i, 1e5, tT);
                phe[faceI] += pY[faceI]*thermo.composition().Hs(i, 1e5, tT);
                pYt[faceI] += pY[faceI];
            }
        }
    }
    // fvScalarMatrix EEqn
    // (
    //     fvm::ddt(rho, he) + mvConvection->fvmDiv(phi, he)
    //   + fvc::ddt(rho, K) + fvc::div(phi, K)
    //   + (
    //         he.name() == "e"
    //       ? fvc::div
    //         (
    //             fvc::absolute(phi/fvc::interpolate(rho), U),
    //             p,
    //             "div(phiv,p)"
    //         )
    //       : -dpdt
    //     )
    //   - fvm::laplacian(turbulence->alphaEff(), he)
    //  ==
    //     rho*(U&g)
    //   + parcels.Sh(he)
    //   + radiation->Sh(thermo, he)
    //   + Qdot
    //   + fvOptions(rho, he)
    // );

    // EEqn.relax();

    // fvOptions.constrain(EEqn);

    // EEqn.solve();

    // fvOptions.correct(he);

    thermo.correct();
    // radiation->correct();

    Info<< "T gas min/max   " << min(T).value() << ", "
        << max(T).value() << endl;

    Info<< "p min/max : " << min(p).value() << ", "
        << max(p).value() << endl;
    Info<< "Y oxidizer min/max : " << min(Y[oxidizerIndex]).value() << ", "
        << max(Y[oxidizerIndex]).value() << endl;
    Info<< "Yt min/max : " << min(Yt).value() << ", "
        << max(Yt).value() << endl;


    flameletTable.find(0.271, 0.0, 0.0);
    Info<< "Find Z0.271 TEMPERATURE: " << flameletTable.lookupT() << endl;
    flameletTable.find(0.01, 0.0, 0.0);
    Info<< "Find Z0.01 YO2: " << flameletTable.lookupY(oxidizerIndex) << endl;
}

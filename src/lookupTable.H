{
    volScalarField& he = thermo.he();

    const scalarField& ZCells = Z.internalField();
    const scalarField& varZCells = varZ.internalField();
    const scalarField& chiCells = chi.internalField();
    const scalarField& pCells = p.internalField();

    scalarField& heCells = he.primitiveFieldRef();

    scalar Zeta;
    scalar tT;

    forAll(Y, i)
    {
        scalarField& YCells = Y[i].primitiveFieldRef();

        forAll(ZCells, cellI)
        {
            if (i == 0)
            {
                Zeta = Foam::sqrt(varZCells[cellI]/Foam::max(ZCells[cellI]*(1 - ZCells[cellI]), small));
                flameletTable.find(ZCells[cellI], Zeta, chiCells[cellI]);
                heCells[cellI] = 0.0;
                tT = flameletTable.lookupT();
            }
            YCells[cellI] = flameletTable.lookupY(i);
            // heCells[cellI] += YCells[cellI]*thermo.composition().Hs(i, pCells[cellI], tT);
            heCells[cellI] += YCells[cellI]*thermo.composition().Hs(i, 1e5, tT);
        }
    }

    forAll(he.boundaryField(), patchI)
    {
        const fvPatchScalarField& pchi = chi.boundaryField()[patchI];
        const fvPatchScalarField& pvarZ = varZ.boundaryField()[patchI];
        const fvPatchScalarField& pZ = Z.boundaryField()[patchI];
        const fvPatchScalarField& pp = p.boundaryField()[patchI];

        fvPatchScalarField& phe = he.boundaryFieldRef()[patchI];

        forAll(Y, i)
        {
            fvPatchScalarField& pY = Y[i].boundaryFieldRef()[patchI];

            forAll(pY, faceI)
            {
                if (i == 0)
                {
                    Zeta = Foam::sqrt(pvarZ[faceI]/Foam::max(pZ[faceI]*(1 - pZ[faceI]), small));
                    flameletTable.find(pZ[faceI], Zeta, pchi[faceI]);
                    phe[faceI] = 0.0;
                    tT = flameletTable.lookupT();
                }
                pY[faceI] = flameletTable.lookupY(i);
                // phe[faceI] += pY[faceI]*thermo.composition().Hs(i, pp[faceI], tT);
                phe[faceI] += pY[faceI]*thermo.composition().Hs(i, 1e5, tT);
            }
        }
    }

    // he = thermo.he(p,T);
    thermo.correct();
    Info<< "T gas min/max : " << min(T).value() << ", "
        << max(T).value() << endl;
    Info<< "p min/max : " << min(p).value() << ", "
        << max(p).value() << endl;
    Info<< "Y oxidizer min/max : " << min(Y[oxidizerIndex]).value() << ", "
        << max(Y[oxidizerIndex]).value() << endl;
}
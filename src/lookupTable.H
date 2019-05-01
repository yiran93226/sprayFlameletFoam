{
    volScalarField& he = thermo.he();
    forAll(Z, i)
    {
        double Zeta = std::sqrt(varZ[i]/std::max(Z[i]*(1-Z[i]), 1e-10));
        flameletTable.find(Z[i], Zeta, chi[i]);
        double tT = flameletTable.lookupT();
        he[i] = 0;
        forAll(Y, j)
        {
            Y[j][i] = flameletTable.lookupY(j);
            he[i] += Y[j][i]*thermo.composition().Hs(j, 1e05, tT);
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
{
    volScalarField& he = thermo.he();
    forAll(Z.internalField(), i)
    {
        double Zeta = std::sqrt(varZ[i]/std::max(Z[i]*(1-Z[i]), 1e-10));
        flameletTable.find(Z[i], Zeta, chi[i]);
        T[i] = flameletTable.lookupT();
        forAll(Y, j)
        {
            Y[j][i] = flameletTable.lookupY(j);
        }
    }
    he = thermo.he(p,T);
    thermo.correct();
    Info<< "T gas min/max : " << min(T).value() << ", "
        << max(T).value() << endl;
    Info<< "Y oxidizer min/max : " << min(Y[oxidizerIndex]).value() << ", "
        << max(Y[oxidizerIndex]).value() << endl;
}
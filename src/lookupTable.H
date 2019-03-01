{
    volScalarField& he = thermo.he();
    forAll(Z.internalField(), i)
    {
        flameletTable.find(Z[i],chi[i]);
        T[i] = flameletTable.lookupT();
        forAll(Y, j)
        {
            if(j!=inertIndex)
                Y[j][i] = flameletTable.lookupY(j);
            else continue;
        }
    }
    he = thermo.he(p,T);
    thermo.correct();
    Info<< "T gas min/max : " << min(T).value() << ", "
        << max(T).value() << endl;
    Info<< "Y O2 min/max : " << min(Y[2]).value() << ", "
        << max(Y[2]).value() << endl;
}
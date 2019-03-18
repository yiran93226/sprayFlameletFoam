#include "tableSolver.H"
using namespace std;

int main()
{
    tableSolver ts;
    std::ofstream outfile("lookupResult.csv");

    size_t nsp = 231;
    std::vector<double> Z;
    std::vector<double> Yc;
    for (size_t i=0; i<1000; i++)
    {
        for (size_t j=0; j<1001; j++)
        {
            Z.push_back(double(j/1000.0));
            Yc.push_back(double(0.0003*i));
        }
    }
    std::vector<double> T(Z.size());
    std::vector<bool> ext(Z.size());
    std::vector<std::vector<double> > Y(nsp);
    for(size_t i=0; i<nsp; i++)
    {
        Y[i].resize(Z.size());
    }
    for(size_t i=0; i<Z.size(); i++)
    {
        ts.find(Z[i], Yc[i]);
        T[i] = ts.lookupT();
        ext[i] = ts.extrapol();
        for(size_t j=0; j<nsp; j++)
        {
            Y[j][i] = ts.lookupY(j);
        }
    }

    for(size_t i=0; i<Z.size(); i++)
    {
        if (ext[i] == false)
        {
            outfile << Z[i] << "," << Yc[i] << "," << T[i] << ",";
            for(size_t j=0; j<nsp; j++)
            {
                outfile << Y[j][i];
                if (j!=nsp-1) outfile << ",";
            }
            outfile << std::endl;
        }
    }

    return 0;
}

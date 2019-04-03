#include "tableSolver.H"
using namespace std;

int main()
{
    tableSolver ts;
    std::ofstream outfile("lookupResult.csv");

    size_t nsp = 231;
    std::vector<double> Z;
    std::vector<double> chi;
    for (size_t i=1; i<1000; i++)
    {
        for (size_t j=0; j<1001; j++)
        {
            Z.push_back(double(j/1000.0));
            chi.push_back(double(i/10.0));
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
        ts.find(Z[i], chi[i]);
        T[i] = ts.lookupT();
        ext[i] = ts.extrapol();
        for(size_t j=0; j<nsp; j++)
        {
            Y[j][i] = ts.lookupY(j);
        }
    }

    for(size_t i=0; i<Z.size(); i++)
    {
        // if (ext[i] == false)
        if (true)
        {
            outfile << Z[i] << "," << chi[i] << "," << T[i] << ",";
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

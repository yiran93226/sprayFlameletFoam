#include "tableSolver.H"
using namespace std;

int main()
{
    std::ifstream sprayFlameFile("sprayTable0.csv");
    table sprayFlame(sprayFlameFile);
    tableSolver ts;
    std::ofstream outfile("lookupResult.csv");

    double nsp = 231;
    std::vector<double> Z = sprayFlame.getZ();
    std::vector<double> Yc = sprayFlame.getYc();
    std::vector<double> T(Z.size());
    std::vector<std::vector<double> > Y(nsp);
    for(size_t i=0; i<nsp; i++)
    {
        Y[i].resize(Z.size());
    }
    for(size_t i=0; i<Z.size(); i++)
    {
        ts.find(Z[i], Yc[i]);
        T[i] = ts.lookupT();
        for(size_t j=0; j<nsp; j++)
        {
            Y[j][i] = ts.lookupY(j);
        }
    }


    outfile << sprayFlame.getFirstLine() << std::endl;
    for(size_t i=0; i<Z.size(); i++)
    {
        outfile << Z[i] << "," << Yc[i] << "," << T[i] << ",";
        for(size_t j=0; j<nsp; j++)
        {
            outfile << Y[j][i];
            if (j!=nsp-1) outfile << ",";
        }
        outfile << std::endl;
    }

    return 0;
}

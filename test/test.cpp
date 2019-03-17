#include "table.H"
using namespace std;

int main()
{
    std::string filename;
    filename = "flameletTable.csv";
    std::ifstream flameletTableFile(filename);
    std::ifstream sprayFlameFile("sprayTable.csv");
    std::ofstream outfile("lookupResult.csv");
    table flameletTable(flameletTableFile);
    table sprayFlame(sprayFlameFile);

    double nsp = flameletTable.getNsp();
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
        flameletTable.find(Z[i]);
        T[i] = flameletTable.lookupT();
        for(size_t j=0; j<nsp; j++)
        {
            Y[j][i] = flameletTable.lookupY(j);
        }
    }

    outfile << flameletTable.getFirstLine() << std::endl;
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

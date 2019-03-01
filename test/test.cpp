#include "table.H"
#include "table.C"
using namespace std;

int main()
{
    double T;
    ifstream flameletTableFile("flameletTable.csv");
    table flameletTable(flameletTableFile);
    flameletTable.find(0.000196492,0.000510477);
    T = flameletTable.lookupT();
    cout << T << endl;
    return 0;
}

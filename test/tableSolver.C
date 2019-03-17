#include "tableSolver.H"

tableSolver::tableSolver()
{
    collectTables();
}

void tableSolver::collectTables()
{
    size_t n=0;
    while(true)
    {
        std::ifstream flameletTableFile("tablesYc/flameletTable_"+std::to_string(n)+".csv");
        if (!flameletTableFile) break;
        table* tTable = new table(flameletTableFile);
        tables_.push_back(tTable);
        n++;
    }
    tableNum_ = tables_.size();
}

void tableSolver::find(double Z, double Yc)
{
    std::vector<double> Ycs(tableNum_);
    std::map<double, size_t> YcMap;
    for (size_t i=0; i<tableNum_; i++)
    {
        tables_[i]->find(Z);
        Ycs[i] = tables_[i]->lookupYc();
        YcMap.insert({Ycs[i], i});
    }
    std::sort(Ycs.begin(), Ycs.end());
    double YcH, YcL;
    if ( (Yc > Ycs[1]) && (Yc < Ycs[tableNum_-1]) )
    {
        for (size_t i=1; i<tableNum_; i++)
        {
            if (Ycs[i] > Yc)
            {
                YcH = Ycs[i];
                YcL = Ycs[i-1];
                break;
            }
        }
        weightL_ = (YcH - Yc) / (YcH - YcL);
        weightH_ = (Yc - YcL) / (YcH - YcL);
    }
    else if (Yc <= Ycs[1])
    {
        YcH = Ycs[0];
        YcL = Ycs[0];
        weightL_ = 0.5;
        weightH_ = 0.5;
    }
    else
    {
        YcH = Ycs[tableNum_-1];
        YcL = Ycs[tableNum_-1];
        weightL_ = 0.5;
        weightH_ = 0.5;
    }
    positionL_ = YcMap[YcL];
    positionH_ = YcMap[YcH];
}
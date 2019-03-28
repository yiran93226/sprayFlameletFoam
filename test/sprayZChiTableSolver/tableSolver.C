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
        std::ifstream flameletTableFile("tablesChi/sprayTable"+std::to_string(n)+".csv");
        if (!flameletTableFile) break;
        table* tTable = new table(flameletTableFile);
        tables_.push_back(tTable);
        n++;
    }
    tableNum_ = tables_.size();
}

void tableSolver::find(double Z, double chi)
{
    std::vector<double> chis(tableNum_-1);
    std::map<double, size_t> chiMap;
    for (size_t i=0; i<tableNum_; i++)
    {
        tables_[i]->find(Z);
        // The last table is an extiguished one
        if (i!=tableNum_-1)
        {
            chis[i] = tables_[i]->lookupChi();
            chiMap.insert({chis[i], i});
        }
    }
    std::sort(chis.begin(), chis.end());
    double chiH, chiL;
    if ( (chi > chis[0]) && (chi < chis[tableNum_-2]) )
    {
        for (size_t i=0; i<tableNum_-1; i++)
        {
            if (chis[i] > chi)
            {
                chiH = chis[i];
                chiL = chis[i-1];
                break;
            }
        }
        weightL_ = (chiH - chi) / (chiH - chiL);
        weightH_ = (chi - chiL) / (chiH - chiL);
        positionL_ = chiMap[chiL];
        positionH_ = chiMap[chiH];
    }
    else if (chi <= chis[0])
    {
        chiH = chis[0];
        chiL = chis[0];
        weightL_ = 0.5;
        weightH_ = 0.5;
        positionL_ = chiMap[chiL];
        positionH_ = chiMap[chiH];
    }
    else
    {
        weightL_ = 0.5;
        weightH_ = 0.5;
        positionL_ = tableNum_-1;
        positionH_ = tableNum_-1;
    }
}
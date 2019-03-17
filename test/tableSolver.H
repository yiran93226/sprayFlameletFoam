#ifndef TABLESOLVER_H
#define TABLESOLVER_H

#include <map>
#include "table.H"

class tableSolver
{
public:
    // Constructor
    tableSolver();

    // Getters

    // Member functions
    void find(double Z, double Yc);

    double lookupT() const
    {
        return tables_[positionH_]->lookupT() * weightH_ +
                tables_[positionL_]->lookupT() * weightL_;
    }

    double lookupY(size_t i) const
    {
        return tables_[positionH_]->lookupY(i) * weightH_ + 
                tables_[positionL_]->lookupY(i) * weightL_;
    }


private:
    // Private member functions
    void collectTables();
 
    // Private member data
    double weightL_;
    double weightH_;
    size_t positionL_;
    size_t positionH_;

    size_t tableNum_;
    std::vector<table*> tables_;
};

#endif
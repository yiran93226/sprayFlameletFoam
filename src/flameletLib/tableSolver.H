#ifndef TABLESOLVER_H
#define TABLESOLVER_H

#include <map>
#include <cmath>
#include "table.H"

class tableSolver
{
public:
    // Constructor
    tableSolver(size_t zetaIndex);

    // Getters

    // Member functions
    void find(double Z, double Yc);

    size_t nsp() const
    {
        return nsp_;
    }

    void bdSet(double Z)
    {
        tables_[0]->find(Z);
        positionL_ = 0;
        positionH_ = 0;
        weightH_ = 0.5;
        weightL_ = 0.5;
    }

    double lookupT() const
    {
        return tables_[positionH_]->lookupT() * weightH_ +
               tables_[positionL_]->lookupT() * weightL_;
    }

    double lookupOmegaYc() const
    {
        return tables_[positionH_]->lookupOmegaYc() * weightH_ +
               tables_[positionL_]->lookupOmegaYc() * weightL_;
    }

    double lookupY(size_t i) const
    {
        return tables_[positionH_]->lookupY(i) * weightH_ +
               tables_[positionL_]->lookupY(i) * weightL_;
    }

private:
    // Private member functions
    void collectTables(size_t zetaIndex);
 
    // Private member data
    double weightL_;
    double weightH_;
    size_t positionL_;
    size_t positionH_;

    size_t nsp_;
    size_t tableNum_;
    std::vector<table*> tables_;
};

#endif

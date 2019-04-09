#ifndef SELECTOR_H
#define SELECTOR_H

#include <map>
#include "tableSolver.H"

class selector
{
public:
    // Constructor
    selector(std::string type);

    // Member functions
    void find(double Z, double Zeta, double chi);

    double lookupT() const
    {
        return tableSolvers_[positionH_]->lookupT() * weightH_ +
                tableSolvers_[positionL_]->lookupT() * weightL_;
    }

    double lookupY(size_t i) const
    {
        return tableSolvers_[positionH_]->lookupY(i) * weightH_ + 
                tableSolvers_[positionL_]->lookupY(i) * weightL_;
    }

private:
 
    // Private member data
    double weightL_;
    double weightH_;
    size_t positionL_;
    size_t positionH_;
    double Z_;
    double chi_;
    double minZeta_;
    double maxZeta_;
    double interZeta_;
    size_t num_;
    std::vector<tableSolver*> tableSolvers_;
    std::vector<double> Zetas_;
};

#endif
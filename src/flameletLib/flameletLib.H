#ifndef FLAMELETLIB_H
#define FLAMELETLIB_H

#include "tableSolver.H"

class flameletLib
{
public:
    // Constructor
    flameletLib(std::string type);

    // Member functions
    void find(double Z, double Zeta, double Yc);

    double lookupT() const
    {
        return tableSolvers_[positionH_]->lookupT() * weightH_ +
                tableSolvers_[positionL_]->lookupT() * weightL_;
    }

    double lookupOmegaYc() const
    {
        return tableSolvers_[positionH_]->lookupOmegaYc() * weightH_ +
                tableSolvers_[positionL_]->lookupOmegaYc() * weightL_;
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
    double Yc_;
    double minZeta_;
    double maxZeta_;
    double interZeta_;
    size_t num_;
    std::vector<tableSolver*> tableSolvers_;
    std::vector<double> Zetas_;
};

#endif
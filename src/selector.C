#include "selector.H"

selector::selector(std::string type)
{
    if (type == "RAS")
    {
        num_ = 10;
    }
    else num_ = 1;
    tableSolvers_.resize(num_);
    Zetas_.resize(num_);
    for (size_t i=0; i<num_; i++)
    {
        Zetas_[i] = i*0.11;
        tableSolver* tTableS = new tableSolver(i);
        std::cout << "#FLAMELETTABLE Constructed" << std::endl;
        tableSolvers_[i] = tTableS;
    }
    minZeta_ = Zetas_[0];
    maxZeta_ = Zetas_[num_-1];
}

void selector::find(double Z, double Zeta, double chi)
{
    if ( (Zeta > minZeta_) && (Zeta < maxZeta_) )
    {
        interZeta_ = Zeta;
        for (size_t i=0; i<num_; i++)
        {
            if (Zetas_[i] > interZeta_)
            {
                positionL_ = i-1;
                positionH_ = i;
                break;
            }
        }
        weightL_ = (Zetas_[positionH_] - interZeta_) / (Zetas_[positionH_] - Zetas_[positionL_]);
        weightH_ = (interZeta_ - Zetas_[positionL_]) / (Zetas_[positionH_] - Zetas_[positionL_]);
    }
    else if(Zeta <= minZeta_)
    {
        interZeta_ = minZeta_;
        positionL_ = 0;
        positionH_ = 0;
        weightL_ = 0.5;
        weightH_ = 0.5;
    }
    else
    {
        interZeta_= maxZeta_;
        positionL_ = num_-1;
        positionH_ = num_-1;
        weightL_ = 0.5;
        weightH_ = 0.5;
    }
    tableSolvers_[positionL_]->find(Z, chi);
    tableSolvers_[positionH_]->find(Z, chi);
}
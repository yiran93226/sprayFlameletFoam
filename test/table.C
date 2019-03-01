#include "table.H"


table::table
(
    std::ifstream& tableFile
)
:
    tableFile_(tableFile)
{
    read();
    positions_.resize(4);
}


double table::lookupT() const
{
    double T;
    double chiL;
    double chiH;
    double TL;
    double TH;
    if(positions_[3]==0)
    {
        T = ZKey_*480.0 + (1.0-ZKey_)*800;
    }
    else
    {
        chiL = ( (Z_[positions_[1]] - ZKey_)*chi_[positions_[0]]
                + (ZKey_ - Z_[positions_[0]])*chi_[positions_[1]] )
                / (Z_[positions_[1]] - Z_[positions_[0]]);
        chiH = ( (Z_[positions_[3]] - ZKey_)*chi_[positions_[2]]
                + (ZKey_ - Z_[positions_[2]])*chi_[positions_[3]] )
                / (Z_[positions_[3]] - Z_[positions_[2]]);
        TL = ( (Z_[positions_[1]] - ZKey_)*T_[positions_[0]]
                + (ZKey_ - Z_[positions_[0]])*T_[positions_[1]] )
                / (Z_[positions_[1]] - Z_[positions_[0]]);
        TH = ( (Z_[positions_[3]] - ZKey_)*T_[positions_[2]]
                + (ZKey_ - Z_[positions_[2]])*T_[positions_[3]] )
                / (Z_[positions_[3]] - Z_[positions_[2]]);
        T = ( (chiH - chiKey_)*TL + (chiKey_ - chiL)*TH ) / (chiH - chiL);
    }
    return T;
}


double table::lookupY(size_t i) const
{
    double Y;
    double chiL;
    double chiH;
    double YL;
    double YH;
    if(positions_[3]==0)
    {
        Y = ZKey_*480.0 + (1.0-ZKey_)*800;
    }
    else
    {
        chiL = ( (Z_[positions_[1]] - ZKey_)*chi_[positions_[0]]
                + (ZKey_ - Z_[positions_[0]])*chi_[positions_[1]] )
                / (Z_[positions_[1]] - Z_[positions_[0]]);
        chiH = ( (Z_[positions_[3]] - ZKey_)*chi_[positions_[2]]
                + (ZKey_ - Z_[positions_[2]])*chi_[positions_[3]] )
                / (Z_[positions_[3]] - Z_[positions_[2]]);
        YL = ( (Z_[positions_[1]] - ZKey_)*Y_[i][positions_[0]]
                + (ZKey_ - Z_[positions_[0]])*Y_[i][positions_[1]] )
                / (Z_[positions_[1]] - Z_[positions_[0]]);
        YH = ( (Z_[positions_[3]] - ZKey_)*Y_[i][positions_[2]]
                + (ZKey_ - Z_[positions_[2]])*Y_[i][positions_[3]] )
                / (Z_[positions_[3]] - Z_[positions_[2]]);
        Y = ( (chiH - chiKey_)*YL + (chiKey_ - chiL)*YH ) / (chiH - chiL);
    }
    return Y;
}


void table::find(double Z, double chi)
{
    std::vector<size_t> interPos;
    double interChi;

    if(Z <= 0) Z = 0;
    if(Z >= 1) Z = 1;
    ZKey_ = Z;
    chiKey_ = chi;
    for(size_t i=0;i<Z_.size();i++)
    {
        if(Z_[i]<=Z)
        {
            interPos.push_back(i);
            for(size_t j=i;j<Z_.size();j++)
                if(Z_[j]>Z) 
                {
                    i=j-1;
                    break;
                }
        }
    }
    for(size_t i=0;i<interPos.size();i++)
    {
        // Linear interpolate chi
        interChi = ( (Z_[interPos[i]] - Z)*chi_[interPos[i]-1]
                    + (Z - Z_[interPos[i]-1])*chi_[interPos[i]] )
                    / (Z_[interPos[i]] - Z_[interPos[i]-1]);
        if(interChi>chi)
        {
            positions_[3] = interPos[i];
            positions_[2] = interPos[i]-1;
            positions_[1] = interPos[i-1];
            positions_[0] = interPos[i-1]-1;
            break;
        }
        if(i==interPos.size()-1)
        {
            positions_[3] = 0;
            positions_[2] = 0;
            positions_[1] = 0;
            positions_[0] = 0;
        }
    }
}


void table::read()
{
    std::string line, str;
    std::getline(tableFile_, str); // The first line
    nColumn_ = std::count(str.begin(), str.end(), ',') + 1;
    std::cout << "#Number of columns:\t" << nColumn_ << std::endl;
    Y_.resize(nColumn_-3);
    while(std::getline(tableFile_,line))
    {
        size_t n = 0;
        std::istringstream buffer(line);
        while(std::getline(buffer, str, ',')) // Comma separated values
        {
            if(n > 2)
                Y_[n-3].push_back(std::stod(str));
            else if(n == 0)
                Z_.push_back(std::stod(str));
            else if(n == 1)
                chi_.push_back(std::stod(str));
            else
                T_.push_back(std::stod(str));
            n++;
        }
    }
    std::cout << "#Number of species:\t" << Y_.size() << std::endl;
    std::cout << "#Length of Z:\t" << Z_.size() << std::endl;
    std::cout << "#Length of Yi:\t" << Y_[0].size() << std::endl;
    std::cout << "#Length of T:\t" << T_.size() << std::endl;
}

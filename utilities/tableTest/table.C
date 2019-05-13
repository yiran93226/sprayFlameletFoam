#include "table.H"

table::table
(
    std::ifstream& tableFile
)
:
    tableFile_(tableFile)
{
    read();
}



void table::find(double Z)
{
    if ( (Z > minZ_) && (Z < maxZ_) )
    {
        interZ_ = Z;
        for (size_t i=0; i<lenZ_; i++)
        {
            if (Z_[i] > interZ_)
            {
                positionL_ = i-1;
                positionH_ = i;
                break;
            }
        }
        weightL_ = (Z_[positionH_] - interZ_) / (Z_[positionH_] - Z_[positionL_]);
        weightH_ = (interZ_ - Z_[positionL_]) / (Z_[positionH_] - Z_[positionL_]);
    }
    else if(Z <= minZ_)
    {
        interZ_ = minZ_;
        positionL_ = 0;
        positionH_ = 0;
        weightL_ = 0.5;
        weightH_ = 0.5;
    }
    else
    {
        interZ_= maxZ_;
        positionL_ = lenZ_-1;
        positionH_ = lenZ_-1;
        weightL_ = 0.5;
        weightH_ = 0.5;
    }
    interChi_ = chi_[positionH_] * weightH_ + chi_[positionL_] * weightL_;
}


void table::read()
{
    std::string line, str;
    std::getline(tableFile_, firstLine_); // The first line
    nColumn_ = std::count(firstLine_.begin(), firstLine_.end(), ',') + 1;
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
    lenZ_ = Z_.size();
    maxZ_ = Z_[lenZ_-1];
    minZ_ = Z_[0];
}

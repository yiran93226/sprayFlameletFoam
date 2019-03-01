#ifndef TABLE_H
#define TABLE_H

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>

class table
{
public:
    // Constructor
    table(std::ifstream& tableFile);

    // Member functions
    double lookupT() const;

    void find(double Z, double chi);

private:
    // Private member functions
    void read();

    // Private member data
    std::vector<double> Z_;
    std::vector<double> chi_;
    std::vector<double> T_;
    std::vector<std::vector<double> > Y_;

    double ZKey_;
    double chiKey_;
    std::vector<int> positions_;

    int nColumn_;
    std::ifstream& tableFile_;
};

#endif
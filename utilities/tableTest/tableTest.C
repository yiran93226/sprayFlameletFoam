#include <iostream>
#include "scalar.H"
#include "List.H"
#include "OFstream.H"
#include "Gamma.h"
#include "table.H"

using namespace Foam;

void betaPDFIntegration(const scalar& Zeta, const List<List<scalar> > singleData_, List<List<scalar> >& integratedData_);
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
int main(int argc, char *argv[])
{
   const List<scalar> Zetas = {0.11, 0.22, 0.33, 0.44, 0.55, 0.66, 0.77, 0.88, 0.99};
    
   size_t nTable = 0;
   while(true)
   {
      std::ifstream inTableFile("pdfTable_pilot/Zeta_0/flameletTable_"+std::to_string(nTable)+".csv");
      if (!inTableFile) break;

      table t_0(inTableFile);
      std::vector<double> Z = t_0.getZ();
      std::vector<double> Chi = t_0.getChi();
      std::vector<double> T = t_0.getT();
      std::vector<std::vector<double> > Y = t_0.getY();
      double nsp = t_0.getNsp();

      std::string firstLine = t_0.getFirstLine();
    
      //  Info<< "species number = " << nsp << endl;

      List<scalar> single(Z.size(),1.0);
      List<List<scalar> > singleData_(3+nsp,single);

      // data assignment from the csv file
      for(int j = 0; j < Z.size(); j++)
      {
         singleData_[0][j] = Z[j];
      }

      for(int j = 0; j < Chi.size(); j++)
      {
         singleData_[1][j] = Chi[j];
      }

      for(int j = 0; j < T.size(); j++)
      {
          singleData_[2][j] = T[j];
      }

      for(int i = 0; i < nsp; i++)
      {
         for(int j = 0; j < Y[i].size(); j++)
         {
            singleData_[i+3][j] = Y[i][j];
         }       
      }

      List<List<scalar> > integratedData_(singleData_);
    
    
      for(int i = 0; i < Zetas.size(); i++)
      {
         const scalar Zeta = Zetas[i];
         betaPDFIntegration(Zeta, singleData_, integratedData_);

         // output csv file
         OFstream outTableFile("pdfTable_pilot/Zeta_"+std::to_string(i+1)+"/flameletTable_"+std::to_string(nTable)+".csv");
         if (!outTableFile) {
            std::cerr << "Unable to save data!\n";
         }
         else
         {
            outTableFile << firstLine << endl;
            for (int j = 0; j < Z.size(); j++)
            {
               for(int i = 0; i < integratedData_.size(); i++)
               {
                  if (i!= integratedData_.size()-1) {
                     outTableFile << integratedData_[i][j] << ',' ;
                  }
                  else
                  {
                     outTableFile << integratedData_[i][j] ;
                  }  
               }      
               outTableFile <<endl;
            }       
         }
      }  

      nTable++;
   }
   Info<< "nTable = " << nTable << endl;

   //  std::ifstream flameletTableFile("flameletTable_0.csv");

    
    
    
    
   //  print on screen
   //  for (int j = 0; j < Z.size(); j++)
   //  {
   //  //    Info<< "i = " << i << endl;
   //     for(int i = 0; i < integratedData_.size(); i++)
   //     {
   //      //    Info<< "j = " << j << endl;
   //         Info<< integratedData_[i][j] << ',' ;
   //     }      
   //     Info<< endl;
   //  }

    
    
}


void betaPDFIntegration(const scalar& Zeta, const List<List<scalar> > singleData_, List<List<scalar> >& integratedData_)
{
if (Zeta != 0)
{
   List<scalar> Z_(integratedData_[0]);
   List<scalar> varZ_(integratedData_[0].size(), 0.0);
   List<scalar> pdfAlpha(integratedData_[0].size(), 0.0);
   List<scalar> pdfBeta(integratedData_[0].size(), 0.0);
   List<scalar> PDF(integratedData_[0].size(), 0.0);

   for (int i=0;i<integratedData_[0].size();i++)
   {
      varZ_[i] = sqr(Zeta) * (Z_[i]*(1.0 - Z_[i]));

      if (varZ_[i] > 1e-7)
      {
          pdfAlpha[i] = Z_[i] * ((Z_[i]*(1.0-Z_[i]))/varZ_[i]-1);
          pdfBeta[i] = (1.0 - Z_[i]) * ((Z_[i]*(1.0-Z_[i]))/varZ_[i]-1);

          // Limit alpha and beta but keep their ratio
          if (pdfAlpha[i] > 500)
          {
             pdfBeta[i] = pdfBeta[i]/pdfAlpha[i] * 500;
             pdfAlpha[i] = 500;
          }
          if (pdfBeta[i] > 500)
          {
             pdfAlpha[i] = pdfAlpha[i]/pdfBeta[i] * 500;
             pdfBeta[i] = 500;
          }

          int    gridPoints = 250;
          List<long double> hZ_(gridPoints, 0.0);
          List<long double> helpZ_(gridPoints, 0.0);
          List<long double> delta_(gridPoints, 0.0);

          if ((pdfAlpha[i] > 1) && (pdfBeta[i] > 1))
          {
             // Allocation of Z for PDF integration
             scalar Zmax = 0;
             int    n1 = 0;
             int    n2 = 0;
             PDF.clear();
             PDF.resize(Z_.size(), 0.0);
             
             for (int j=0;j<Z_.size();j++)
             {
                PDF[j] = std::pow(Z_[j],(pdfAlpha[i]-1.0)) * std::pow((1.0 - Z_[j]),(pdfBeta[i]-1.0)) * Gamma(pdfAlpha[i] + pdfBeta[i]) / (min(Gamma(pdfAlpha[i]),1e17)*min(Gamma(pdfBeta[i]),1e17));
                if (PDF[j] > PDF[max(j-1, 0)]) Zmax = Z_[j];
             }

             if(pdfAlpha[i]/pdfBeta[i] <= 0.5)
             {
                n1 = 0.2*gridPoints;
                n2 = 0.8*gridPoints+1;
             }
             else if (pdfAlpha[i]/pdfBeta[i] >= 2)
             {
                n1 = 0.8*gridPoints;
                n2 = 0.2*gridPoints+1;
             }
             else
             {
                n1 = 0.5*gridPoints;
                n2 = 0.5*gridPoints+1;
             }

             //  Allocate Z for 0 < Z < Zmax
             scalar ex1 = 0.9;
             delta_[0] = (1.0 - ex1)/(1.0 - pow(ex1,(n1-1)));

             for (int j=0; j<n1; j++)
             {
                delta_[j] = pow(ex1,j) * delta_[0];
                hZ_[j+1] = hZ_[j] + delta_[j];
             }
             for (int j=1; j<n1; j++)
             {
                hZ_[j] *= Zmax;
             }

             // Allocate Z for Zmax < Z < 1
             scalar ex2 = 1.1;
             delta_[0] = (1.0-ex2)/(1.0-pow(ex2,(n2-1)));

             for (int j=0; j<n2-1; j++)
             {
               delta_[j] = pow(ex2,j) * delta_[0];
               helpZ_[j+1] = helpZ_[j] + delta_[j];
             }
             for (int j=0;j<n2;j++)
             {
                helpZ_[j] *= (1.0 - Zmax);
             }
             for (int j=0;j<n2-1;j++)
             {
                hZ_[n1+j] = hZ_[n1-1]+helpZ_[j];
             }

             // Scaling
             for (int j=0;j<gridPoints;j++)
             {
                hZ_[j] /= hZ_[gridPoints-1];
             }

             // Calculate BetaPDF
             PDF.clear();
             PDF.resize(hZ_.size(), 0.0);

             for (int j=0;j<hZ_.size();j++)
             {
                PDF[j] = (std::pow(hZ_[j],(pdfAlpha[i]-1e0))) * std::pow((1e0 - hZ_[j]),(pdfBeta[i]-1e0)) * Gamma(pdfAlpha[i] + pdfBeta[i])/(min(Gamma(pdfAlpha[i]),1e17)*min(Gamma(pdfBeta[i]),1e17));
             }
          }

          else if ((pdfAlpha[i] <= 1) && (pdfBeta[i] > 1))
          {
             // PDF Singularity at Z = 0
             // Allocation of Z for PDF integration
             int    n1 = 0;
             int    n2 = 0;

             if (pdfAlpha[i]/pdfBeta[i] > 0.5)
             {
                scalar Zmax = 0.5;
                scalar ex1 = 1.1;
                n1 = 0.7*gridPoints;
                delta_[0] = (1.0 - ex1)/(1.0 - pow(ex1,(n1-1)));

                // Allocate Z for 0 < Z < Zmax
                for (int j=0; j<n1; j++)
                {
                   delta_[j] = pow(ex1,j) * delta_[0];
                   hZ_[j+1] = hZ_[j] + delta_[j];
                }
                for (int j=1; j<n1; j++)
                {
                   hZ_[j] *= Zmax;
                }

                // Allocate Z for Zmax < Z < 1
                scalar ex2 = 1.1;
                n2 = 0.3*gridPoints+1;
                delta_[0] = (1.0-ex2)/(1.0-pow(ex2,(n2-1)));

                for (int j=0; j<n2-1; j++)
                {
                   delta_[j] = pow(ex2,j) * delta_[0];
                   helpZ_[j+1] = helpZ_[j] + delta_[j];
                }
                for (int j=0;j<n2;j++)
                {
                   helpZ_[j] *= (1.0 - Zmax);
                }
                for (int j=0;j<n2-1;j++)
                {
                   hZ_[n1+j] = hZ_[n1-1]+helpZ_[j];
                }
             }

             else
	         {
                scalar ex2 = 1.05;
                delta_[0] = (1.0 - ex2)/(1.0 - pow(ex2,(gridPoints-1)));
                for (int j=0; j<gridPoints-1; j++)
                {
                   delta_[j] = pow(ex2,j) * delta_[0];
                   hZ_[j+1] = hZ_[j] + delta_[j];
                }
             }

             // Scaling
             for (int j=0;j<gridPoints;j++)
             {
                hZ_[j] /= hZ_[gridPoints-1];
             }

             // Calculate BetaPDF
             PDF.clear();
             PDF.resize(hZ_.size(), 0.0);
             for (int j=1;j<hZ_.size();j++)
             {
                PDF[j] = std::pow(hZ_[j],(pdfAlpha[i]-1e0)) * std::pow((1e0 - hZ_[j]),(pdfBeta[i]-1.0)) * min(Gamma(pdfAlpha[i] + pdfBeta[i]),1e17)/(min(Gamma(pdfAlpha[i]),1e17)*min(Gamma(pdfBeta[i]),1e17));
             }
             PDF[0] = 1.5 * PDF[1] / pdfAlpha[i];

          }

          else if ((pdfAlpha[i] > 1) && (pdfBeta[i] <= 1))
          {
          // PDF Singularity at Z = 1
          // Allocation of Z for PDF integration

          int    n1 = 0;
          int    n2 = 0;

          if (pdfAlpha[i]/pdfBeta[i] < 2)
          {
             scalar Zmax = 0.5;
             scalar ex1 = 1.1;
             n1 = 0.3*gridPoints;
             delta_[0] = (1.0 - ex1)/(1.0 - pow(ex1,(n1-1)));

             // Allocate Z for 0 < Z < Zmax
             for (int j=0; j<n1; j++)
             {
                delta_[j] = pow(ex1,j) * delta_[0];
                hZ_[j+1] = hZ_[j] + delta_[j];
             }
             for (int j=1; j<n1; j++)
             {
                hZ_[j] *= Zmax;
             }

             // Allocate Z for Zmax < Z < 1
             scalar ex2 = 0.9;
             n2 = 0.7*gridPoints+1;
             delta_[0] = (1.0-ex2)/(1.0-pow(ex2,(n2-1)));

             for (int j=0; j<n2-1; j++)
             {
                delta_[j] = pow(ex2,j) * delta_[0];
                helpZ_[j+1] = helpZ_[j] + delta_[j];
             }
             for (int j=0;j<n2;j++)
             {
                helpZ_[j] *= (1.0 - Zmax);
             }
             for (int j=0;j<n2-1;j++)
             {
                hZ_[n1+j] = hZ_[n1-1]+helpZ_[j];
             }
          }

          else
          {
             scalar ex1 = 0.95;
             delta_[0] = (1.0 - ex1)/(1.0 - pow(ex1,(gridPoints-1)));

             for (int j=0; j<gridPoints-1; j++)
             {
                 delta_[j] = pow(ex1,j) * delta_[0];
                 hZ_[j+1] = hZ_[j] + delta_[j];
             }
          }

          // Scaling
          for (int j=0;j<gridPoints;j++)
          {
             hZ_[j] /= hZ_[gridPoints-1];
          }

          // Calculate BetaPDF
          PDF.clear();
          PDF.resize(hZ_.size(), 0.0);

          for (int j=0;j<hZ_.size()-1;j++)
          {
              PDF[j] = std::pow((hZ_[j]),(pdfAlpha[i]- 1e0)) * std::pow((1e0 - hZ_[j]),(pdfBeta[i]-1.0)) * min(Gamma(pdfAlpha[i] + pdfBeta[i]),1e17)/(min(Gamma(pdfAlpha[i]),1e17)*min(Gamma(pdfBeta[i]),1e17));
          }
          PDF[gridPoints - 1] = 1.5 * PDF[gridPoints - 2] / pdfBeta[i];

       }

       else if ((pdfAlpha[i] <= 1) && (pdfBeta[i] <= 1))
       {
          // PDF Singularity at Z = 1 and Z = 0
          // Allocation of Z for PDF integration
          int    n1 = 0;
          int    n2 = 0;
          scalar Zmax = 0.5;
          scalar ex1 = 1.1;
          n1 = 0.5*gridPoints;
          delta_[0] = (1.0 - ex1)/(1.0 - pow(ex1,(n1-1)));

          // Allocate Z for 0 < Z < Zmax
          for (int j=0; j<n1; j++)
          {
              delta_[j] = pow(ex1,j) * delta_[0];
              hZ_[j+1] = hZ_[j] + delta_[j];
          }
          for (int j=1; j<n1; j++)
          {
             hZ_[j] *= Zmax;
          }

          // Allocate Z for Zmax < Z < 1
          scalar ex2 = 0.9;
          n2 = 0.5*gridPoints+1;
          delta_[0] = (1.0-ex2)/(1.0-pow(ex2,(n2-1)));

          for (int j=0; j<n2-1; j++)
          {
              delta_[j] = pow(ex2,j) * delta_[0];
              helpZ_[j+1] = helpZ_[j] + delta_[j];
          }
          for (int j=0;j<n2;j++)
          {
              helpZ_[j] *= (1.0 - Zmax);
          }
          for (int j=0;j<n2-1;j++)
          {
             hZ_[n1+j] = hZ_[n1-1]+helpZ_[j];
          }

          // Scaling
          for (int j=0;j<gridPoints;j++)
          {
             hZ_[j] /= hZ_[gridPoints-1];
          }

          // Calculate BetaPDF
          PDF.clear();
          PDF.resize(hZ_.size(), 0.0);

          for (int j=1;j<hZ_.size()-1;j++)
          {
             PDF[j] = std::pow(hZ_[j],(pdfAlpha[i]- 1e0)) * std::pow((1e0 - hZ_[j]),(pdfBeta[i]-1.0)) * min(Gamma(pdfAlpha[i] + pdfBeta[i]),1e17)/(min(Gamma(pdfAlpha[i]),1e17)*min(Gamma(pdfBeta[i]),1e17));
          }
          PDF[gridPoints - 1] = 1.5 * PDF[gridPoints - 2] / pdfBeta[i];
          PDF[0] = 1.5 * PDF[1] / pdfAlpha[i];
       }

       // Calculate the area of the PDF for scaling
       scalar intPDF = 0;
       for (int j=1;j<hZ_.size();j++)
       {
    	  intPDF += (hZ_[j-1] - hZ_[j]) * (PDF[j-1] + PDF[j])/2;
       }

       // Interpolate singleData entries to the new mixture fraction space
       List<scalar> hY_(gridPoints, 0.0);
       scalar intY = 0;

       for (int j=0;j<integratedData_.size();j++)
       {
          hY_ = 0.0;
          hY_[0] = singleData_[j][0];
          hY_[hY_.size()-1] = singleData_[j][singleData_[j].size()-1];
          intY = 0;

          for (int k=1;k<hZ_.size()-1;k++)
          {
             int ubZ = 0;
             for (int l=0;l<Z_.size();l++)
             {
                ubZ = l;
                if (hZ_[k] < Z_[l])
                break;
             }
             int lbZ = ubZ -1;

             // Interpolation to hZ space
             hY_[k] = (singleData_[j][ubZ] - singleData_[j][lbZ])/max(Z_[ubZ] - Z_[lbZ], SMALL) * (hZ_[k] - Z_[lbZ]) + singleData_[j][lbZ];
             // PDF Integration using the trapezoidal rule
             intY += (hZ_[k-1] - hZ_[k]) * (hY_[k-1]*PDF[k-1] + hY_[k]*PDF[k])/(2.0 * intPDF);
          }

          // Special treatment for the boundaries
          intY += (hZ_[hZ_.size()-2] - hZ_[hZ_.size()-1]) * (hY_[hZ_.size()-2]*PDF[hZ_.size()-2] + hY_[hZ_.size()-1]*PDF[hZ_.size()-1])/(2.0 * intPDF);
          if (i != 0 && i != integratedData_[0].size()-1 && j != 0)
             integratedData_[j][i] = intY;
       }
     }
   }
 }
}
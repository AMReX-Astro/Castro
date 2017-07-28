#include<iostream>
#include<cmath>
#include<fstream>
#include<iomanip>
using namespace std;

int main(){
  ofstream data_output;
  double P1, Pc,dP, Lad, Lin;
  double R , z , den1 , grav , T1, Td;
  double exp1, exp2;
  double X = 1.0;  
  double P2, dz, den2, T2;
  int alpha=1, beta=0;//, count = 0;
 
 
  //Declaration Block
  Td = 1500.0;
  Pc = 1.0e9;
  dP = 2e3;
  Lad = (2.0/7.0);
  Lin = (1.0/2.0);
  grav = -1.0e3;
  R = (1.38e-16)/(2.34*1.66e-24); 

  exp1=1 + alpha;
  exp2=(1.0/(4.0-beta));

  P1 = 1.0e9;
  z=0;
 
  //Open data file
  data_output.open ("newmodelcpp.hse",ios::out);
 
  data_output << "# npts = 499999" << '\n' << "# num of variables = 4" << '\n' << "# density" << '\n' << "# temperature" <<'\n' << "# pressure" << '\n' << "# X " << '\n';

  //Iterate over values until density is 0 using constant dP
  do{

    P2 = P1-dP;
   
    T1 = Td * pow(1.0 + (Lad / (Lin - Lad))*pow((P1/Pc),exp1),exp2);
    T2 = Td * pow(1.0 + (Lad / (Lin - Lad))*pow((P2/Pc),exp1),exp2);
    
    den1 = (P1/(R*T1));
    den2 = (P2/(R*T2));

    dz = - (P1-P2)* (1.0/grav)*(1.0/2.0)*((1.0/den1)+(1.0/den2));

    
    data_output << setprecision(10) << z << ' ' << den1 << ' ' << T1 << ' ' << P1 << ' '<<  X << '\n';

    z = z + dz;
    P1 = P2;
    /* cout << count << '\n';
       count++;*/
 
  }while(den1 > 1e-7);

 
  //Close data file
  data_output.close();
  
  return 0;
}

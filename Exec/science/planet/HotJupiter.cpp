//July 25

#include<iostream>
#include<cmath>
#include<fstream>
#include<iomanip>

using namespace std;

int main(){
  ofstream data_output,plot_output;
  double P1, Pc,Pd,dP, Lad, Lin,P_char,T_char,Tc,P_bottom,P_top,den_top,z_top;
  double R , z , den1 , grav , T1, Td, z0;
  double exp1, exp2,exp3,dP_factor,number_cell;
  double X = 1.0,buffer_height,T_buffer,den_buffer,dP_factor2;
  double P2, dz, den2, T2,error_expected,error,P2_temp;
    int alpha=1, beta=0, count = 0,count_cell,count2=0 , count3=0;
 // static double P[10000000],T[10000000],den[10000000],y[10000000];
  double P[5000],T[5000],den[5000],y[5000];
  //std::cout.precision(20);
 
  //Declaration Block
  Td = 1500.0e0;       //T_deep : the temperature at the top of the atmosphere[K]
 // Pc = 1.0e9;        //the pressure at the convective-radiation boundary [10^-6 bar]
  P_char = 1.0e6;
  T_char = 250.0e0;
  dP = 1.0e+3;         //dP defaul value
  Lad = (2.0e0/7.0e0); //the adiabatic index of an ideal diatomic gas
  Lin = (1.0e0/2.0e0);
  grav = -1.0e3;       //constant gravity [cm/s^2]
  R = (1.3806488e-16)/2.34e0*6.02214129e23;   //ideal gas constant (defined with mean molecular mass)

  exp1=1 + alpha;
  exp2=(1.0e0/(4.0e0-beta));
  exp3=exp1/exp2;
  //Defining the computation box
  z_top=1.0e10;                         // the size of the computation box
  number_cell=512.0e0;                  // the number of cells

  
  
  //Defining the properties of the planet atmosphere
  den_top=1.e-8 ;                       // the density limit at the top of the atmosphere
  P_top=den_top*R*Td ;                  // the pressure limit at the density limit
  Tc = Td * pow(Lin/(Lin-Lad),exp2);                 //the temperature at the radiative-convective boundary (RCB)
  Pc = P_char * pow(Tc/T_char,1.0e0/Lad);            //the pressure at the RCB
  Pd = pow((Lin-Lad)/(2.0e0*Lin-Lad),1.0e0/exp1)*Pc; //the characteristic pressure defining the radiative zone

  
  
  //Defining the properties of the buffer
  den_buffer=1.e-10;                    // the density in the buffer above the atmosphere
  T_buffer=1.0e2;                       // the temperature in the buffer above the atmosphere
  buffer_height=z_top*1.0e0/3.0e0;      // the height at which the buffer and the top of the atmosphere meet

  //Initializing
  P1=P_top;
  dP_factor=1.e5;
  dP_factor2=1.e-10*dP_factor;
  z=z_top;//+0.5*z_top/number_cell;
  P1=P_top;
  count=0;
  error_expected=1.e2/z_top;
  z0=100.0e0;
  
  
  
  data_output.open ("newmodelcpp.hse",ios::out);
  plot_output.open ("initialmodel_1500_T1_250.txt",ios::out);
  

  do{
    if(z>buffer_height){
      den1 = den_top;
      den2 = den_top;
      T1   = T_buffer;
      T2   = T_buffer;
      P1   = R*T_buffer*den1;
      P2   = R*T_buffer*den2;
      T[count2] = T1;
      P[count2] = P1;
      den[count2] = den1;
      y[count2] = z;
      z = z_top-z_top/number_cell*(float(count_cell)+0.5e0);
      count2++;
      count++;
      count_cell++;
      P1 = R*Td*den1;
      P2 = R*Td*den2;
      cout<<"In buffer" << ' '<< count_cell << ' ' << z<<' '<< dz/(z_top/number_cell) << ' ' << error << '\n';
    }
    else{
      if(Pc<P1){
        dP_factor = 1.e12;
        dP_factor2 = 1.e-12*dP_factor;
      }
        dP = P1/dP_factor;
        P2 = P1+dP;
        P2_temp = 0;
      
      if(Pc>P1){
        T1 = Td * pow(1.0e0 + (Lad / (Lin - Lad))*pow((P1/Pc),exp1),exp2);
      }
      else{
        T1 = T_char * pow(P1/P_char,Lad);
      }
      den1 = (P1/(R*T1));
    do{
          if(Pc>P2){
            T2 = Td * pow(1.0e0 + (Lad / (Lin - Lad))*pow((P2/Pc),exp1),exp2);
          }
          else{
            T2 = T_char * pow(P2/P_char,Lad);
          }
        
      
        
      
        den2 = (P2/(R*T2));
        dz = abs ((-P1+P2)* (1.0/grav)/((1.0/2.0)*(den1+den2)));
        error = dz-z_top/number_cell;//*float(count_cell);
        count++;
      if(count % 10000000 == 0){
        if(Pc>P1){
        cout << setprecision(15)<<"In atmos. (radiation)" <<setprecision(20) <<' ' <<count_cell<< ' ' << abs(abs((-P1+P2)/dz)+(grav)*((1.0e0/2.0e0)*(den1+den2)))<< ' ' << error << ' ' <<dz/(z_top/number_cell)<<' '<<P2<< '\n';
        }
        else{
        cout << setprecision(15)<<"In atmos.(convection)" <<setprecision(20) <<' ' <<count_cell<< ' ' << abs(abs((-P1+P2)/dz)+(grav)*((1.0e0/2.0e0)*(den1+den2)))<< ' ' << error << ' ' <<dz/(z_top/number_cell)<<' '<<P2<< '\n';
        }
          
      }
      if(abs(error) > z_top/number_cell){
        dP = dP*0.5e0;
        if(error > 0.0e0){
          P2 = P2 - dP;
        }
        else{
          P2 = P2 + dP;
        }
      cout<<"Too large dP"<<' '<< dP<< ' ' <<P2<<' '<<dz<<' '<< error << ' ' <<z_top/number_cell<<'\n';
      }
      else if(abs(error) > error_expected){
            if(error > 0.0e0){
              
              if((P2 + P2_temp)/2.0e0 < P1){
                do{
                  count3++;
                  if(count3 == int(z0)){
                    count3 = 0;
                    z0 = z0*10.0e0;
                  }
                }while((P2 + P2_temp)/(2.0e0 - float(count3)/z0) <= P1);
                P2 = (P2 + P2_temp)/(2.0e0 - float(count3)/z0);
                count3 = 0;
              }
              else{
                P2 = (P2 + P2_temp)/2.0e0;
              }
            }
            else if(error<0){
              P2_temp = max(P2,0.0e0);
              P2 = P2 + P2/dP_factor*dP_factor2*pow(abs(error/error_expected),0.16e0);
            }
          
          }
      
          else if(abs(error)<error_expected){
  cout<<"RECORD               "<< setprecision(20) <<' ' <<count_cell<< ' ' << abs(abs((-P1+P2)/dz)+(grav)*((1.0e0/2.0e0)*(den1+den2)))<< ' ' << error << ' ' <<dz/(z_top/number_cell)<<' '<<P2<< '\n';
          z = z - z_top/number_cell;
            if(z < 0e0){
              break;
            }
          T[count2] = T2;
          P[count2] = P2;
          den[count2] = den2;
          y[count2] = z;
          count2++;
          count++;
          P1 = P2;
          count_cell++;
          }
      if (dP<1.e-20){
        return 0;
      }
    }while(P1!=P2);
    }
  }while(z>=0);
  
    
    
  data_output << "# npts = " << count2 <<'\n' << "# num of variables = 4" << '\n' << "# density" << '\n' << "# temperature" <<'\n' << "# pressure" << '\n' << "# X " << '\n';
 
  do{
    --count2;
    data_output << setprecision(15) << y[count2] << ' ' << den[count2] << ' ' << T[count2] << ' ' << P[count2] << ' '<<  X << '\n';
  }while(count2>0);
    
    
  cout<<"count= " << count2<<'\n';
  cout<<"Pc[dyne/cm^{2}]= "<< ' ' << Pc << '\n';
  cout<<"Tc[K]= "<< ' ' << Tc << '\n';
  data_output.close();
  plot_output.close();

  return 0;
  
  
  
}

  
  
  
  


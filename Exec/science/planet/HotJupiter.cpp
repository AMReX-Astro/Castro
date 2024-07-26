//Aug 25

#include<iostream>
#include<cmath>
#include<fstream>
#include<iomanip>

using namespace std;
long double difference(long double dz,long double delta_z);    //declare the function for the given equation
long double difference(long double dz,long double delta_z)    //define the function here, ie give the equation
{
  long double difference=dz-delta_z;
  return difference;
}

long double Temp(long double P1,long double Pc,long double P_char,long double T_char,long double Td,  long double Lin, long double Lad, long double exp1, long double exp2);    //declare the function for the given equation
long double Temp(long double P1,long double Pc,long double P_char,long double T_char,long double Td,  long double Lin, long double Lad, long double exp1, long double exp2)    //define the function here, ie give the equation
{
  long double Temp;

  if(Pc>P1){
    Temp = Td * std::pow(1.0e0 + (Lad / (Lin - Lad)) * std::pow((P1/Pc),exp1),exp2);
  }
  else{
    Temp = T_char * std::pow(P1/P_char,Lad);
  }
  return Temp;
}

long double density(long double P1,long double T1,long double R);    //declare the function for the given equation
long double density(long double P1,long double T1,long double R)    //define the function here, ie give the equation
{
  long double density=P1/(R*T1);
  return density;
}


int main(){
  ofstream data_output,plot_output;
  long double P1, dP,P_top,den_top,z_top;
  long double R , z , den1 ,  T1,  delta_z;
  long double dP_factor,P_at_1opticaldepth,den_at_1opticaldepth,T_at_1opticaldepth,z_at_1opticaldepth;
  double X = 1.0,buffer_height,T_buffer,den_buffer,P_buffer;
  long double P2, P3,dz,dz2,dz3, den2,den3, T2,T3,error_expected,error;
  long double P_sol,T_sol,den_sol,dz_sol,kappa0;
  long double Pc,Pd, Lad, Lin,P_char,T_char,Tc,grav ,Td,exp1, exp2,number_cell;
  double P2_temp,optical_depth_buffer,optical_depth_atm,optical_depth_radiative;
  int alpha=1, beta=0, count = 0,count_cell,count2=0;
  bool continuous_P, continuous_rho;
  long double P[10000],T[10000],den[10000],y[10000];

  //Declaration Block
  Td = 1500.0;                              //T_deep : the temperature at the top of the atmosphere[K]
  P_char = 1.0e6;
  T_char = 250.0;
  dP = 1.0e+3;                              //dP default value
  Lad = (2.0/7.0);                          //the adiabatic index of an ideal diatomic gas
  Lin = (1.0/2.0);
  grav = -1.0e3;                            //constant gravity [cm/s^2]
  R = (1.3806488e-16)/2.34*6.02214129e23;   //ideal gas constant (defined with mean molecular mass)

  exp1=1.0 + alpha;
  exp2=(1.0/(4.0-beta));


  //Defining the computation box
  z_top=4.0e9;                               // the size of the computation box
  number_cell=1024.0;                        // the number of cells


  //Initializing
  kappa0 = 6.35e-3 ;
  dP_factor=1.e5;
  optical_depth_atm=0.0;
  optical_depth_radiative=0.0;
  delta_z=z_top/number_cell;
  z=z_top+1.5*delta_z;
  den_top=1.e-20;
  count=0;
  error_expected=1.e-10/z_top;
  T_at_1opticaldepth=0.0;
  P_at_1opticaldepth=0.0;
  den_at_1opticaldepth=0.0;
  z_at_1opticaldepth=0.0;
  data_output.open ("newmodelcpp.hse",ios::out);

  //Choosing either continuous pressure or density at the boundary between the buffer and the atmosphere
  continuous_P = true;
  continuous_rho = false;

  if((continuous_P==true && continuous_rho==true) || (continuous_P==false && continuous_rho==false)){
    cout<<"P and rho can not be continuous at the same time." <<'\n';
  }

  //Defining the properties of the planet atmosphere
                            // the density limit at the top of the atmosphere
  P_top=den_top*R*Td ;                       // the pressure limit at the density limit
  Tc = Td * std::pow(Lin/(Lin-Lad),exp2);                 //the temperature at the radiative-convective boundary (RCB)
  Pc = P_char * std::pow(Tc/T_char,1.0/Lad);              //the pressure at the RCB
  Pd = std::pow((Lin-Lad)/(2.0*Lin-Lad),1.0/exp1)*Pc;     //the characteristic pressure defining the radiative zone



  //Defining the properties of the buffer
  den_buffer=1.3e-15;                          // the density in the buffer above the atmosphere
  T_buffer=1.0e-2;                             // the temperature in the buffer above the atmosphere (not continuous at RCB). -> needs sponge
  P_buffer=den_buffer*R*T_buffer;
  if(continuous_rho==true){
    den_top=den_buffer;
  }
    else if(continuous_P==true){
  }

  buffer_height=z_top*2.25/4.0;               // the height at which the buffer and the top of the atmosphere meet (not continuous at RCB).-> needs sponge
  optical_depth_buffer=6.35e-3*T_buffer*std::pow(den_buffer,2.0)*(z_top-buffer_height);// the optical depth in the buffer region [dimensionless]



  do{
    if(z>buffer_height){
        den1 = den_buffer;
        den2 = den_buffer;
        T1   = T_buffer;
        T2   = T_buffer;
        P1   = P_buffer;
        P2   = P_buffer;
        T[count] = T1;
        P[count] = P1;
        den[count] = den1;
        y[count] = z;
        count_cell++;
        cout<<"In buffer" << ' '<< count_cell << ' ' << z<<' '<< dz/(z_top/number_cell) << ' '<<"Pressure[P/Pc]="<<' ' << P1/Pc << ' '<<"optical depth ="<< ' '<<optical_depth_buffer<<'\n';

        z = z-delta_z;
        count++;

        if(continuous_rho==true){
          P1 = R*Td*den_buffer;
          P2 = R*Td*den_buffer;
        }
      }
    else{

      if(Pc/10.0<P1){
        dP_factor = 1.e5;
      }

      T1=Temp(P1,Pc,P_char,T_char, Td,Lin,Lad, exp1,exp2);
      den1 = density(P1,T1,R);

      dP = P1/dP_factor;
      P2 = P1;
      P3 = P1;
      P2 = P2+dP;

      T2 = Temp(P2, Pc, P_char, T_char, Td, Lin, Lad, exp1, exp2);
      den2 = density(P2, T2, R);
      T3 = T1;
      den3 = den1;
      dz2 = std::abs ((-P1+P2)* (1.0/grav)/((1.0/2.0)*(den1+den2)));
      dz3 = std::abs ((-P1+P3)* (1.0/grav)/((1.0/2.0)*(den1+den3)));

    a:  if(difference(dz2, delta_z)*difference(dz3, delta_z)<0.0){
        P_sol = (P2+P3)/2.0;
        T_sol = Temp(P_sol, Pc, P_char, T_char, Td, Lin, Lad, exp1, exp2);
        den_sol = density(P_sol,T_sol, R);
        dz_sol = std::abs ((-P1+P_sol)* (1.0/grav)/((1.0/2.0)*(den1+den_sol)));
        count2++;
        }
        else{
          if(difference(dz2, delta_z)>0.0){
          P3 = P3 - P3/dP_factor/1.0;
          T3 = Temp(P3, Pc, P_char, T_char, Td, Lin, Lad, exp1, exp2);
          den3 = density(P3, T3, R);
          dz3 = std::abs ((-P1+P3)* (1.0/grav)/((1.0/2.0)*(den1+den3)));
          }
        else if(difference(dz2, delta_z)<0.0){
          P2 = P2 + P2/dP_factor/1.0;
          T2 = Temp(P2, Pc, P_char, T_char, Td, Lin, Lad, exp1, exp2);
          den2 = density(P2,T2,R);
          dz2 = std::abs ((-P1+P2)* (1.0/grav)/((1.0/2.0)*(den1+den2)));
          }
          goto a;
        }

      if(difference(dz_sol, delta_z)*difference(dz2, delta_z)<0.0){
        P3 = P_sol;
      }
      else{
        P2 = P_sol;
      }
      if(count2>100000){
        goto b;
      }
      if(std::abs(difference(dz_sol, delta_z)) >= error_expected){
        goto a;
      }

    b:  cout<<count_cell<<' ' <<setprecision(5) << "Final value, "<<' '<<"height[km]="<< ' '<< z/1.e5<<' ' <<"pressure[bar]=" <<' ' << P_sol/1.e6  <<' '<< "error="<<' ' << std::abs(difference(dz_sol, delta_z))<<' ' << std::abs ((-P1+P_sol)/dz_sol+ grav*((1.0/2.0)*(den1+den_sol))) << '\n';
      T[count] = T_sol;
      P[count] = P_sol;
      den[count] = den_sol;
      y[count] = z;
      count++;
      count2 = 0;
      z = z - delta_z;
      P1 = P_sol;
      count_cell++;
      if(P_sol < Pc){
      optical_depth_radiative = optical_depth_radiative+kappa0*std::pow(den_sol,2.0)*T_sol*delta_z;
      }
      optical_depth_atm = optical_depth_atm+kappa0*std::pow(den_sol,2.0)*T_sol*delta_z;
      if(optical_depth_atm+optical_depth_buffer>=1.0 && P_at_1opticaldepth==0.0 ){
        P_at_1opticaldepth=P_sol;
        T_at_1opticaldepth=T_sol;
        den_at_1opticaldepth=den_sol;
        z_at_1opticaldepth=y[count-1];
      }

    }



  }while(z > -2.0*delta_z);



  //output data
  cout<<"count= " << count-1<<'\n';
  cout<<"At bottom, "<< ' ' <<"density[g/cm^3]= " << ' ' << den[count-1] << ' ' << "   T[K]= "<< ' ' << T[count-1] << ' ' << "   P[bar]= " << ' ' << P[count-1]/1.e6  << '\n';
  cout<<"At optical depta=1, "<< ' ' <<"z[km]=" << ' ' << z_at_1opticaldepth/1e5 << ' ' <<"density[g/cm^3]= " << ' ' << den_at_1opticaldepth << ' ' << "   T[K]= "<< ' ' << T_at_1opticaldepth << ' ' << "   P[bar]= " << ' ' << P_at_1opticaldepth/1.e6  << '\n';

  if (optical_depth_buffer > 1.0){
    cout<<"Optically [thick] buffer : " << ' ' << "optical depth = " << ' '<<optical_depth_buffer << '\n';
  }
  else{
    cout<<"Optically [thin ] buffer : " << ' ' << "optical depth = " << ' '<<optical_depth_buffer << '\n';

  }


  data_output << "# npts = " << count <<'\n' << "# num of variables = 4" << '\n' << "# density" << '\n' << "# temperature" <<'\n' << "# pressure" << '\n' << "# X " << '\n';

  do{
    --count;
    data_output << setprecision(25) << y[count] << ' ' << den[count] << ' ' << T[count] << ' ' << P[count] << ' '<<  X << '\n';
  }while(count>0);


  //summary
  cout<<"Pc                      [dyne/cm^{2}] = "<< ' ' << Pc << '\n';
  cout<<"Tc                                [K] = "<< ' ' << Tc << '\n';
  cout<<"density_c                    [g/cm^3] =" <<' ' << Pc/R/Tc << '\n';
  cout<<"optical depth        (radiative zone) =" <<' ' << optical_depth_radiative << '\n';
  cout<<"optical depth (radiative+ convective) =" <<' ' << optical_depth_atm << '\n';
  data_output.close();
  return 0;



}







#include <iostream>
#include <math.h>
#define pi acos(-1)
#include"string"
using namespace std;
 
// Abstract base class
class EquationSolver
{
   public:
   virtual double Solve()=0;
    
    void setfun(string f)
    {
      func=f;
      
    }
    
    void setfun(string f1,string f2)
    {
      func=f1;
      func1=f2;
    }
    void setnum(double l,double h,double D,double B1){
      cin>>l>>h>>D>>B1;
       A=l*sin(B1);
       B=l*cos(B1);
       C=(h+0.5D)*sin(B1)-0.5*D*tan(B1);
       E=(h+0.5*D)*cos(B1)-0.5*D;}
       
    double func_1(double x)
   {return (1/x)-tan(x);}
   
    double func_2(double x)
    {return 1/x-pow(2,x);}
    
    double func_3(double x)
    {return 1/(pow(2,x))+exp(x)+2*cos(x)-6;}
    
    double func_4(double x)
    {return (pow(x,3)+4*pow(x,2)+3*x+5)/(2*pow(x,3)-9*pow(x,2)+18*x-2);}
    
    double func_5(double x)
    {return x-tan(x);}
    
    double derivfunc_5(double x)
    {return 1-pow(1/cos(x),2);}
    
    double func_6(double x)
    {return sin(x/2)-1;}
    
    double func_7(double x)
    {return exp(x)-tan(x);}
    
    double func_8(double x)
    {return pow(x,3)-12*pow(x,2)+3*x+1;}
    
    double func_9(double x)
    {return 5*pi-10*asin(x)-10*x*sqrt(1-pow(x,2))-12.4;}
    
    double derivfunc_9(double x)
    {return (-10/sqrt(1-pow(x,2)))-10*sqrt(1-pow(x,2))+(10*pow(x,2))/sqrt(1-pow(x,2));}
    
    double func_10(double x)
    { return A*sin(x)*cos(x)+B*pow(sin(x),2)-C*cos(x)-E*sin(x);}
      
    double derivfunc_10(double x)
    {return A*pow(cos(x),2)-A*pow(sin(x),2)+2*B*sin(x)*cos(x)+C*sin(x)-E*cos(x);}
     
    double function(double x,string function)
    {
    if(function == "func_1") {
        return func_1(x);
  } else if(function == "func_2") {
        return func_2(x);} 
    else if(function == "func_5") {
        return func_5(x);}
    else if(function == "derivfunc_5") {
        return derivfunc_5(x);}
    else if(function == "func_6") {
        return func_6(x);}
    else if(function == "func_7") {
        return func_7(x);}
    else if(function == "func_8") {
        return func_8(x);}
    else if(function == "func_9") {
        return func_9(x);}
    else if(function == "derivfunc_9")
    {return derivfunc_9(x);}
    else if(function == "func_10") {
        return func_10(x);}
     else if(function == "derivfunc_10")
    {return derivfunc_10(x);}
        
 
    else if(function == "func_3"){
      return func_3(x);}
    else{
        return func_4(x);
    }
  }
   void seta(double a)
   {
      LeftInterval = a;
   }
   void setb(double b)
   {
      RightInterval = b;
   }
   protected:
   string func;
   string func1;
   double LeftInterval;
   double RightInterval;
   double A;
   double B;
   double C;
   double E;
   };
   
// Derived class
class BisectionMethod: public EquationSolver
{
   public:
   double Solve()
   {
     double u=function(LeftInterval,func);
     double v=function(RightInterval,func);
      if ((u*v)>=0)
      {
         cout << "There are no roots in this interval.";
         return 0;
      }
      
      double h=LeftInterval;
      while ((RightInterval-LeftInterval)>=0.0001)
      {
         h=(LeftInterval+RightInterval)/2;
         if (function(h,func)==0.0)
         break;

         else if (function(h,func)*function(LeftInterval,func)<0)
         RightInterval=h;

         else
         LeftInterval=h;
      }
      return h;
   }
};


class NewtonMethod: public EquationSolver
{
   public:
   double Solve()
   { int m=50;
     double x0=LeftInterval;
     
     for(int i=0;i<=m;i++) 
     {
       double u=function(x0,func);
       if(abs(u)<0.001)
       {break;}
       else{
         double h=function(x0,func)/function(x0,func1);
         double x1=x0-h;
         x0=x1;
       }
     }
     return x0;
    }
};


class SecantMethod: public EquationSolver
{
   public:
   double Solve()
   
   { double x_1=RightInterval;
     double x_0=LeftInterval;
     double u=function(x_1,func);
     double v=function(x_0,func);
     for(int i=2;i<=50;i++)
     {
       if (abs(u)>abs(v))
       {
         swap(x_1,x_0);
         swap(u,v);
        }
        double s=((x_1)-(x_0))/(u-v);
        x_0=x_1;
        v=u;
        x_1=x_1-u*s;
        u=function(x_1,func);
        if(abs((x_1)-(x_0))<0.00001||abs(u)<0.00001)
        {break;}
     }
    return x_1;
    }
    
};

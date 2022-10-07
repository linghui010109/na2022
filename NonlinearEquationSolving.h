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


int main(void)
{
  //Question B(i)
   BisectionMethod Bis1;
   Bis1.seta(0);
   Bis1.setb(pi/2);
   string a="func_1";
   Bis1.setfun(a);
   
   cout<<"B(i).The value of root is:"<<Bis1.Solve()<<endl;
   
   
   //Question B(ii)
   BisectionMethod Bis2;
   Bis2.seta(0);
   Bis2.setb(1);
   string b="func_2";
   Bis2.setfun(b);
   
   cout<<"B(ii).The value of root is:"<<Bis2.Solve()<<endl;
   
   
   //Question B(iii)
   BisectionMethod Bis3;
   Bis3.seta(1);
   Bis3.setb(3);
   string c="func_3";
   Bis3.setfun(c);
   
   cout<<"B(iii).The value of root is:"<<Bis3.Solve()<<endl;
   
   
   //Question B(iv)
   BisectionMethod Bis4;
   Bis4.seta(0);
   Bis4.setb(4);
   string d="func_4";
   Bis4.setfun(d);
   
   cout<<"B(iv).The value of root is:"<<Bis4.Solve()<<endl;
   
   //Question C
   NewtonMethod New1;
   New1.seta(4.5);
   string e="func_5";
   string f="derivfunc_5";
   New1.setfun(e,f);
   
    cout<<"C(i).The value of root is:"<<New1.Solve()<<endl;
    
   NewtonMethod New2;
   New2.seta(7.7);
   string w="func_5";
   string u="derivfunc_5";
   New2.setfun(w,u);
   
    cout<<"C(ii).The value of root is:"<<New2.Solve()<<endl;
   
   //Question D(i)
   SecantMethod Sec1;
   Sec1.seta(0);
   Sec1.setb(pi/2);
   string g="func_6";
   Sec1.setfun(g);
   
   cout<<"D(i).The value of root is:"<<Sec1.Solve()<<endl;
   
   //Question D(ii)
   SecantMethod Sec2;
   Sec2.seta(1);
   Sec2.setb(1.4);
   string h="func_7";
   Sec2.setfun(h);
   
   cout<<"D(ii).The value of root is:"<<Sec2.Solve()<<endl;
   
   //Question D(iii)
   SecantMethod Sec3;
   Sec3.seta(0);
   Sec3.setb(-0.5);
   string i="func_8";
   Sec3.setfun(i);
   
   cout<<"D(iii).The value of root is:"<<Sec3.Solve()<<endl;
   
   //Question E
   //Test1
   BisectionMethod Bis5;
   Bis5.seta(0);
   Bis5.setb(1);
   string j="func_9";
   Bis5.setfun(j);
   
   cout<<"E.The depth of water(BisectionMethod) ="<<1-Bis5.Solve()<<endl;
   
   //Test2
   NewtonMethod New3;
   New3.seta(0);
   string k="func_9";
   string l="derivfunc_9";
   New3.setfun(k,l);
   
    cout<<"E.The depth of water(NewtonMethod) ="<<1-New3.Solve()<<endl;
    
    //Test3
   SecantMethod Sec4;
   Sec4.seta(0);
   Sec4.setb(1);
   string m="func_9";
   Sec4.setfun(m);
   
   cout<<"E.The depth of water(SecantMethod) ="<<1-Sec4.Solve()<<endl;
   
    //Question F(a)
   NewtonMethod New4;
   New4.setnum(80,49,55,11.5);
   New4.seta(30);
   string o="func_10";
   string n="derivfunc_10";
   New4.setfun(o,n);
   
    cout<<"F(a).The value a is:"<<New4.Solve()<<endl;
	
	//Question F(b)
   NewtonMethod New5;
   New5.setnum(80,49,30,11.5);
   New5.seta(33);
   string p="func_10";
   string q="derivfunc_10";
   New5.setfun(p,q);

   cout<<"F(b).The value a is:"<<New5.Solve()<<endl;
	
   //Question F(c)
   SecantMethod Sec5;
   Sec5.setnum(80,49,55,11.5);
   Sec5.seta(0);
   Sec5.setb(50);
   string r="func_10";
   Sec5.setfun(r);
   cout<<"F(c).The value a is:"<<Sec5.Solve()<<endl;
   
   return 0;
   
  
}
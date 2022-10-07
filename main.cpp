#include <iostream>
#include <math.h>
#define pi acos(-1)
#include"string"
#include"NonlinearEquationSolving.h"
using namespace std;

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

#include <iostream>
#include <iostream>
#include <bitset>
#include <cmath>
using namespace std;

int main()
{
    int beta=2,p=3,L=-1,U=1;
    double UFL=pow(beta,L);
    double OFL=pow(beta,U)*(beta-pow(beta,1-p));
    cout << "UFL = "<<UFL<<endl;
    cout << "OFL = "<<OFL<<endl;

    int count=1;
    cout<<"Normal number in F:"<<endl;
    for (int i=L;i<=U;i++)
    {
        for (int j=0;j<=1;j++)
        {
            for (int k=0;k<=1;k++)
            {
                
                cout<<"1."<<j<<k<<"*2^"<<i<<endl;
                cout<<"-1."<<j<<k<<"*2^"<<i<<endl;
                count+=2;
            }
        }
    }
    cout<<"#F = "<<count<<endl;
    
    cout<<"Subnormal number in F:"<<endl;
    for (int j=0;j<=1;j++)
    {     
        for (int k=0;k<=1;k++)
        {
            if (j==0 & k==0)
            {continue;}
            else
            {
                cout<<"0."<<j<<k<<"*2^-1"<<endl;
                cout<<"-0."<<j<<k<<"*2^-1"<<endl;
            }
        }
    }
    
}

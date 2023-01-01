#include <iostream>
#include <cmath>
#include <vector>
#include <cstring>
#include <fstream>
#include "Eigen/Dense"
#include "Eigen/Core"
#include "Splines.h"

#define pi acos(-1)
using namespace std;
using namespace Eigen;

const double _epsL = 10 * std::numeric_limits<double>::epsilon();
/*class f_1 :public Function
{
public:
	double operator()(const double &x)
	const{
		return log(x)  ;
	}
};*/
class f_2 :public Function
{
public:
	double operator()(const double& x)
		const {
		return 1/(1+25*x*x);
	}
};
/*void setpoint(vector<double> x, int N, double L, double R)
{
	for (int i = 0; i < N; i++)
	{	
		x.push_back(L + (R - L) * i / (N - 1));
	}
}*/


int main(void)
{
	/*vector<double> x1(5), y1(5), val(5);
	x1 = { 1.0,2.0,3.0,4.0,6.0 };
	val = { 1.1,3.4,2.6,4.5,5.7,6.1};
	f_1 f1;
	for (int i = 0; i < 5; i++)
	{
		//y1[i] = f1(x1[i]);
	}
	//PPForm pp1(val,x1,y1,f1,"complete");
	*/
	//vector<double> x2(6), y2(6),val2(6);
	
	//f_2 f2;
	//for (int i = 0; i < 6; i++)
	//{
		//x2[i] = -1.0 + (2.0 * i / 5.0);
		//y2[i] = f2(x2[i]);
		//val2[i] = -1.0 + (2.0 * i / 5.0);
		//cout << x2[i] << " ";
		//cout << y2[i] << endl;
	//}
	//val2 = { 0.01,-0.7 };
	//PPForm pp2(val2,x2, y2, f2, "complete");

	/*vector<double> x3(3), y3(3);
	x3 = { -1,0,1 };
	y3 = { 0,1,0 };
	PPForm pp2(x3, y3, f2, "natural");
  vector<double> x2, y2,val2;*/
  
	int n;
	f_2 f2;
	cout << "Please input the number of N:"<<endl;
	cin >> n;
	vector<double> x2(n), y2(n);
	for (int i = 0; i < n; i++)
	{
		x2[i] = -1.0 + (2.0 * i / (n-1));
		y2[i] = f2(x2[i]);
		//cout << x2[i] << endl;
	}
	int max;
	cout << "Please input the number of max sample:" << endl;
	cin >> max;
	vector<double> val2(max);
	for (int i = 0; i < max; i++)
	{
	val2[i] = -1.0 + (2.0 * i / (max - 1));
	//cout << val2[i] << endl;
	}
	
	PPForm pp1(val2, x2, y2, f2, "not a knot");

	return 0;
}

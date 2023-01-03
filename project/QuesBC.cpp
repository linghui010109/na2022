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
class f_1 :public Function
{
public:
	double operator()(const double& x)
		const {
		return 1/(1+pow(x,2));
	}
};

int main()
{
	vector <double> t1(11);
	for (int i = 0.0; i < 11; i++)
	{
		t1[i] = -6.0 + i+1.0;
		//cout << t1[i] << endl;
	}
	f_1 f1;
	int max;
	cout << "Please input the number of max sample:" << endl;
	cin >> max;
	vector<double> val(max);
	for (int i = 0; i < max; i++)
	{
		val[i] = -5.0 + (10.0 * i / (max - 1));
	}
	//vector<double> val = { -5,-4,-3,-2,-1,0,1.0,2.0,3.0,4.0,5.0};
	cout << "complete cubic cardinal:" << endl;
	B_Spline_c B1(t1, f1, val, "complete cubic cardinal");
	
	vector <double> t2(10);
	for (int i = 0.0; i < 10; i++)
	{
		t2[i] = i+1.0 - (11.0 / 2.0);
		//cout << t2[i] << endl;
	}
	cout << "complete quadratic cardinal:" << endl;
	B_Spline_c B2(t2, f1, val, "complete quadratic cardinal");
	return 0;
}

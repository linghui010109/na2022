#pragma once
#include <iostream>
#include <cmath>
#include <vector>
#include <limits>
#include "Eigen/Dense"
#include "Eigen/Core"
#define pi acos(-1)
using namespace std;
using namespace Eigen;

/*class SetPoint
{
public:
    MatrixXf m;
    int r;
    vector <double> mu, lambda;
    vector <double> x, y;
    void set_point(vector<double> xx, vector<double> yy)
    {
        x = xx;
        y = yy;
    }
};*/

class SolveMatrix
{
public:
    MatrixXf x;
    MatrixXf solve_matrix(MatrixXf A, MatrixXf b)
    {
        int s = A.size();
        x = A.inverse()*b;
        return x;
    }
};

class polynomail
{
public:
    polynomail(vector<double> &x,vector<double> &y,MatrixXf &M)
    {
        int N = size(M);
        vector<double> a(N - 1), b(N - 1), c(N - 1), d(N - 1);
        for (int i = 0; i < N - 1; i++)
        {
            a[i] = y[i];
            b[i] = ((y[i + 1] - y[i]) / (x[i + 1] - x[i])) - (M(i + 1, 0) + 2 * M(i, 0))*(x[i + 1] - x[i]) / 6;
            c[i] =M(i,0) / 2;
            d[i] = (M(i + 1, 0) - M(i, 0)) / ((x[i + 1] - x[i])*6);
        }
        for (int i = 0; i < N-1; i++)
        {
            cout << "g_" << i << "(X)=" << a[i] << "+" << b[i] << "(x-" << x[i] << ")" << "+" << c[i] << "(x-" << x[i] << ")^2+" << d[i] << "(x-" << x[i] << ")^3" << endl;
        }
    }
};
class getcoeff
{
public:
    getcoeff(vector<double> &value,MatrixXf &M, vector<double>& x, vector<double>& y)
    {
        int N = size(M);
        int num = size(value);
        vector<double> a(N - 1), b(N - 1), c(N - 1), d(N - 1);
        double ans;
        ofstream output;
        string filename;
        cout << "Write the filename" << endl;
        cin >> filename;
        output.open(filename);
        for (int i = 0; i < N - 1; i++)
        {
            a[i] = y[i];
            b[i] = ((y[i + 1] - y[i]) / (x[i + 1] - x[i])) - (M(i + 1, 0) + 2 * M(i, 0)) * (x[i + 1] - x[i]) / 6;
            c[i] = M(i, 0) / 2;
            d[i] = (M(i + 1, 0) - M(i, 0)) / ((x[i + 1] - x[i]) * 6);
        
            for (int j = 0; j < num; j++)
            {
                if (x[i] <= value[j] && value[j] < x[i + 1])
                {
                    ans = a[i] + b[i] * (value[j] - x[i]) + c[i] * std::pow((value[j] - x[i]),2) + d[i] * std::pow((value[j] - x[i]), 3);
                    //cout << value[j] << " " << ans <<  endl;
                    output << value[j]<<" "<<ans << endl;
                }
                else
                {
                    continue;
                }
            }
            
        }
    //output.close();
    }
};

class Function
{
    const double _epsL = 10 * numeric_limits<double>::epsilon();
public:
    virtual double operator () (const double& x) const = 0;
    virtual double diff(const double& x) const 
    {
        return ((*this)(x + _epsL) - (*this)(x - _epsL)) / (2*_epsL);
    }
    virtual double diff2(const double& x) const 
    {
        return ((*this)(x + 2 * _epsL) + (*this)(x - 2 * _epsL) - 2 * (*this)(x)) / (4 * _epsL);
    }
};


class PPForm
{
public:
    PPForm(vector<double> &val,vector<double> &x,vector<double> &y,Function &func,const string& type_name)
    {
        int N = x.size();
        //cout << N << endl;
        vector <double> h(N - 1);
        vector <double> mu(N+1);
        vector <double> lambda(N+1 );
        MatrixXf A(N, N);
        MatrixXf b(N, 1);
        MatrixXf M(N, 1);
        SolveMatrix s;
        if (type_name == "complete")
        {
            for (int i = 0; i < N - 1; i++)
            {
                h[i] = (y[i + 1] - y[i]) / (x[i + 1] - x[i]);
                //cout << h[i] << endl;
            }
            for (int i = 1; i < N - 1; i++)
            {
                mu[i + 1] = (x[i] - x[i - 1]) / (x[i + 1] - x[i - 1]);
                lambda[i + 1] = (x[i + 1] - x[i]) / (x[i + 1] - x[i - 1]);
            }
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    if (i == j + 1)
                    {
                        A(i, j) = mu[j + 2];
                    }
                    else if (i == j - 1)
                    {
                        A(i, j) = lambda[j];
                    }
                    else if (i == j)
                    {
                        A(i, j) = 2;
                    }
                    else
                    {
                        A(i,j)=0;
                    }
                }
            }
            A(0, 0) = A(N - 1, N - 1) = 2;
            A(0, 1) = A(N - 1, N - 2) = 1;

            for (int i = 1; i < N - 1; i++)
            {
                b(i, 0) = 6*(h[i] - h[i - 1]) / (x[i + 1] - x[i - 1]);
            }
            b(0, 0) = 6*(h[0] - func.diff(x[0])) / (x[1] - x[0]);
            //cout << func.diff(x[0]) << endl;
            //cout << func.diff(x[N-1]) << endl;
            b(N - 1, 0) = 6*(func.diff(x[N - 1]) - h[N - 2]) / (x[N - 1] - x[N - 2]);
            //cout << A <<endl;
            //cout << b <<endl;
            M=s.solve_matrix(A, b);
            //cout << M << endl;
            //polynomail po1(x, y, M);
            getcoeff get(val, M, x, y);
            
        }
        else if (type_name == "not a knot")
        {
            for (int i = 0; i < N - 1; i++)
            {
                h[i] = (y[i + 1] - y[i]) / (x[i + 1] - x[i]);
                //cout << h[i] << endl;
            }
            for (int i = 1; i < N - 1; i++)
            {
                mu[i + 1] = (x[i] - x[i - 1]) / (x[i + 1] - x[i - 1]);
                lambda[i + 1] = (x[i + 1] - x[i]) / (x[i + 1] - x[i - 1]);
            }
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    if (i == j + 1)
                    {
                        A(i, j) = mu[j + 2];
                    }
                    else if (i == j - 1)
                    {
                        A(i, j) = lambda[j];
                    }
                    else if (i == j)
                    {
                        A(i, j) = 2;
                    }
                    else
                    {
                        A(i, j) = 0;
                    }

                }
            }
            A(0, 0) = -lambda[2];
            A(0, 1) = A(N - 1, N - 2) = 1;
            A(0, 2) = -mu[2];
            A(N - 1, N - 3) = -lambda[N - 1];
            A(N - 1, N - 1) = -mu[N - 1];

            for (int i = 1; i < N - 1; i++)
            {
                b(i, 0) = 6 * (h[i] - h[i - 1]) / (x[i + 1] - x[i - 1]);
            }
            b(0, 0) =0;
            b(N - 1, 0) = 0;
            //cout << A << endl;
            //cout << b << endl;
            M = s.solve_matrix(A, b);
            //cout << M << endl;
            getcoeff get(val, M, x, y);
        }
        if (type_name == "natural")
        {
            for (int i = 0; i < N - 1; i++)
            {
                h[i] = (y[i + 1] - y[i]) / (x[i + 1] - x[i]);
                //cout << h[i] << endl;
            }
            for (int i = 1; i < N - 1; i++)
            {
                mu[i + 1] = (x[i] - x[i - 1]) / (x[i + 1] - x[i - 1]);
                lambda[i + 1] = (x[i + 1] - x[i]) / (x[i + 1] - x[i - 1]);
            }
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    if (i == j + 1)
                    {
                        A(i, j) = mu[j + 2];
                    }
                    else if (i == j - 1)
                    {
                        A(i, j) = lambda[j];
                    }
                    else if (i == j)
                    {
                        A(i, j) = 2;
                    }
                    else
                    {
                        A(i, j) = 0;
                    }

                }
            }
            A(0, 0) = A(N - 1, N - 1) = 1;
            A(0, 1) = A(N - 1, N - 2) = 0;

            for (int i = 1; i < N - 1; i++)
            {
                b(i, 0) = 6 * (h[i] - h[i - 1]) / (x[i + 1] - x[i - 1]);
            }
            b(0, 0) = 0;
            b(N - 1, 0) = 0;
            //cout << A << endl;
            //cout << b << endl;
            M = s.solve_matrix(A, b);
            //cout << M << endl;
            getcoeff get(val, M, x, y);
        }
        
    };
    ~PPForm() {};
};






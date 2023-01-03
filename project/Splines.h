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

class B_Spline_c
{
public:
    B_Spline_c(vector<double> &x,Function &func,vector<double> val, const string& type_name)
    {
        int N = x.size();
        MatrixXf A(N, N);
        MatrixXf b(N, 1);
        MatrixXf a(N, 1);
        SolveMatrix s;
        if (type_name == "complete cubic cardinal")
        {
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    if (i == j + 1)
                    {
                        A(i, j) = 1;
                    }
                    else if (i == j - 1)
                    {
                        A(i, j) = 1;
                    }
                    else if (i == j)
                    {
                        A(i, j) = 4;
                    }
                    else
                    {
                        A(i, j) = 0;
                    }

                }
            }
            A(0, 0) = A(N - 1, N - 1) = 2;
            A(0, 1) = A(N - 1, N - 2) = 1;
            for (int i = 1; i < N - 1; i++)
            {
                b(i, 0) = 6.0 * func(x[i]);
            }
            b(0, 0) = 3.0 * func(x[0]) + func.diff(x[0]);
            b(N - 1, 0) = 3.0 * func(x[N-1]) - func.diff(x[N-1]);
            //cout << A << endl;
            //cout << A.rows() << endl;
            //cout << b << endl;
            //cout << b.size() << endl;
            a = s.solve_matrix(A, b);
            //cout << a << endl;
            double a_1 = a(1, 0) - 2 * func.diff(x[0]);
            double a_N = a(N - 2, 0) + 2 * func.diff(x[N - 1]);
            ofstream output;
            string filename;
            cout << "Write the filename" << endl;
            cin >> filename;
            output.open(filename);
            for (int j = 0; j < val.size(); j++)
            {
                if (val[j] < 0)
                {
                    output << val[j]; 
                    //cout << val[j];
                    val[j] = -val[j];
                }
                else
                {
                    output << val[j];
                    //cout << val[j];
                }
                double final_val = a_1 * B(-1, 3, val[j]) + a_N * B(N, 3, val[j]);
                //cout << B(-1, 3, val[j]) << endl;
                for (int i = 0; i < N; i++)
                {
                    double ans = a(i, 0) * B(i, 3, val[j]);
                    final_val = ans + final_val;
                }
                output<< " " << final_val << endl;
                //cout <<" "<< final_val << endl;
            }
            output.close();
            
        }
        else if (type_name == "complete quadratic cardinal")
        {
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    if (i == j + 1)
                    {
                        A(i, j) = 1.0;
                    }
                    else if (i == j - 1)
                    {
                        A(i, j) = 1.0;
                    }
                    else if (i == j)
                    {
                        A(i, j) = 6.0;
                    }
                    else
                    {
                        A(i, j) = 0.0;
                    }

                }
            }
            A(0, 0) = A(N - 1, N - 1) = 5.0;
            A(0, 1) = A(N - 1, N - 2) = 1.0;
            for (int i = 1; i < N-1 ; i++)
            {
                b(i, 0) = 8.0 * func(x[i]);
            }
            b(0, 0) = 8.0 * func(x[0]) - 2.0 * func(1.0);
            b(N - 1, 0) = 8.0 * func(x[N-1]) - 2.0 * func(N+1);
            a = s.solve_matrix(A, b);
            double a_0 = 2 * func(1.0) - a(0, 0);
            double a_N = 2 * func(N + 1) - a(N - 1, 0);
            //cout << A << endl;
            //cout << b << endl;
            //cout << a << endl;
            ofstream output;
            string filename;
            cout << "Write the filename" << endl;
            cin >> filename;
            output.open(filename);
            for (int j = 0; j < val.size(); j++)
            {
                if (val[j] < 0)
                {
                    output << val[j];
                    //cout << val[j];
                    val[j] = -val[j];
                }
                else
                {
                    output << val[j];
                    //cout << val[j];
                }
                double final_val = a_0*B(0,2,val[j]) + a_N * B(N + 1, 2, val[j]);
                for (int i = 0; i < N; i++)
                {

                    double ans = a(i, 0) * B(i + 1, 2, val[j]);
                    final_val = ans + final_val;
                }
                output << " " << final_val << endl;
               //cout <<" "<< final_val << endl;
            }
            output.close();
        }
        else
        {
            cout << "Wrong type name!" << endl;
        }
        
    }
    /*double B(int i, int n, double val, vector<double>& t)
    {
        vector<double>
        if (n == 0)
        {
            if (val > t[i + 1] && val <= t[i + 2])
            {
                return 1.0;
            }
            else
            {
                return 0;
            }
        }
        else
        {
            return (val - t[i + 1]) / (t[i + n + 1] - t[i + 1]) * B(n - 1, i, val,t) + (t[i + n + 2] - val) / (t[i + n + 2] - t[i + 2]) * B(n - 1, i + 1, val,t);
        }
    }*/
    double B(int i, int n, double val)
    {
        //int temp = n;
        if (n == 0)
            {
                if (val > i-1 && val <= i)
                {
                    return 1.0;
                }
                else
                {
                    return 0.0;
                }
            }
            else
            {
                //double ans;
                double x1 = B(i, n - 1, val);
                double x2 = B(i + 1, n - 1, val);
                return ((val - i+1) * x1) / n + ((i+n - val) * x2) / n;
            
            }
    }
    
};

class B_Spline_n
{
public:
    B_Spline_n(double val, vector<double>& t, Function& func, const string& type_name)
    {
        int N = t.size();
        SolveMatrix s1;
        MatrixXf a1(N, 1);
        MatrixXf A(N + 2, N + 2), b(N + 2, 1), a2(N + 2, 1);
        double dis = (t[1] - t[0]);
        double t_0 = t[0] - dis;
        double t_1 = t[0] - 2 * dis;
        double t_2 = t[0] - 3 * dis;
        double t_N_1 = t[N - 1] + dis;
        double t_N_2 = t[N - 1] + 2 * dis;

        if (type_name == "complete cubic cardinal")
        {
            for (int i = 1; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    if (i == j + 1)
                    {
                        A(i, j) = 1;
                    }
                    else if (i == j - 1)
                    {
                        A(i, j) = 1;
                    }
                    else if (i == j)
                    {
                        A(i, j) = 4;
                    }
                    else
                    {
                        A(i, j) = 0;
                    }

                }
            }
            A(0, 0) = A(N - 1, N - 1) = 2;
            A(0, 1) = A(N - 1, N - 2) = 1;
            for (int i = 1; i < N - 1; i++)
            {
                b(i, 0) = 6.0 * func(i + 1);
            }
            b(0, 0) = 3.0 * func(1) + func.diff(1);
            b(N - 1, 0) = 3.0 * func(N) + func.diff(N);
        }
        if (type_name == "linear")
        {
            for (int j = 0; j < N; j++)
            {
                a1(j, 0) = func(t[j]);
            }
        }

        if (type_name == "complete")
        {
            for (int i = 1; i < N + 1; i++)
            {
                for (int j = 0; j < N + 2; j++)
                {
                    if (i == j + 1)
                    {
                        A(i, j) = B(i - 2, 3, t[i - 1], t);
                    }
                    else if (i == j - 1)
                    {
                        A(i, j) = B(i, 3, t[i - 1], t);
                    }
                    else if (i == j)
                    {
                        A(i, j) = B(i - 1, 3, t[i - 1], t);
                    }
                    else
                    {
                        A(i, j) = 0;
                    }

                }
            }
            A(0, 0) = -3 * B(0, 2, t[0], t) / (t[1] - t_1);
            A(0, 1) = (3 * B(0, 2, t[0], t) / (t[1] - t_1)) - (3 * B(1, 2, t[0], t) / (t[2] - t_0));
            A(0, 2) = 3 * B(1, 2, t[0], t) / (t[1] - t_0);
            A(N - 1, N - 3) = -3 * B(N - 1, 2, t[N - 1], t) / (t_N_1 - t[N - 3]);
            A(N - 1, N - 2) = (3 * B(N - 1, 2, t[N - 1], t) / (t_N_1 - t[N - 3])) - (3 * B(N, 2, t[N - 1], t) / (t_N_2 - t[N - 2]));
            A(N - 1, N - 1) = 3 * B(N, 2, t[N - 1], t) / (t_N_2 - t[N - 2]);
            for (int i = 1; i < N + 1; i++)
            {
                b(i, 0) = func(t[i - 1]);
            }
            b(0, 0) = func.diff(t[0]);
            b(N + 1, 0) = func.diff(t[N - 1]);
            cout << A << endl;
            cout << b << endl;
            a2 = s1.solve_matrix(A, b);
        }
        if (type_name == "natural")
        {
            for (int i = 1; i < N + 1; i++)
            {
                for (int j = 0; j < N + 2; j++)
                {
                    if (i == j + 1)
                    {
                        A(i, j) = B(i - 2, 3, t[i - 1], t);
                    }
                    else if (i == j - 1)
                    {
                        A(i, j) = B(i, 3, t[i - 1], t);
                    }
                    else if (i == j)
                    {
                        A(i, j) = B(i - 1, 3, t[i - 1], t);
                    }
                    else
                    {
                        A(i, j) = 0;
                    }

                }
            }

        }

    }

    double B(int i, int n, double val, vector<double>& t)
    {
        int temp = n;
        if (temp == 0)
        {
            if (val > t[i + 1] && val <= t[i + 2])
            {
                return 1.0;
            }
            else
            {
                return 0;
            }
        }
        else
        {
            double ans;
            ans = ((val - t[i + 1]) * B(i, temp - 1, val, t)) / (t[i + temp + 1] - t[i + 1]) + ((t[i + temp + 2] - val) * B(i + 1, temp - 1, val, t)) / (t[i + temp + 2] - t[i + 2]);
            temp--;
            return ans;
        }
    }
    void get_result(MatrixXd& M, int number_of_N, double val, vector<double>& t)
    {
        int N = number_of_N;
        double xx = 0.0;
        for (int i = 0; i < N; i++)
        {
            double x = M(i, 0) * B(i, 3, val, t);
            xx = xx + x;
        }
    }
};





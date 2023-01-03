# Splines

### -架构

##### PPform
1.PPForm(vector<double> &val,vector<double> &x,vector<double> &y,Function &func,const string& type_name):val为输入值，x和y分别为节点的值，func为函数，type_name为函数的类别，分别有complete、not a knot、natural以及linear，此函数直接返回val的拟合函数值的文件。  
2.SolveMatrix（MatrixXf A, MatrixXf b）：解线性方程组。  
3.polynomail(vector<double> &x,vector<double> &y,MatrixXf &M):返回PPFrom的拟合函数。  
4.getcoeff(vector<double> &value,MatrixXf &M, vector<double>& x, vector<double>& y)：返回拟合函数的值。  
5.Function：用来管理函数以及关于函数一阶微分与二阶微分的定义。  

##### Bspline
1.B_Spline_c(vector<double> &x,Function &func,vector<double> val, const string& type_name):x为插值节点，func为函数，val为输入值，type_name为函数的类型，分别有complete cubic cardinal以及complete quadratic cardinal。  
2.B(int i, int n, double val)：cardinal的基函数。  
3.B_Spline_n（double val, vector<double>& t, Function& func, const string& type_name, const double ii）：val为输入值，t为插值节点，func为函数，type_name为函数类型，分别为linear，conplete以及natural。（但我来不及完成，只完成一半）  
  
 
### -测试过程
A.通过用户输入N，即均分点数以及输入拟合数据的大小和数据的文件名称来得到函数拟合值，函数拟合值将直接保存在目录下，拟合数据也是均分的。（我是通过在区间[-1,1]上的5000次均分值得到拟合值来画图。）  
C.通过用户输入拟合数据的大小以及文件名称，拟合数据将直接保存在目录下。（我一样是通过在区间[-5,5]上的5000次均分值得到拟合值画图，但是得到的数据有误，但我觉得我的a，b矩阵定义没有错，解出的答案也没有错误，错误在B的定义。）  


### -结果
本次作业，我觉得难度较大，其中我在PPForm的函数值返回花了大量时间来一一检查我的错误，但也因为自己时间管理不妥当，导致作业无法如期完成。

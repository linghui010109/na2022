# NonlinearEquationSolving

### -架构

#### 抽象类基类(Abstract base class EquationSolver)
1.void setfun():输入函数的名字，以便调用函数。  
2.void setnum():用户输入以设l,h,D,B1的值，并取得函数A，B，C，E得之。  
3.double function():输入x的指以及对应函数的名字，以查找对应的函数并得到其返回值。  
4.void seta():设置初始的左区间或者初始值。  
5.void setb():设置初始的右区间或者初始值。  

#### 派生类（derived class of EquationSolver)
1.BisectionMethod(二分法）：通过求端点与中点的函数值，判断根的存在区间，使区间长度不断缩减一半，最终收敛到根的求根方法。  

2.Newton Method（牛顿法）：利用迭代点x0处的一阶导数(梯度)和二阶导数(Hessen矩阵)对目标函数进行二次函数近似，然后把二次模型的极小点作为新的迭代点，并不断重复这一过程，直至求得满足精度的近似极小值。  

3.Secant Method(割线法）：通过弦的斜率近似代替目标函数的切线斜率，并用割线与横轴交点的横坐标作为方程式的根的近似的方法。  

### -测试过程

我将所有目标函数都命名为func_i(i=1,2,3...)目标函数的导数命名为derivfunc_i(i-1,2..),顺序是根据课本的题目的顺序。      

对于二分法，在int mian（）中设置一个用来调用BisectionMethod的一个Bisx（x=1，2...),再通过Bisx.seta（）和Bisx.setb（）设置初始区间,以及Bisx.setfun（）设置目标函数的名字。最后通过Bisx.slove()得到函数的根。      

对于牛顿法，在int main（）中设置一个用来调用NewtonMethod的一个Newx（x=1，2...),再通过Bisx.seta（）设置初始值，以及Newx.setfun（）设置目标函数以及目标函数导数的名字。最后通过Bisx.slove()得到函数的根。     
  
对于割线法，在int mian（）中设置一个用来调用SecantMethod的一个Secx（x=1，2...),再通过Secx.seta（）和Secx.setb（）设置初始的两个点,以及Secx.setfun（）设置目标函数的名字。最后通过Secx.slove()得到函数的根。       

# Polynomial Interpolation

### -架构

##### 基类(class Interpolation)
1.void getdata():用户输入节点的个数以及输入各个节点。


##### 派生类(NewtonInterpolation)
1.getNewtonTable():输出差商表。  
2.getInterpolation():输出多项式在某个点的函数值。  
3.getInterpolationPolynomail():输出插值多项式。  
  
 
##### 派生类(HermiteInterpolation)
1.getHermiteTable():输出差商表。  
2.getInterpolation():输出多项式在某个点的函数值。  
3.getInterpolationPolynomail():输出插值多项式。  

### -测试过程

通过用户输入节点的个数以及节点来计算出插值多项式。

对于作业题B，C用户只需输入n即可得到对应的差商表以及多项式。


### -结果

对于本次作业的Hermite插值多项式，我认为我的错误在于只能使用一次微分的节点，对于多次的微分的节点，我的差商表就会出现错误，这是我需要改进的地方。

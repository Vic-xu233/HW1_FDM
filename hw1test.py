import numpy as np

%%一阶微分 一阶精度 向前差分
void 
def u(x):
    u=x+x^2+x^3
    return u

def f_diff(x,dx): %%前向差分，一阶精度
    u_fdiff=(u(x+dx)-u(x))/dx
    return u_diff 

def c_diff(x,dx): %%中间差分，二阶精度
    u_cdiff=(u(x+dx)-u(x-dx))/(2*x)
    return u_cdiff

def tri_c_diff(x,dx):
    xi_1=x+dx
    xi_2=x+2*dx
    xi_n1=x-dx
    u_tcdiff=-u(xi_2)+6*u(xi_1)-3*u(x)-2*u(xi_n1)
    u_tcdiff=u_tcdiff/(6*dx)
    return u_tcdiff


    


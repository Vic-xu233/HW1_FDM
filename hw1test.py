import numpy as np

##一阶微分 一阶精度 向前差分
def u(x):
    u=3*x+4*x**2+5*x**3
    #u=sin(x)
    return u
def u_rdiff(x):
    u=3+8*x+15*x**2
    #u=cos(x)
    return u
def u_2rdiff(x):
    u=8+30*x
    #u=-sin(x)
    return u

def f_diff(x,dx):                   ##前向差分一阶精度
    u_fdiff=(u(x+dx)-u(x))/dx
    return u_diff 

def c_diff(x,dx):                    ##中间差分，二阶精度
    u_cdiff=(u(x+dx)-u(x-dx))/(2*x)
    return u_cdiff

def tri_c_diff(x,dx):         ##四格点  三阶精度
    xi_1=x+dx
    xi_2=x+2*dx
    xi_n1=x-dx
    u_tcdiff=-u(xi_2)+6*u(xi_1)-3*u(x)-2*u(xi_n1)
    u_tcdiff=u_tcdiff/(6*dx)
    return u_tcdiff

def double_diff(x,dx):         ##三格点 二阶精度
    u_double_diff=(u(x+dx)-2*u(x)+u(x-dx))/(dx*dx)
    return u_double_diff

x=1.0
deltax=np.linspace(1e-9,0.1,1000)
errors_diff=[]
errors_2diff=[]

for dx in deltax:
    errors_diff.append(c_diff(x,dx)-u_rdiff(x))
    errors_2diff.append(double_diff(x,dx)-u_2rdiff(x))
    



    


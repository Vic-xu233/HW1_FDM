import numpy as np
import matplotlib.pyplot as plt
flag=1  ##选择对哪个函数进行差分
k=2

def u(x):##函数
    match flag:
        case 1:
            u=np.sin(k*x)
        case 2:
            u=x**2
        case 3:
            u=x**3
        case _:
            u=x
    return u
def u_rdiff(x):##函数解析一阶导数
    match flag:
        case 1:
            u=k*np.cos(k*x)
        case 2:
            u=2*x
        case 3:
            u=3*x**2
        case _:
            u=1
    return u
def u_2rdiff(x):##函数解析二阶导数
    match flag:
        case 1:
            u=-k*k*np.sin(k*x)
        case 2:
            u=2
        case 3:
            u=6*x
        case _:
            u=0
    return u

def f_diff(x,dx):                   ##前向差分一阶精度
    u_fdiff=(u(x+dx)-u(x))/dx
    return u_fdiff 

def c_diff(x,dx):                    ##中间差分，二阶精度
    u_cdiff=(u(x+dx)-u(x-dx))/(2*dx)
    return u_cdiff

def double_diff(x,dx):         ##三格点 二阶精度
    u_double_diff=(u(x+dx)-2*u(x)+u(x-dx))/(dx*dx)
    return u_double_diff

def dou_tri_diff(x,dx):         ##四格点  1阶精度
    xi_1=x+dx
    xi_2=x+2*dx
    xi_n1=x-dx
    u_dtcdiff=-u(xi_2)+4*u(xi_1)-5*u(x)+2*u(xi_n1)
    u_dtcdiff=u_dtcdiff/(dx*dx)
    return u_dtcdiff

def cal(x, deltax):
    errors_cdiff = []
    errors_fdiff = []
    errors_2cdiff = []
    errors_2fdiff = []

    for dx in deltax:
        u1_real = u_rdiff(x)
        u2_real = u_2rdiff(x)
        errors_cdiff.append(c_diff(x, dx) - u1_real)               #2阶
        errors_fdiff.append(f_diff(x, dx) - u1_real)               #1阶
        errors_2cdiff.append(double_diff(x, dx) - u2_real)         #2阶
        errors_2fdiff.append(dou_tri_diff(x, dx) - u2_real)        #1阶

    return errors_cdiff, errors_fdiff, errors_2cdiff, errors_2fdiff



x=1.0
deltax=np.linspace(1e-5,0.1,10000)
errors_cdiff, errors_fdiff, errors_2cdiff, errors_2fdiff = cal(x, deltax)

data = np.column_stack((deltax, errors_cdiff, errors_fdiff, errors_2cdiff, errors_2fdiff))
np.savetxt('errors_output_with_dx.txt', data,
           header='deltax   errors_cdiff   errors_fdiff   errors_2cdiff   errors_2fdiff')


x0, x1 = deltax[4], deltax[-5]
y0, y1 = errors_cdiff[4], errors_cdiff[-5]
y2, y3 = errors_2fdiff[4], errors_2fdiff[-5]

plt.rcParams['font.family'] = 'SimSun' ##中文

#二阶精度误差图
plt.figure()  # 显式创建新图
plt.plot(deltax, errors_cdiff, label='1st C-Derivative Error', color='blue')
plt.plot(deltax, errors_2cdiff, label='2nd C-Derivative Error', color='red')
# 添加图例、标题、标签等
plt.title('Error vs. Δx')
plt.xlabel('Δx (Step size)')
plt.ylabel('Error')
plt.legend()
plt.grid(True)
plt.savefig("二阶精度error_comparison.svg", bbox_inches='tight')
# 一阶精度误差图
plt.figure()  # 显式创建新图

plt.plot(deltax, errors_2fdiff, label='2nd F-Derivative Error(1阶精度的二阶微分误差)', color='orange')
plt.plot(deltax, errors_fdiff, label='1st F-Derivative Error(1阶精度的一阶微分误差)', color='green')
plt.plot(deltax, errors_cdiff, label='1st C-Derivative Error(2阶精度一阶微分)', color='blue')
plt.plot(deltax, errors_2cdiff, label='2nd C-Derivative Error(2阶精度二阶微分)', color='red')
plt.plot([x0, x1], [y2, y3], '--', color='gray', label='辅助直线1')
plt.plot([x0, x1], [errors_fdiff[4], errors_fdiff[-5]], '--', color='gray', label='辅助直线2')
# 添加图例、标题、标签等
plt.title('Error2 vs. Δx')
plt.xlabel('Δx (Step size)')
plt.ylabel('Error')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("error_comparison.svg", bbox_inches='tight')
plt.show()
    


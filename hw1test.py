import numpy as np
import matplotlib.pyplot as plt
flag=4  ##选择对哪个函数进行差分
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
        errors_cdiff.append(c_diff(x, dx) - u1_real)
        errors_fdiff.append(f_diff(x, dx) - u1_real)
        errors_2cdiff.append(double_diff(x, dx) - u2_real)
        errors_2fdiff.append(dou_tri_diff(x, dx) - u2_real)

    return errors_cdiff, errors_fdiff, errors_2cdiff, errors_2fdiff



x=1.0
deltax=np.linspace(1e-7,0.1,100)
errors_cdiff, errors_fdiff, errors_2cdiff, errors_2fdiff = cal(x, deltax)

data = np.column_stack((deltax, errors_cdiff, errors_fdiff, errors_2cdiff, errors_2fdiff))
np.savetxt('errors_output_with_dx.txt', data,
           header='deltax   errors_cdiff   errors_fdiff   errors_2cdiff   errors_2fdiff')


x=np.float32(1.0)
deltax=np.linspace(1e-7,0.1,100,dtype=np.float32)
errors_cdiff, errors_fdiff, errors_2cdiff, errors_2fdiff = cal(x, deltax)

data = np.column_stack((deltax, errors_cdiff, errors_fdiff, errors_2cdiff, errors_2fdiff))
np.savetxt('errors_output_with_dx_float32.txt', data,
           header='deltax   errors_cdiff   errors_fdiff   errors_2cdiff   errors_2fdiff')
print("精度=",deltax.dtype,errors_cdiff[3].dtype,errors_2fdiff[2].dtype)


# 转换为 NumPy 数组（如果你用的是列表）
errors_diff = np.array(errors_cdiff)
errors_2diff = np.array(errors_2cdiff)



# 第一阶导数误差图
'''
plt.plot(deltax, errors_diff, label='1st Derivative Error', color='blue')

# 第二阶导数误差图
#plt.plot(deltax, errors_2diff, label='2nd Derivative Error', color='red')



# 添加图例、标题、标签等
plt.title('Error vs. Δx')
plt.xlabel('Δx (Step size)')
plt.ylabel('Error')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
'''

    


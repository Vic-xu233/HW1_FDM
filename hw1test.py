import numpy as np
import matplotlib.pyplot as plt
flag=1
##一阶微分 一阶精度 向前差分
def u(x):
    match x:
        case 1:
            
        case 2:
            print("这是2")
        case _:
            print("默认情况")
    u=x**3
    #u=sin(x)
    return u
def u_rdiff(x):
    u=3*x**2
    #u=cos(x)
    return u
def u_2rdiff(x):
    u=6*x
    #u=-sin(x)
    return u

def f_diff(x,dx):                   ##前向差分一阶精度
    u_fdiff=(u(x+dx)-u(x))/dx
    return u_diff 

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
    u_dtcdiff=-u(xi_2)+4*u(xi_1)-5*u(x)-2*u(xi_n1)
    u_dtcdiff=u_dtcdiff/(dx*dx)
    return u_dtcdiff

x=1.0
deltax=np.linspace(1e-7,0.1,10000)
errors_diff=[]
errors_2diff=[]

for dx in deltax:
    print(c_diff(x,dx))
    errors_diff.append(c_diff(x,dx)-u_rdiff(x))
    errors_2diff.append(double_diff(x,dx)-u_2rdiff(x))
    
#print(errors_diff)


# 转换为 NumPy 数组（如果你用的是列表）
errors_diff = np.array(errors_diff)
errors_2diff = np.array(errors_2diff)



# 第一阶导数误差图
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

    


'''
所有变量均为32位精度
'''


import numpy as np
import matplotlib.pyplot as plt
flag=1  ##选择对哪个函数进行差分
k=np.float32(2.0)

def u(x):##函数
    match flag:
        case 1:
            u=np.sin(k*x)       #均有截断误差
        case 2:
            u=x**np.float32(2.0)              #前向一次微分有截断误差，其余没有
        case 3:
            u=x**np.float32(3.0)              #中心二次微分无截断误差，其余有
        case _:
            u=x                  #均无截断误差
    return u
def u_rdiff(x):##函数解析一阶导数
    match flag:
        case 1:
            u=k*np.cos(k*x)
        case 2:
            u=np.float32(2.0)*x
        case 3:
            u=np.float32(3.0)*x**np.float32(2.0)
        case _:
            u=np.float32(1.0)
    return u
def u_2rdiff(x):##函数解析二阶导数
    match flag:
        case 1:
            u=-k*k*np.sin(k*x)
        case 2:
            u=np.float32(2.0)
        case 3:
            u=np.float32(6.0)*x
        case _:
            u=np.float32(0.0)
    return u

def f_diff(x,dx):                   ##前向差分一阶精度
    u_fdiff=(u(x+dx)-u(x))/dx
    return u_fdiff 

def c_diff(x,dx):                    ##中间差分，二阶精度
    u_cdiff=(u(x+dx)-u(x-dx))/(np.float32(2.0)*dx)
    return u_cdiff

def double_diff(x,dx):         ##三格点 二阶精度
    u_double_diff=(u(x+dx)-np.float32(2.0)*u(x)+u(x-dx))/(dx*dx)
    return u_double_diff

def dou_tri_diff(x,dx):         ##四格点  1阶精度
    xi_1=x+dx
    xi_2=x+np.float32(2.0)*dx
    xi_n1=x-dx
    u_dtcdiff=-u(xi_2)+np.float32(4.0)*u(xi_1)-np.float32(5.0)*u(x)+np.float32(2.0)*u(xi_n1)
    u_dtcdiff=u_dtcdiff/(dx*dx)
    return u_dtcdiff

def cal(x, deltax):          #计算误差并存储在列表中
    errors_cdiff = []        #中间差分误差
    errors_fdiff = []       #前向差分误差
    errors_2cdiff = []        #中心二阶差分误差
    errors_2fdiff = []        #一阶精度二阶差分误差

    for dx in deltax:
        u1_real = u_rdiff(x)
        u2_real = u_2rdiff(x)
        x1=c_diff(x, dx) - u1_real
        errors_cdiff.append(x1)                                    #2阶
        errors_fdiff.append(f_diff(x, dx) - u1_real)               #1阶
        errors_2cdiff.append(double_diff(x, dx) - u2_real)         #2阶 2次导
        errors_2fdiff.append(dou_tri_diff(x, dx) - u2_real)        #1阶 2次导

    errors_cdiff = np.abs(errors_cdiff)     #取绝对值便于画图
    errors_fdiff = np.abs(errors_fdiff)
    errors_2cdiff = np.abs(errors_2cdiff)
    errors_2fdiff = np.abs(errors_2fdiff)
    '''
    print("errors_2fdiff==",errors_2fdiff.dtype)
    print("fdiff=",errors_2fdiff.dtype)
    print("errors_2cdiff==",errors_2fdiff.dtype)
    print("cdiff==",errors_2fdiff.dtype)
    '''
    return errors_cdiff, errors_fdiff, errors_2cdiff, errors_2fdiff

x=np.float32(1.0)
deltax=np.linspace(1e-7,1e-3,10000,dtype=np.float32)
start=deltax[0]
stop=deltax[-1]
dx=deltax[1]-deltax[0]
# 创建一个 2x2 的子图布局
fig, axs = plt.subplots(2, 2, figsize=(14, 10))
plt.rcParams['font.family'] = 'SimSun'  # 中文字体

for idx, f in enumerate(range(1, 5)):
    flag = f
    errors_cdiff, errors_fdiff, errors_2cdiff, errors_2fdiff = cal(x, deltax)
    
    # 保存对应数据
    filename =f"新的errors_output_with_dx{dx}_flag{flag}.txt"
    data = np.column_stack((deltax, errors_cdiff, errors_fdiff, errors_2cdiff, errors_2fdiff))
    np.savetxt(filename, data,
               header='deltax   errors_cdiff   errors_fdiff   errors_2cdiff   errors_2fdiff')

    # 找到对应子图位置
    ax = axs[idx // 2][idx % 2]
    ax.loglog(deltax, errors_2fdiff, label='1阶精度 二阶微分', color='orange')
    ax.loglog(deltax, errors_fdiff, label='1阶精度 一阶微分', color='green')
    ax.loglog(deltax, errors_cdiff, label='2阶精度 一阶微分', color='blue')
    ax.loglog(deltax, errors_2cdiff, label='2阶精度 二阶微分', color='red')
    
    ax.set_title(f'flag = {flag}')
    ax.set_xlabel('Δx')
    ax.set_ylabel('errors')
    ax.legend()
    ax.grid(True)

plt.suptitle(f"单精度差分误差随 Δx({start},{stop}) 变化（4类函数对比)dx ≈ {dx:.2e}", fontsize=16)
plt.tight_layout(rect=[0, 0, 1, 0.96])  # 给标题留空间
plt.savefig(f"n单精度误差变化图,Δx({start},{stop}),dx ≈ {dx:.2e}.png",dpi=600, bbox_inches='tight')

plt.show()


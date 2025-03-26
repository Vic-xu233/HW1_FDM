import numpy as np
import matplotlib.pyplot as plt
flag=1  ##选择对哪个函数进行差分
k=2

def u(x):##函数
    match flag:
        case 1:
            u=np.sin(k*x)       #均有截断误差
        case 2:
            u=x**2              #前向一次微分有截断误差，其余没有
        case 3:
            u=x**3              #中心二次微分无截断误差，其余有
        case _:
            u=x                  #均无截断误差
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

    errors_cdiff = np.abs(errors_cdiff)
    errors_fdiff = np.abs(errors_fdiff)
    errors_2cdiff = np.abs(errors_2cdiff)
    errors_2fdiff = np.abs(errors_2fdiff)

    return errors_cdiff, errors_fdiff, errors_2cdiff, errors_2fdiff



x=1.0
deltax=np.linspace(1e-7,1e-3,int(1e4))
start
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
    '''
    ax.plot(deltax, errors_2fdiff, label='1阶精度 二阶微分', color='orange')
    ax.plot(deltax, errors_fdiff, label='1阶精度 一阶微分', color='green')
    ax.plot(deltax, errors_cdiff, label='2阶精度 一阶微分', color='blue')
    ax.plot(deltax, errors_2cdiff, label='2阶精度 二阶微分', color='red')
    '''
    ax.loglog(deltax, errors_2fdiff, label='1阶精度 二阶微分', color='orange')
    ax.loglog(deltax, errors_fdiff, label='1阶精度 一阶微分', color='green')
    ax.loglog(deltax, errors_cdiff, label='2阶精度 一阶微分', color='blue')
    ax.loglog(deltax, errors_2cdiff, label='2阶精度 二阶微分', color='red')
    ax.set_title(f'flag = {flag}')
    ax.set_xlabel('Δx')
    ax.set_ylabel('errors')
    ax.legend()
    ax.grid(True)

plt.suptitle(f"差分误差随 Δx 变化（4类函数对比)dx ≈ {dx:.2e}", fontsize=16)
plt.tight_layout(rect=[0, 0, 1, 0.96])  # 给标题留空间
plt.savefig(f"误差变化图,dx ≈ {dx:.2e}.svg", bbox_inches='tight')
plt.show()

'''
##32位精度
x=np.float32(1.0)
deltax=np.linspace(1e-7,0.1,100,dtype=np.float32)
errors_cdiff, errors_fdiff, errors_2cdiff, errors_2fdiff = cal(x, deltax)

data = np.column_stack((deltax, errors_cdiff, errors_fdiff, errors_2cdiff, errors_2fdiff))
np.savetxt('errors_output_with_dx_float32.txt', data,
           header='deltax   errors_cdiff   errors_fdiff   errors_2cdiff   errors_2fdiff')
print("精度=",deltax.dtype,errors_cdiff[3].dtype,errors_2fdiff[2].dtype)

'''
# 转换为 NumPy 数组（如果你用的是列表）
#errors_diff = np.array(errors_cdiff)
#errors_2diff = np.array(errors_2cdiff)
'''
x0, x1 = deltax[4], deltax[-5]
y0, y1 = errors_cdiff[4], errors_cdiff[-5]
y2, y3 = errors_2fdiff[4], errors_2fdiff[-5]

plt.rcParams['font.family'] = 'SimSun' ##中文
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
plt.show()
'''
    


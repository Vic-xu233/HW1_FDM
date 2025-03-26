# 计算流体力学HW1 数值微分误差分析

本项目用于分析不同差分格式（前向差分、中心差分、中心二阶差分等）在一阶和二阶导数近似中的误差随步长 Δx 的变化情况，支持多种函数类型和双/单精度计算，可视化展示截断误差与舍入误差的交汇规律。

---

## 📁 文件结构

| 文件名 | 功能说明 |
|--------|----------|
| `hw1test.py` | 对应作业第二题，使用中心/前向/多点差分方法计算不同精度的一阶与二阶导数误差，并绘图保存结果（双精度） |
| `HW1Q3.py` | 对应作业第三题，批量绘制 log-log 图，分析误差阶数随 Δx 变化的趋势（含不同函数类型 flag=1~4） |
| `HW1Q4.py` | 对应第四题，所有变量强制转为 `float32`，分析单精度误差与数值波动差异 |

---

## 📌 支持的函数类型

通过 `flag` 参数控制所用的函数 `u(x)`：

- `flag = 1`：u(x) = sin(kx)
- `flag = 2`：u(x) = x²
- `flag = 3`：u(x) = x³
- `flag = 4`：u(x) = x

---

## 📐 差分方法说明

- `f_diff(x, dx)`：前向差分（1阶精度）
- `c_diff(x, dx)`：中心差分（2阶精度）
- `double_diff(x, dx)`：中心二阶导数（2阶精度）
- `dou_tri_diff(x, dx)`：四点格式二阶导数（1阶精度）

---

## 🛠️ 调试说明（Debug Notes）

每次手动修改以下代码行：

```python
deltax = np.linspace(1e-7, 1e-3, 10000)
用于验证在不同 Δx 范围和不同步长（dx）下，各种差分方法的误差行为。
```

### 🔍 调整目标：

- **控制 Δx 的范围**：例如 `(1e-7, 1e-3)`、`(1e-5, 0.1)` 等
- **控制 Δx 的密度**：通过第三个参数修改点数，例如 `10000`、`5000`
- **验证误差趋势是否符合理论精度阶数**：比如是否呈现一阶或二阶收敛行为



## 📊 输出结果

程序将自动生成以下内容：

- `.txt`：误差数据保存为表格（包含 Δx 和各类差分误差）
- `.svg` / `.png`：误差曲线图（普通坐标与 log-log 坐标），比较不同差分格式的数值精度表现

---

## ✅ 使用方法

1. 安装依赖：
   ```bash
   pip install numpy matplotlib
   ```

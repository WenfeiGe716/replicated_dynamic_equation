import sympy as sp
# 参数	含义
# Ps	供应商不进行绿色生产经营的正常收益
# Pm	制造商不进行绿色生产经营的正常收益
# Cs	供应商投入绿色成本
# Cm 	制造商投入绿色成本
# r	表示消费者绿色度偏好系数
# Q1	供应商实施绿色化生产的产品绿色度
# Q2	制造商实施绿色化生产的产品绿色度
# u	供应商实施绿色化对最终产品绿色度的贡献程度
# 1-u	制造商实施绿色化对最终产品绿色度的贡献程度
# EPs	供应商实施绿色化生产的超额收益
# EPm	制造商实施绿色化生产的超额收益
# Es	供应商因制造商实施绿色化生产导致声誉提升、品牌溢价等收益
# Em	制造商因供应商实施绿色化生产导致声誉提升、品牌溢价等收益
# M	政府对进行绿色生产经营的供应链企业实施补贴

print("双方演化博弈：")
print("定义变量：x,y,Ps,Pm,Cs,Cm,r,Q1,Q2,u,EPs,EPm,Es,Em,M")
# 定义变量
x,y,Ps,Pm,Cs,Cm,r,Q1,Q2,u,EPs,EPm,Es,Em,M,Gx,Gy,x1,y1 = sp.symbols('x,y,Ps,Pm,Cs,Cm,r,Q1,Q2,u,EPs,EPm,Es,Em,M,Gx,Gy,x1,y1')
# 供应商参与的收益
# U1s = y*(Ps + EPs + r*u*Q1*Ps + r*(1-u)*Q2*Ps -Cs + u*Q1*M) + (1-y)*(Ps+EPs + r*u*Q1*Ps -Cs + u*Q1*M)
U1s = y*(Ps + EPs + r*u*Q1*Ps + r*(1-u)*Q2*Ps -Cs + u*Q1*M) + (1-y)*(Ps+EPs + r*u*Q1*Ps -Cs + u*Q1*M)
U1s = sp.simplify(U1s)
# 供应商不参与的收益
U2s = y*(Ps + Es + r*(1-u)*Q2*Ps + u*Q1*M) + (1-y)*Ps
U2s = sp.simplify(U2s)
# 供应商的平均收益
Us = x*U1s + (1-x)*U2s
Us = sp.simplify(Us)
# 供应商的复制动态方程，格式：f(x) = x(U1s - U1)
fx = x*(U1s-Us)
fx = sp.simplify(fx)
# 因式化简
fx = sp.factor(fx)
print("F(x)=", fx)

# 制造商参与的收益
U1m = x*(Pm + EPm + r*u*Q1*Pm + r*(1-u)*Q2*Pm -Cm + (1-u)*Q2*M) + (1-x)*(Pm + EPm + r*(1-u)*Q2*Pm -Cm + (1-u)*Q2*M)
# 制造商不参与的收益
U2m = x*(Pm + Em + r*u*Q1*Pm) + (1-x)*Pm
# 制造商的平均收益
Um = y*U1m + (1-y)*U2m
# 制造商的复制动态方程，格式：f(x) = x(U1m - U1)
fy = y*(U1m-Um)
fy = sp.simplify(fy)
# 因式化简
fy = sp.factor(fy)
print("F(y)=", fy)

print("输出雅可比矩阵：")
jacobian = []
for f in [fx, fy]:
    temp = []
    for i in [x, y]:
        # 求导
        fxx = sp.diff(f, i)
        # 化简
        fxx = sp.simplify(fxx)
        temp.append(fxx)
    jacobian.append(temp)
# 转化为雅可比矩阵矩阵
jacobian_matrix = sp.Matrix(jacobian)
print(jacobian_matrix)
# 求解平衡点
points = [{"x": x, "y": y} for x, y in sp.solve([fx, fy], (x, y))]
print("points:", points)
# 求解行列式与迹
# 为了方便计算，令
# Gy = Cs - EPs + Es*y + M*Q1*u*y - M*Q1*u - Ps*Q1*r*u
# Gx = Cm - EPm + Em*x + M*Q2*u - M*Q2 + Pm*Q2*r*u - Pm*Q2*r
# Gy = 0 时， y = y1
# Gx = 0 时， x = x1
jacobian_matrix_replace = sp.Matrix([[(2*x - 1)*Gy, x*(Es + M*Q1*u)*(x - 1)], [Em*y*(y - 1), (2*y - 1)*Gx]])
points[len(points)-1]["x"] = x1
points[len(points)-1]["y"] = y1
print("points:", points)
# 初始化字典数组
results = []
for point in points:
    result_matrix = jacobian_matrix_replace.subs(point)
    # 计算雅可比矩阵的行列式
    determinant = sp.factor(result_matrix.det())
    # 计算雅可比矩阵的迹
    trace = result_matrix.trace()
    # 将结果添加到字典数组中
    results.append({'Point': point, 'detJ': determinant, 'trJ': trace})
print("results:", results)

if __name__ == '__main__':
    print("")

# See PyCharm help at https://www.jetbrains.com/help/pycharm/

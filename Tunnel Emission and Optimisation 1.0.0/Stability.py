import json
import os
import numpy as np
import sympy as sp
import matplotlib.pyplot as plt

# 目标存储路径
CACHE_DIR = os.path.join(os.path.dirname(__file__), ".cache")
JSON_FILE = os.path.join(CACHE_DIR, "input_data.json")

# 读取 JSON 数据
if os.path.exists(JSON_FILE):
    with open(JSON_FILE, "r") as f:
        input_data = json.load(f)
else:
    raise FileNotFoundError(f"Failed to Read JSON File: {JSON_FILE}")

# 解析 JSON 参数
sigma3 = input_data["sigma3"]
r0 = input_data["tunnelRadius"]
c = input_data["c"]
phi = input_data["phi"]
Em = input_data["deformationModulus"]
v = input_data["poissonRatio"]
Ec = input_data["concreteModulus"]
vc = input_data["concretePoisson"]
tc = input_data["liningThickness"]
sigma_cc = input_data["concreteUCS"]
xi0 = input_data["Xi0"]

# **计算参数**
sigma_cm = 2 * c * np.cos(np.deg2rad(phi)) / (1 - np.sin(np.deg2rad(phi)))
k = (1 + np.sin(np.deg2rad(phi))) / (1 - np.sin(np.deg2rad(phi)))
po = sigma3
Pcr = (2 * po - sigma_cm) / (1 + k)

p_i = sp.symbols('p_i')
uie = r0 * (1 + v) * (po - p_i) / Em
rp = r0 * (2 * (po * (k - 1) + sigma_cm) / ((1 + k) * ((k - 1) * p_i + sigma_cm)))**(1/(k-1))
uip = r0 * (1 + v) * (2 * (1 - v) * (po - Pcr) * (rp/r0)**2 - (1 - 2*v) * (po - p_i)) / Em

rpf = float(rp.subs(p_i, 0))
u_lamda1 = float(uip.subs(p_i, 0))
u_lamda1 = min(u_lamda1, r0)
uf = (1 / 3) * np.exp(-0.15 * rpf / r0) * u_lamda1

# **计算 LDP**
u = sp.symbols('u')
xa = -2 * rpf / 3 * sp.log((1 - u / u_lamda1) / (1 - uf / u_lamda1))
ui0 = float(sp.solve(sp.Eq(xa, xi0), u)[0])

kc = Ec * (r0**2 - (r0 - tc)**2) / ((1 + vc) * ((1 - 2 * vc) * r0**2 + (r0 - tc)**2))
psc_max = sigma_cc / 2 * (1 - (r0 - tc)**2 / r0**2)
strain_max = psc_max / kc
ui_max = ui0 + strain_max * r0

# 计算 u_pcr
u_pcr = float(uie.subs(p_i, Pcr).evalf())

# **计算 pi_ui_max**
if ui_max > u_pcr:
    eqn_pi_ui_max = sp.Eq(uip, ui_max)
else:
    eqn_pi_ui_max = sp.Eq(uie, ui_max)

pi_ui_max = float(sp.nsolve(eqn_pi_ui_max, p_i, Pcr))
print(pi_ui_max)

# **计算 ui_pcr**
ui_pcr = u_pcr - Pcr / kc * r0

# **计算 p_mob**
eqn_p_mob = sp.Eq(uip, p_i * r0 / kc + ui0)
p_mob = float(sp.nsolve(eqn_p_mob, p_i, 0))
u_mob = float(uip.subs(p_i, p_mob).evalf())
print(u_mob)
print(p_mob, Pcr)
FoS = psc_max / p_mob if p_mob > 0 else np.inf
print(FoS)

# **提取 x=10 时 LDP 中的 u 值**
try:
    u_10_solution = sp.solve(sp.Eq(xa, 10), u)
    u_at_x10 = float(u_10_solution[0])
except Exception:
    u_at_x10 = u_lamda1

# **GRC 画图**
p_vals = np.linspace(Pcr, po, 100)
uie_vals = [uie.subs(p_i, p) for p in p_vals]
q_vals = np.linspace(0, Pcr, 100)
uip_vals = [uip.subs(p_i, p) for p in q_vals]

plt.figure(figsize=(6, 5))
plt.plot(uie_vals, p_vals, linewidth=4, label="Ground Reaction Curve")
plt.plot(uip_vals, q_vals, linewidth=4, label="Plastic Zone Response")

SCC_ui = sp.Symbol('SCC_ui')
SCC_p_i = (SCC_ui - ui0) / r0 * kc
SCC_maxed = psc_max + (SCC_ui * 0)
u_vals_1 = np.linspace(ui0, ui_max, 100)
SCC_p_vals = [float(SCC_p_i.subs(SCC_ui, val).evalf()) for val in u_vals_1]
u_vals_2 = np.linspace(ui_max, u_lamda1, 100)
SCC_max_vals = [float(SCC_maxed.subs(SCC_ui, val).evalf()) for val in u_vals_2]

plt.plot(u_vals_1, SCC_p_vals, linewidth=4, label="Support Characteristic Curve_p_i")
plt.plot(u_vals_2, SCC_max_vals, linewidth=4, label="Support Characteristic Curve_maxed", linestyle="dashed")
plt.plot(u_mob, p_mob, "*", markersize=10, label="Intersection Point", color="purple")
plt.text(u_mob, p_mob, f"({u_mob:.3f}, {p_mob:.3f})", fontsize=10, verticalalignment='bottom', horizontalalignment='right')
plt.text(0.1, 0.9, f"FoS: {FoS:.4f}", transform=plt.gca().transAxes, fontsize=12, bbox=dict(facecolor='white', alpha=0.6))
plt.xlabel("Tunnel Deformation, m",fontsize=14)
plt.ylabel("Ground reaction pressure, MPa",fontsize=14)
plt.title("Ground Reaction Curve",fontsize=14)
plt.xlim(0, u_at_x10)  # 修改这里：设置横坐标最大值为x=10对应的u
plt.legend()
plt.grid()
plt.savefig(os.path.join(CACHE_DIR, "stability_plot1.png"))
print("stability_plot1.png Saved Successfully")
plt.close()

# **LDP 画图**
plt.figure(figsize=(6, 5))
eps = 1e-6

# 限定 xa 范围在 -10 到 10
v_vals = np.linspace(uf, u_lamda1 - eps, 300)
xa_vals_filtered = []
v_vals_filtered = []
for val in v_vals:
    try:
        x_val = float(xa.subs(u, val))
        if -10 <= x_val <= 10:
            xa_vals_filtered.append(x_val)
            v_vals_filtered.append(val)
    except Exception:
        continue

# 限定 xb 范围在 -10 到 10
xb = sp.log(u / uf) * r0
u_vals = np.linspace(uf * 1e-3, uf, 300)
xb_vals_filtered = []
u_vals_filtered = []
for val in u_vals:
    try:
        x_val = float(xb.subs(u, val))
        if -10 <= x_val <= 10:
            xb_vals_filtered.append(x_val)
            u_vals_filtered.append(val)
    except Exception:
        continue

plt.plot(xb_vals_filtered, u_vals_filtered, linewidth=4, label="Before tunnel")
plt.plot(xa_vals_filtered, v_vals_filtered, linewidth=4, label="After tunnel")

plt.xlabel("Location from tunnel face, m",fontsize=14)
plt.ylabel("Normalised Deformation, m",fontsize=14)
plt.title("Longitudinal Deformation Profile",fontsize=14)
plt.grid()
plt.legend()
plt.xlim(-12, 12)  # 设置显示范围为 -12 到 12
plt.savefig(os.path.join(CACHE_DIR, "stability_plot2.png"))
print("stability_plot2.png Saved Successfully")
plt.close()

import numpy as np
import sympy as sp
import matplotlib.pyplot as plt


# 解析 JSON 参数
sigma3 = 5.4
r0 = 3.64
c = 0.01
phi = 30.0
Em = 500
v = 0.25
Ec = 30000
vc = 0.2
tc = 0.252
sigma_cc = 20
xi0 = 0

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
print(p_mob,Pcr)
FoS = psc_max / p_mob if p_mob > 0 else np.inf
print(FoS)

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

plt.plot(u_vals_1, SCC_p_vals, linewidth=4, label="SCC_p_i")
plt.plot(u_vals_2, SCC_max_vals, linewidth=4, label="SCC_maxed", linestyle="dashed")
plt.plot(ui_max, p_mob, "*", markersize=10, label="Intersection Point", color="purple")
plt.text(ui_max, p_mob, f"({ui_max:.3f}, {p_mob:.3f})", fontsize=10, verticalalignment='bottom', horizontalalignment='right')
plt.text(0.1, 0.9, f"FoS: {FoS:.4f}", transform=plt.gca().transAxes, fontsize=12, bbox=dict(facecolor='white', alpha=0.6))
plt.xlabel("Tunnel convergence (radial displacement), m")
plt.ylabel("Ground reaction pressure, MPa")
plt.title("Ground Reaction Curve")
plt.legend()
plt.grid()

# **LDP 画图**
plt.figure(figsize=(6, 5))
eps = 1e-6
v_vals = np.linspace(uf, u_lamda1 - eps, 100)
xa_vals = [float(xa.subs(u, val)) for val in v_vals]

xb = sp.log(u / uf) * r0
u_vals = np.linspace(1e-6, uf, 100)
xb_vals = [float(xb.subs(u, val)) for val in u_vals]

plt.plot(xb_vals, u_vals, linewidth=4, label="Before tunnel")
plt.plot(xa_vals, v_vals, linewidth=4, label="After tunnel")

plt.xlabel("Location from tunnel face, m")
plt.ylabel("Normalized convergence (radial displacement), m")
plt.title("Lining Design Parameters")
plt.show()
plt.grid()



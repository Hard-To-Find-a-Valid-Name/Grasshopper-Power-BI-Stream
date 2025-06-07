import math
import pyomo.environ as pyo
import numpy as np
import sympy as sp
from pymoo.core.problem import Problem
from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.optimize import minimize
from pymoo.operators.sampling.rnd import FloatRandomSampling
from pymoo.operators.crossover.sbx import SBX
from pymoo.operators.mutation.pm import PM
from pymoo.termination import get_termination
import matplotlib.pyplot as plt
import json
import os

# 目标存储路径
print("HEY, THIS IS REALLY RUNNING Opt2D!")
CACHE_DIR = os.path.join(os.path.dirname(__file__), ".cache")
JSON_FILE = os.path.join(CACHE_DIR, "input_data.json")
OUTPUT_JSON_FILE = os.path.join(CACHE_DIR, "output_data.json")

try:
    with open(JSON_FILE, "r", encoding="utf-8") as f:
        input_data = json.load(f)
    print("Parameter Read Successfully:", input_data)
except Exception as e:
    print(f"Parameter Read Failed: {e}")
    exit()

# 读取 JSON 数据
if os.path.exists(JSON_FILE):
    with open(JSON_FILE, "r") as f:
        input_data = json.load(f)
    print("Parameter Read Successfully:", input_data)
else:
    raise FileNotFoundError(f"JSON File {JSON_FILE}")



# 解析参数
r0 = float(input_data["tunnelRadius"])
sigma3 = float(input_data["sigma3"])
c = float(input_data["c"])
phi = float(input_data["phi"])
Em = float(input_data["deformationModulus"])
v = float(input_data["poissonRatio"])
Ec = float(input_data["concreteModulus"])
vc = float(input_data["concretePoisson"])
sigma_cc = float(input_data["concreteUCS"])
sigma_cc_2 = float(input_data["concreteUCS_2"])
tc = float(input_data["liningThickness"])
tc_2 = float(input_data["liningThickness_2"])
xi0 = float(input_data["Xi0"])
xi0_2 = float(input_data["Xi0_2"])

for key, value in input_data.items():
    print(f"{key}: {value}, Type: {type(value)}")


model = pyo.ConcreteModel()
model.tc = pyo.Var(bounds=(tc, tc_2))
model.sigma_cc = pyo.Var(bounds=(sigma_cc / 10, sigma_cc_2 / 10))

print("CCM_CO2_fun() :")
print(f", xi0={xi0}, r0={r0}, sigma3={sigma3}")
print(f"c={c}, phi={phi}, Em={Em}, v={v}, Ec={Ec}, vc={vc}")

def CCM_CO2_fun(tc, xi0, sigma_cc, r0, sigma3, c, phi, Em, v, Ec, vc):
    sigma_cc = sigma_cc * 10
    q = np.pi * (r0 ** 2 - (r0 - tc) ** 2)
    concrete_A13 = 5.5474 * sigma_cc + 123.96
    steelreinforceratio = 0.02
    s_A13 = 15622
    EC_A13 = q * (concrete_A13 + steelreinforceratio * s_A13)
    sigma_cm = (2 * c * np.cos(np.deg2rad(phi))) / (1 - np.sin(np.deg2rad(phi)))
    k = (1 + np.sin(np.deg2rad(phi))) / (1 - np.sin(np.deg2rad(phi)))
    po = sigma3
    Pcr = (2 * po - sigma_cm) / (1 + k)

    # GRC
    p_i = sp.Symbol('p_i', real=True, positive=True)
    uie = r0 * (1 + v) * (po - p_i) / Em
    rp = r0 * ((2 * (po * (k - 1) + sigma_cm)) / ((1 + k) * ((k - 1) * p_i + sigma_cm))) ** (1 / (k - 1))
    uip = r0 * (1 + v) * (2 * (1 - v) * (po - Pcr) * (rp / r0) ** 2 - (1 - 2 * v) * (po - p_i)) / Em

    # LDP
    rpf = float(rp.subs(p_i, 0))
    u_lamda1 = float(uip.subs(p_i, 0))
    if u_lamda1 > r0:
        u_lamda1 = r0
    uf = (1 / 3) * np.exp(-0.15 * rpf / r0) * u_lamda1

    u = sp.Symbol('u', real=True)
    xa = -2 * rpf / 3 * sp.log((1 - u / u_lamda1) / (1 - uf / u_lamda1))

    # SCC
    kc = Ec * (r0 ** 2 - (r0 - tc) ** 2) / ((1 + vc) * ((1 - 2 * vc) * r0 ** 2 + (r0 - tc) ** 2))
    psc_max = sigma_cc / 2 * (1 - (r0 - tc) ** 2 / r0 ** 2)  # 最大支护压力

    # Solve u xi0
    eqn_u = sp.Eq(xa, xi0)
    ui0_sol = sp.solve(eqn_u, u)
    ui0 = float(ui0_sol[0]) if ui0_sol else 0  # 提取解

    # Solve p_mob
    eqn_pmob = sp.Eq(uip, p_i * r0 / kc + ui0)

    print("Starting Solving p_mob...")
    print("Solve eqn_pmob:", eqn_pmob)

    try:
        p_mob_sol = sp.nsolve(eqn_pmob, p_i, 0)  # 0 是初始猜测值
        print(f"  Solving Successfully: p_mob_sol = {p_mob_sol}")
        p_mob = float(p_mob_sol.as_real_imag()[0]) if p_mob_sol.as_real_imag()[0] > 0 else 0
    except Exception as e:
        print(f"  Solving Failed: {repr(e)}")

    # u_mob
    u_mob = float(uip.subs(p_i, p_mob))

    # FoS
    FoS = psc_max / p_mob if p_mob > 0 else np.inf

    return FoS, u_mob, EC_A13
class TunnelOptimization(Problem):
    print(" _evaluate()")
    def __init__(self):
        super().__init__(n_var=2,             # 变量数量 (tc, sigma_cc)
                         n_obj=2,             # 目标数量 (FoS, EC_A13)
                         n_constr=0,          # 约束数量（这里没有约束）
                         xl=np.array([tc, sigma_cc / 10]),   # 变量最小值 (tc, sigma_cc)
                         xu=np.array([tc_2, sigma_cc_2 / 10]))    # 变量最大值 (tc, sigma_cc)

    def _evaluate(self, X, out, *args, **kwargs):
        tc_vals = X[:, 0]
        sigma_cc_vals = X[:, 1]

        FoS_list = []
        EC_A13_list = []

        for tc, sigma_cc in zip(tc_vals, sigma_cc_vals):
            FoS, u_mob, EC_A13 = CCM_CO2_fun(tc, xi0, sigma_cc, r0, sigma3, c, phi, Em, v, Ec, vc)
            FoS_list.append(-FoS)  # 负值以最小化
            EC_A13_list.append(EC_A13)

        out["F"] = np.column_stack([FoS_list, EC_A13_list])

plt.ioff()


class MyDisplay:
    def __init__(self):
        # 只在初始化时创建图像
        self.fig, self.ax = plt.subplots(figsize=(8, 6))
        self.save_path = os.path.join(CACHE_DIR, "pareto_front.png")  # 图片存储路径
        self.latest_population = None  # 记录每个优化问题的最新一代种群
        self.latest_pareto = None  # 记录每个优化问题的最新 Pareto 前沿

    def update(self, algorithm, gen=None):
        """ 每个 generation 生成一个新图，并覆盖旧图，只保留最后一代的点 """
        self.ax.clear()
        self.ax.set_xlabel("Factor of Safety (FoS)",fontsize=14)
        self.ax.set_ylabel("Carbon Emissions (A1-A3) (kgCO2e/m)",fontsize=14)

        if gen is not None:
            self.ax.set_title(f"Pareto Front - Generation {gen}",fontsize=14)
        else:
            self.ax.set_title("Pareto Front",fontsize=14)

        # 获取当前种群
        pop_F = algorithm.pop.get("F")
        self.latest_population = (-pop_F[:, 0], pop_F[:, 1])

        # 获取 Pareto Front
        pareto_F = algorithm.opt.get("F") if algorithm.opt is not None else None
        if pareto_F is not None and pareto_F.size > 0:
            self.latest_pareto = (-pareto_F[:, 0], pareto_F[:, 1])

        # 只绘制最新一代的数据
        if self.latest_population:
            self.ax.scatter(*self.latest_population, color='black', s=5, label="Current Population")

        if self.latest_pareto:
            self.ax.scatter(*self.latest_pareto, color='red', edgecolors='black', s=40, label="Pareto Front")

        self.ax.legend()

        # 保存图像
        print("Saving Pareto Front Figure")
        self.fig.savefig(self.save_path, dpi=300)
        print(f"Figure has been saved: {self.save_path}")


algorithm = NSGA2(
    pop_size=80,
    sampling=FloatRandomSampling(),
    crossover=SBX(prob=0.9, eta=15),
    mutation=PM(eta=20),
    eliminate_duplicates=True
)

display = MyDisplay()
termination = get_termination("n_gen", 20)
problem = TunnelOptimization()
print("Starting Optimisation...")
res = minimize(problem, algorithm, termination, seed=1, verbose=True, callback=display.update)
print("Optimisation Finished!")

# **提取 Pareto 前沿**
FoS_pareto = -res.F[:, 0]  # 负号转换回正值
EC_pareto = res.F[:, 1]
tc_pareto = res.X[:, 0]
sigma_cc_pareto = res.X[:, 1]

# **存储优化结果 JSON**
output_data = {
    "FoS_list": FoS_pareto.tolist(),
    "sigma_cc_vals": sigma_cc_pareto.tolist(),
    "tc_vals": tc_pareto.tolist()
}
with open(OUTPUT_JSON_FILE, "w", encoding="utf-8") as f:
    json.dump(output_data, f, indent=4)


fig, ax1 = plt.subplots(figsize=(8, 6))
ax2 = ax1.twinx()
ax1.scatter(FoS_pareto, sigma_cc_pareto * 10, color='black', marker='s', label="Concrete Strength")
ax2.scatter(FoS_pareto, tc_pareto * 1000, color='red', marker='s', label="Lining Thickness")
ax1.set_xlabel("FoS",fontsize=14)
ax1.set_ylabel("Concrete Strength (MPa)", color='black',fontsize=14)
ax2.set_ylabel("Lining Thickness (mm)", color='red',fontsize=14)
ax1.tick_params(axis='y', colors='black')
ax2.tick_params(axis='y', colors='red')
ax1.legend(loc="upper left",fontsize=14)
ax2.legend(loc="upper right",fontsize=14)
plt.title("(c) FoS vs. Concrete Strength & Lining Thickness",fontsize=14)
plt.grid(True)

#  保存图像
pareto_plot_path = os.path.join(CACHE_DIR, "pareto_fos_vs_concrete.png")
plt.savefig(pareto_plot_path, dpi=300, bbox_inches = "tight")
print(f"Saved (c) FoS vs. Concrete Strength & Lining Thickness at: {pareto_plot_path}")

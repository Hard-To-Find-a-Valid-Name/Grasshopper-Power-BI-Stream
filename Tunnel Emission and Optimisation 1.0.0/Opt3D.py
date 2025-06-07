import os
import json
import pyomo.environ as pyo
import numpy as np
import sympy as sp
import matplotlib.pyplot as plt

# pymoo related imports
from pymoo.core.problem import Problem
from pymoo.algorithms.moo.nsga3 import NSGA3
from pymoo.optimize import minimize
from pymoo.operators.sampling.rnd import FloatRandomSampling
from pymoo.operators.crossover.sbx import SBX
from pymoo.operators.mutation.pm import PM
from pymoo.termination import get_termination
import pandas as pd
from pandas.plotting import parallel_coordinates
from sklearn.preprocessing import MinMaxScaler, StandardScaler
from sklearn.decomposition import PCA
from pymoo.util.ref_dirs import get_reference_directions
from scipy.spatial import ConvexHull
from matplotlib.patches import Polygon
import alphashape
from descartes import PolygonPatch
# ================ 0. Read JSON file ==================
CACHE_DIR = os.path.join(os.path.dirname(__file__), ".cache")
JSON_FILE = os.path.join(CACHE_DIR, "input_data.json")
OUTPUT_JSON_FILE = os.path.join(CACHE_DIR, "output_data.json")

if not os.path.exists(JSON_FILE):
    raise FileNotFoundError(f"Can't find JSON File: {JSON_FILE}")

try:
    with open(JSON_FILE, "r", encoding="utf-8") as f:
        input_data = json.load(f)
    print("Successfully Read JSON Parameters:", input_data)
except Exception as e:
    print(f"Failed to Read JSON Parameters: {e}")
    exit()

# Extract fixed parameters from JSON (convert all values to float)
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

print("Parameters loaded:")
for k, v_ in input_data.items():
    print(f"  - {k}: {v_} (type={type(v_)})")

model = pyo.ConcreteModel()
model.tc = pyo.Var(bounds=(tc, tc_2))
model.sigma_cc = pyo.Var(bounds=(sigma_cc / 10, sigma_cc_2 / 10))
model.xi0 = pyo.Var(bounds=(xi0, xi0_2))


def filter_middle_percent(x, y, percent=90):
    low = (100 - percent) / 2
    high = 100 - low
    y_min, y_max = np.percentile(y, [low, high])
    mask = (y >= y_min) & (y <= y_max)
    return x[mask], y[mask]


def add_convex_hull(ax, x, y, color, alpha=0.2, zorder=0, filter_percent=90):
    if len(x) < 4:
        return

    # 筛选中间百分比的点，去掉极端值
    x, y = np.array(x), np.array(y)
    x_filtered, y_filtered = filter_middle_percent(x, y, percent=filter_percent)

    if len(x_filtered) < 3:
        print(f"[SKIP] Not enough points after filtering ({len(x_filtered)} remaining)")
        return

    points = np.column_stack((x_filtered, y_filtered))
    hull = ConvexHull(points)
    polygon = Polygon(points[hull.vertices], closed=True,
                      facecolor=color, alpha=alpha, edgecolor=None,
                      zorder=zorder, transform=ax.transData)
    ax.add_patch(polygon)
# ================ 1. Define CCM_CO2_fun function (core calculation) ==================
def CCM_CO2_fun(tc, xi0, sigma_cc, r0, sigma3, c, phi, Em, v, Ec, vc):
    sigma_cc = sigma_cc * 10
    Em = 9760.9 * (sigma_cc) ** 0.319
    print(f"Em = {Em}")
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
    print(f"xa = {xa}")
    # SCC
    kc = Ec * (r0 ** 2 - (r0 - tc) ** 2) / ((1 + vc) * ((1 - 2 * vc) * r0 ** 2 + (r0 - tc) ** 2))
    psc_max = sigma_cc / 2 * (1 - (r0 - tc) ** 2 / r0 ** 2)  # 最大支护压力

    # Solve u xi0
    eqn_u = sp.Eq(xa, xi0)
    ui0_sol = sp.solve(eqn_u, u)
    print(f"ui0_sol = {ui0_sol}")
    ui0 = float(ui0_sol[0]) if ui0_sol else 0  # 提取解

    # Solve p_mob
    eqn_pmob = sp.Eq(uip, p_i * r0 / kc + ui0)
    print(f"eqn_pmob = {eqn_pmob}")


    try:
        p_mob_sol = sp.nsolve(eqn_pmob, p_i, 0.1)
        print(f"p_mob_sol = {p_mob_sol}")
        p_mob = float(p_mob_sol.as_real_imag()[0]) if p_mob_sol.as_real_imag()[0] > 0 else 0
    except Exception as e:
        print(f"  Solving Failed: {repr(e)}")
    # u_mob
    u_mob = float(uip.subs(p_i, p_mob))
    print(f"u_mob = {u_mob}")
    # FoS
    FoS = psc_max / p_mob if p_mob > 0 else np.inf
    print(f"FoS = {FoS}")
    return FoS, u_mob, EC_A13



class TunnelOptimization(Problem):
    """
    Problem class for tunnel optimization with 3 objectives and 3 variables.

    Objectives:
    - Maximize Factor of Safety (FoS)
    - Minimize tunnel convergence (u_mob)
    - Minimize carbon emissions (EC_A13)

    Variables:
    - tc: lining thickness
    - sigma_cc: concrete strength
    - xi0: installation location
    """

    def __init__(self):
        super().__init__(
            n_var=3,  # Number of variables (tc, sigma_cc, xi0)
            n_obj=3,  # Number of objectives (FoS, u_mob, EC_A13)
            n_constr=0,  # No constraints
            xl=np.array([tc, sigma_cc / 10, xi0]),  # Lower bounds
            xu=np.array([tc_2, sigma_cc_2 / 10, xi0_2])  # Upper bounds
        )

    def _evaluate(self, X, out, *args, **kwargs):
        FoS_list = []
        u_mob_list = []
        EC_list = []

        for i in range(X.shape[0]):
            tc       = X[i, 0]
            sigma_cc = X[i, 1]
            xi0_val  = X[i, 2]

            try:
                # 调用你的目标函数
                FoS, u_mob, EC_A13 = CCM_CO2_fun(
                    tc, xi0_val, sigma_cc,
                    r0, sigma3, c, phi, Em, v, Ec, vc
                )
            except Exception as e:
                # 如果 CCM_CO2_fun 里抛异常(比如 nsolve失败/log负数等)
                print(f"[ERROR] Exception at i={i}: (tc={tc}, sigma_cc={sigma_cc}, xi0={xi0_val})")
                print("        Exception detail:", repr(e))
                print(">>> Force stopping optimization for debug!")
                raise  # 直接抛出，让程序停止

            # 如果没有抛异常，则检查是否出现NaN/Inf
            if any(np.isnan(x) or np.isinf(x) for x in [FoS, u_mob, EC_A13]):
                print(f"[ERROR] Invalid solution at i={i}: (tc={tc}, sigma_cc={sigma_cc}, xi0={xi0_val})")
                print(f"        FoS={FoS}, u_mob={u_mob}, EC={EC_A13}")
                print(">>> Force stopping optimization for debug!")
                raise RuntimeError("NaN or Inf encountered in objectives.")

            # 如果都正常，就打印下解的情况（可选）
            print(f"[OK] i={i}: (tc={tc:.4f}, sigma_cc={sigma_cc:.4f}, xi0={xi0_val:.4f}) "
                  f"-> FoS={FoS:.4f}, u_mob={u_mob:.6f}, EC={EC_A13:.2f}")

            # pymoo 里你是最小化目标，但你想最大化 FoS，因此加负号

            FoS_list.append(-FoS)
            u_mob_list.append(u_mob)
            EC_list.append(EC_A13)

        out["F"] = np.column_stack([FoS_list, u_mob_list, EC_list])



plt.ioff()  # Turn off interactive plotting


class MyDisplay:
    """
    Display class for visualization of optimization progress and results.
    Creates and updates 3D plots of the Pareto front.
    """

    def __init__(self):
        self.fig = plt.figure(figsize=(8, 6))
        self.save_path = os.path.join(CACHE_DIR, "pareto_front.png")

    def update(self, algorithm, gen=None):
        """Update the visualization with the current generation data"""
        self.fig.clf()  # Clear the figure

        # Create a 3D subplot
        ax = self.fig.add_subplot(111, projection='3d')
        ax.set_xlabel("Factor of Safety",fontsize=16)
        ax.set_ylabel("Tunnel Convergence, m",fontsize=16)
        ax.set_zlabel("Carbon Emissions, kgCO2e",fontsize=16)

        # Set viewing angle
        ax.view_init(elev=30, azim=45)

        if gen is not None:
            ax.set_title(f"Pareto Front - Generation {gen}",fontsize=25)
        else:
            ax.set_title("Pareto Front",fontsize=25)

        # Get current population data
        pop_F = algorithm.pop.get("F")

        # Plot current population (-F[:,0] to convert back to positive FoS)
        x, y, z = -pop_F[:, 0], pop_F[:, 1], pop_F[:, 2]
        ax.scatter(x, y, z, color='black', s=5, alpha=0.5, label="Population")

        # Plot Pareto front if available
        if algorithm.opt is not None:
            pareto_F = algorithm.opt.get("F")
            if pareto_F is not None and pareto_F.size > 0:
                pareto_x = -pareto_F[:, 0]  # Convert back to positive FoS
                pareto_y = pareto_F[:, 1]
                pareto_z = pareto_F[:, 2]

                # Plot Pareto front points
                ax.scatter(pareto_x, pareto_y, pareto_z,
                           color='red', edgecolors='black', s=40, label="Pareto Front")

                # Create 2D projections (optional)
                # xy projection (FoS vs u_mob)
                ax.scatter(pareto_x, pareto_y, min(pareto_z),
                           color='blue', alpha=0.3, s=20)

                # xz projection (FoS vs EC_A13)
                ax.scatter(pareto_x, min(pareto_y), pareto_z,
                           color='green', alpha=0.3, s=20)

                # yz projection (u_mob vs EC_A13)
                ax.scatter(min(pareto_x), pareto_y, pareto_z,
                           color='purple', alpha=0.3, s=20)

        ax.legend(fontsize=15)
        print("Saving 3D Pareto Front Figure...")
        self.fig.tight_layout()
        self.fig.savefig(self.save_path, dpi=600,)
        print(f"Figure saved at: {self.save_path}")
n_obj = 3
ref_dirs = get_reference_directions("das-dennis", n_obj, n_partitions=12)
# Initialize the NSGA2 algorithm
algorithm = NSGA3(
    pop_size=60,
    sampling=FloatRandomSampling(),
    crossover=SBX(prob=0.9, eta=10),  # Use eta=10 for better exploration
    mutation=PM(eta=10),  # Use eta=10 for better exploration
    eliminate_duplicates=True,
    ref_dirs = ref_dirs
)

display = MyDisplay()
termination = get_termination("n_gen", 40)
problem = TunnelOptimization()


res = minimize(
    problem,
    algorithm,
    termination,
    seed=1,
    verbose=True,
    callback=display.update
)


# Extract optimization results
FoS_pareto = -res.F[:, 0]  # Convert back to positive values
u_mob_pareto = res.F[:, 1]
EC_pareto = res.F[:, 2]
tc_pareto = res.X[:, 0]
sigma_cc_pareto = res.X[:, 1]
xi0_pareto = res.X[:, 2]

# Sort results by FoS for better visualization
sorted_idx = np.argsort(FoS_pareto)
FoS_pareto = FoS_pareto[sorted_idx]
u_mob_pareto = u_mob_pareto[sorted_idx]
EC_pareto = EC_pareto[sorted_idx]
tc_pareto = tc_pareto[sorted_idx]
sigma_cc_pareto = sigma_cc_pareto[sorted_idx]
xi0_pareto = xi0_pareto[sorted_idx]

# Save results to JSON
output_data = {
    "FoS_list": FoS_pareto.tolist(),
    "u_mob_list": u_mob_pareto.tolist(),
    "EC_A13_list": EC_pareto.tolist(),
    "tc_vals": tc_pareto.tolist(),
    "sigma_cc_vals": sigma_cc_pareto.tolist(),
    "xi0_vals": xi0_pareto.tolist()
}

with open(OUTPUT_JSON_FILE, "w", encoding="utf-8") as f:
    json.dump(output_data, f, indent=4)

print(f"Results saved to: {OUTPUT_JSON_FILE}")

# Create additional visualization plots
# 1. 2D projections of the Pareto front
fig, ax1 = plt.subplots(figsize=(8, 6))
ax2 = ax1.twinx()
ax3 = ax1.twinx()
ax3.spines["right"].set_position(("outward", 60))
ax1.scatter(FoS_pareto, sigma_cc_pareto * 10, color='black', marker='s', label="Concrete Strength")
ax2.scatter(FoS_pareto, tc_pareto * 1000, color='red', marker='s', label="Lining Thickness")
ax3.scatter(FoS_pareto, xi0_pareto, color='blue', marker='s', label="Installation Location")
add_convex_hull(ax1, FoS_pareto, sigma_cc_pareto * 10, color='black', alpha=0.15)
add_convex_hull(ax2, FoS_pareto, tc_pareto * 1000, color='red', alpha=0.15)
add_convex_hull(ax3,FoS_pareto, xi0_pareto, color='blue', alpha=0.15)
ax1.set_xlabel("Factor of Safety",fontsize=14)
ax1.set_ylabel("Concrete Strength (MPa)", color='black',fontsize=14)
ax2.set_ylabel("Lining Thickness (mm)", color='red',fontsize=14)
ax3.set_ylabel("Installation Location (m)", color='blue',fontsize=14)
ax1.tick_params(axis='y', colors='black')
ax2.tick_params(axis='y', colors='red')
ax3.tick_params(axis='y', colors='blue')
ax1.legend(loc='upper left', bbox_to_anchor=(0.0, 1.0),fontsize=14)   # 左上角
ax2.legend(loc='upper left', bbox_to_anchor=(0.0, 0.9),fontsize=14)   # 略往下
ax3.legend(loc='upper left', bbox_to_anchor=(0.0, 0.8),fontsize=14)   # 再往下
plt.title("FoS vs. Concrete Strength & Lining Thickness & Installation Location",fontsize=20)
plt.grid(True)
# Save plot
params_plot_path = os.path.join(CACHE_DIR, "pareto_fos_vs_concrete.png")
plt.savefig(params_plot_path, dpi=600, bbox_inches="tight")
print(f"Design variables plot saved at: {params_plot_path}")

# 归一化数据（不改变原始值）
scale_min = 0
scale_max = 25
data = np.column_stack([u_mob_pareto, FoS_pareto, EC_pareto])
normalized_data = (data - data.min(axis=0)) / (data.max(axis=0) - data.min(axis=0)) * (scale_max - scale_min) + scale_min
print("Normalise Finished")

df_focus = pd.DataFrame(normalized_data, columns=["u_mob", "FoS", "EC_A13"])


df_focus = df_focus.copy()
df_focus["Group"] = np.where(
    (df_focus["FoS"] > 17.5) & (df_focus["u_mob"] < 10) & (df_focus["EC_A13"] < 10),
    "Balanced High-FoS",
    "Other"
)

plt.figure(figsize=(8, 6))
parallel_coordinates(df_focus, class_column="Group", colormap=plt.get_cmap("Set1"), linewidth=1.5)
plt.title("High-FoS, Low-Deformation & Low-Carbon Solutions (Normalized to [0–25])",fontsize=14)
plt.grid(True)
plt.savefig(os.path.join(CACHE_DIR, "scaled_parallel_coordinates_0_25.png"), dpi=600, bbox_inches="tight")
save_path = os.path.join(CACHE_DIR, "scaled_parallel_coordinates_0_25.png")
plt.savefig(save_path, dpi=600, bbox_inches="tight")

print(f"Scaled (0-25) parallel coordinates plot saved at: {save_path}")


X_all = np.column_stack([
    FoS_pareto,
    u_mob_pareto,
    EC_pareto,
    tc_pareto * 1000,           # mm
    sigma_cc_pareto * 10,       # MPa
    xi0_pareto                  # m
])

# 标准化
X_scaled = StandardScaler().fit_transform(X_all)

# PCA 投影
pca = PCA(n_components=2)
X_pca = pca.fit_transform(X_scaled)

# 构建 DataFrame（加入 PC1）
df_pareto = pd.DataFrame({
    "u_mob": u_mob_pareto,
    "Lining Thickness (mm)": tc_pareto ,
    "Concrete UCS (MPa)": sigma_cc_pareto ,
    "Support Installation (m)": xi0_pareto,
})


scaler = MinMaxScaler()
df_scaled = pd.DataFrame(scaler.fit_transform(df_pareto.drop(columns=["u_mob"])),
                         columns=df_pareto.columns.drop("u_mob"))
df_scaled["u_mob"] = df_pareto["u_mob"]


df_sorted = df_scaled.sort_values("u_mob").reset_index(drop=True)
df_sorted["Group"] = pd.cut(df_sorted["u_mob"], bins=3, labels=["Low", "Mid", "High"])

custom_colors = ["gold", "teal", "purple"]  # 对应 Low, Mid, High
plt.figure(figsize=(8, 6))
parallel_coordinates(
    df_sorted.drop(columns=["u_mob"]),
    class_column="Group",
    color=custom_colors,
    linewidth=1.5
)
plt.title("Parallel Coordinates Group by Normalised Deformation",fontsize=20)
plt.xticks(fontsize = 14)
plt.grid(True)
plt.legend(fontsize=14)


parcoord_path = os.path.join(CACHE_DIR, "pca_parallel_coordinates_normalized_u_mob.png")
plt.savefig(parcoord_path, dpi=600, bbox_inches="tight")



df_pareto = pd.DataFrame({
    "FoS": FoS_pareto,
    "Lining Thickness (mm)": tc_pareto ,
    "Concrete UCS (MPa)": sigma_cc_pareto ,
    "Support Installation (m)": xi0_pareto,
})


scaler = MinMaxScaler()
df_scaled = pd.DataFrame(scaler.fit_transform(df_pareto.drop(columns=["FoS"])),
                         columns=df_pareto.columns.drop("FoS"))
df_scaled["FoS"] = df_pareto["FoS"]


df_sorted = df_scaled.sort_values("FoS").reset_index(drop=True)
df_sorted["Group"] = pd.cut(df_sorted["FoS"], bins=3, labels=["Low", "Mid", "High"])

custom_colors = [  "purple","teal" ,"gold"]  # 对应 Low, Mid, High
plt.figure(figsize=(8, 6))
parallel_coordinates(
    df_sorted.drop(columns=["FoS"]),
    class_column="Group",
    color=custom_colors,
    linewidth=1.5
)
plt.title("Parallel Coordinates Group by Normalised FoS",fontsize=20)
plt.xticks(fontsize = 10)
plt.grid(True)
plt.legend(fontsize=14)


parcoord_path = os.path.join(CACHE_DIR, "pca_parallel_coordinates_normalized_FoS.png")
plt.savefig(parcoord_path, dpi=600, bbox_inches="tight")

df_pareto = pd.DataFrame({
    "EC": EC_pareto,
    "Lining Thickness (mm)": tc_pareto ,
    "Concrete UCS (MPa)": sigma_cc_pareto ,
    "Support Installation (m)": xi0_pareto,
})


scaler = MinMaxScaler()
df_scaled = pd.DataFrame(scaler.fit_transform(df_pareto.drop(columns=["EC"])),
                         columns=df_pareto.columns.drop("EC"))
df_scaled["EC"] = df_pareto["EC"]


df_sorted = df_scaled.sort_values("EC").reset_index(drop=True)
df_sorted["Group"] = pd.cut(df_sorted["EC"], bins=3, labels=["Low", "Mid", "High"])

custom_colors = ["gold", "teal", "purple"]  # 对应 Low, Mid, High
plt.figure(figsize=(8, 6))
parallel_coordinates(
    df_sorted.drop(columns=["EC"]),
    class_column="Group",
    color=custom_colors,
    linewidth=1.5
)
plt.title("Parallel Coordinates Group by Normalised Emission Carbon",fontsize=20)
plt.xticks(fontsize = 10)
plt.grid(True)
plt.legend(fontsize=14)


parcoord_path = os.path.join(CACHE_DIR, "pca_parallel_coordinates_normalized_EC.png")
plt.savefig(parcoord_path, dpi=600, bbox_inches="tight")

print(f"Normalized PCA-based parallel coordinates plot saved at: {parcoord_path}")
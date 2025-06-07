import sys
import subprocess
import json
import os
import threading
from PyQt6 import QtWidgets
from PyQt6.QtWidgets import QLabel, QVBoxLayout
from PyQt6.QtGui import QPixmap
from PyQt6.QtCore import QTimer, Qt
from ui_edit1 import Ui_MainWindow
from PyQt6.QtWidgets import QMessageBox
import glob
import psutil
import numpy as np



# 目标存储路径
CACHE_DIR = os.path.join(os.path.dirname(__file__), ".cache")
JSON_FILE = os.path.join(CACHE_DIR, "input_data.json")
JSON_FILE1 = os.path.join(CACHE_DIR, "output_data.json")

def clear_cache():
    """清空 .cache 目录下的 JSON 文件和 PNG 图片"""
    if os.path.exists(CACHE_DIR):
        # 删除所有 PNG 图片
        image_files = glob.glob(os.path.join(CACHE_DIR, "*.png"))
        for img_file in image_files:
            try:
                os.remove(img_file)
                print(f"PNG Deleted: {img_file}")
            except Exception as e:
                print(f"PNG Deleted Failed: {img_file}: {e}")

        # 清空 JSON 文件
        if os.path.exists(JSON_FILE):
            try:
                with open(JSON_FILE, "w", encoding="utf-8") as f:
                    json.dump({}, f)  # 写入空 JSON
            except Exception as e:
                print(f"JSON File Deleted Failed: {e}")
        if os.path.exists(JSON_FILE):
            try:
                with open(JSON_FILE1, "w", encoding="utf-8") as f:
                    json.dump({}, f)  # 写入空 JSON
            except Exception as e:
                print(f"JSON File Deleted Failed: {e}")
    else:
        os.makedirs(CACHE_DIR, exist_ok=True)
        print(f"Create .cache Directory: {CACHE_DIR}")


# **在程序启动时清空缓存**
clear_cache()


if hasattr(sys, "_MEIPASS"):  # 检查 `_MEIPASS` 是否存在
    BASE_DIR = sys._MEIPASS
else:
    BASE_DIR = os.path.dirname(os.path.abspath(__file__))



class MainWindow(QtWidgets.QMainWindow,Ui_MainWindow):
    def __init__(self):
        super().__init__()
        self.setupUi(self)
        self.optimized_data = {"FoS_list": []}
        self.input_fields = [...]
        self.input_fields = [self.lineEdit_sigma3, self.lineEdit_r0, self.lineEdit_c,
                             self.lineEdit_phi, self.lineEdit_Em, self.lineEdit_v,
                             self.lineEdit_Ec, self.lineEdit_vc, self.lineEdit_tc,
                             self.lineEdit_sigmacc, self.lineEdit_xi0, self.lineEdit_tc_2,
                             self.lineEdit_sigmacc_2, self.lineEdit_xi0_2]
        self.input_widgets = [self.Optt, self.StabilityPlot1, self.StabilityPlot2, self.Optt1]

        # 绑定 `Optt1` 并添加 QLabel 显示 `(c) FoS vs. Concrete Strength & Lining Thickness`
        self.optt1_widget = self.findChild(QtWidgets.QWidget, "Optt1")
        self.optt1_label = QLabel(self.optt1_widget)
        self.optt1_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.optt1_label.setText("Waiting for optimisation...")
        layout_optt1 = QVBoxLayout(self.optt1_widget)
        layout_optt1.addWidget(self.optt1_label)

        # 定时器刷新 `pareto_fos_vs_concrete.png`
        self.timer_fos_vs_concrete = QTimer(self)
        self.timer_fos_vs_concrete.timeout.connect(self.update_pareto_fos_vs_concrete_image)
        self.timer_fos_vs_concrete.start(5000)  # 每 1 秒刷新一次



         # **绑定输入框控件*
        self.input_sigma3 = self.findChild(QtWidgets.QLineEdit, "lineEdit_sigma3")  # σ3
        self.input_radius = self.findChild(QtWidgets.QLineEdit, "lineEdit_r0")  # Tunnel Radius
        self.input_c = self.findChild(QtWidgets.QLineEdit, "lineEdit_c")  # c
        self.input_phi = self.findChild(QtWidgets.QLineEdit, "lineEdit_phi")  # φ
        self.input_deform_mod = self.findChild(QtWidgets.QLineEdit, "lineEdit_Em")  # Deformation Modulus
        self.input_poisson = self.findChild(QtWidgets.QLineEdit, "lineEdit_v")  # Poisson Ratio
        self.input_concrete_mod = self.findChild(QtWidgets.QLineEdit, "lineEdit_Ec")  # Modulus of Concrete
        self.input_concrete_poisson = self.findChild(QtWidgets.QLineEdit, "lineEdit_vc")  # Concrete Poisson Ratio
        self.input_xi0 = self.findChild(QtWidgets.QLineEdit, "lineEdit_xi0")  # xi0
        self.input_tc = self.findChild(QtWidgets.QLineEdit, "lineEdit_tc")  # Lining Thickness
        self.input_sigmacc = self.findChild(QtWidgets.QLineEdit, "lineEdit_sigmacc")  # Concrete UCS
        self.input_xi0_2 = self.findChild(QtWidgets.QLineEdit, "lineEdit_xi0_2")  # xi0
        self.input_tc_2 = self.findChild(QtWidgets.QLineEdit, "lineEdit_tc_2")  # Lining Thickness
        self.input_sigmacc_2 = self.findChild(QtWidgets.QLineEdit, "lineEdit_sigmacc_2")  # Concrete UCS

        # **绑定优化选择**
        self.radio_2D = self.findChild(QtWidgets.QRadioButton, "Opt2D")
        self.radio_3D = self.findChild(QtWidgets.QRadioButton, "Opt3D")
        self.radio_2D3 = self.findChild(QtWidgets.QRadioButton, "Opt2D3")


        # **绑定按钮**
        self.start_button = self.findChild(QtWidgets.QPushButton, "Start")
        self.start_button.clicked.connect(self.run_optimization)
        self.clean_button = self.findChild(QtWidgets.QPushButton, "Clean")
        self.clean_button.clicked.connect(self.clean_cache_and_inputs)

        # **绑定 `Optt` 并添加 QLabel 显示 Pareto Front**
        self.optt_widget = self.findChild(QtWidgets.QWidget, "Optt")
        self.image_label = QLabel(self.optt_widget)
        self.image_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.image_label.setText("Waiting for the optimisation...")
        layout = QVBoxLayout(self.optt_widget)
        layout.addWidget(self.image_label)

        # **绑定 `StabilityPlot1` 和 `StabilityPlot2`**
        self.stability_plot1 = self.findChild(QtWidgets.QWidget, "StabilityPlot1")
        self.stability_plot2 = self.findChild(QtWidgets.QWidget, "StabilityPlot2")

        self.label_plot1 = QLabel(self.stability_plot1)
        self.label_plot2 = QLabel(self.stability_plot2)

        layout1 = QVBoxLayout(self.stability_plot1)
        layout2 = QVBoxLayout(self.stability_plot2)
        layout1.addWidget(self.label_plot1)
        layout2.addWidget(self.label_plot2)

        # **定时器刷新 `pareto_front.png`**
        self.timer = QTimer(self)
        self.timer.timeout.connect(self.update_pareto_image)
        self.timer.start(1000)
        self.refresh_timer = QTimer(self)
        self.refresh_timer.timeout.connect(self.load_optimized_data)
        self.refresh_timer.start(3000)  # 每 1s 刷新一次
        print("JSON Auto-Refresh Started by refresh_timer")
        self.optimized_data = self.load_optimized_data()
        # **找到 Slicer Bar**
        self.slicer_bar = self.findChild(QtWidgets.QSlider, "slicerBar")
        self.slicer_bar.setMinimum(0)
        self.slicer_bar.setMaximum(len(self.optimized_data["FoS_list"]) - 1)
        # **绑定 Slicer Bar 滑动事件**
        self.plot_button = self.findChild(QtWidgets.QPushButton, "Plot")
        self.plot_button.clicked.connect(self.update_stability_inputs)
        self.plot_button.clicked.connect(self.run_stability)

        #**绑定 展示的Optimisation Path**
        self.FoS_relation = self.findChild(QtWidgets.QRadioButton, "FoS")
        self.umob_relation = self.findChild(QtWidgets.QRadioButton, "umob")
        self.EC_relation = self.findChild(QtWidgets.QRadioButton, "EC")
        self.update_button = self.findChild(QtWidgets.QPushButton, "Update")
        self.update_button.clicked.connect(self.update_pareto_fos_vs_concrete_image)

    def clean_cache_and_inputs(self):
        """清空 .cache 目录（JSON 文件、优化图片），清空所有输入框，并重置所有 Widget"""

        # **清空 .cache 目录的 PNG 图片和 JSON 文件**
        if os.path.exists(CACHE_DIR):
            # 删除所有 PNG 图片
            image_files = glob.glob(os.path.join(CACHE_DIR, "*.png"))
            for img_file in image_files:
                try:
                    os.remove(img_file)
                    print(f"PNG Deleted: {img_file}")
                except Exception as e:
                    print(f"PNG Deleted Failed: {img_file}: {e}")

            # 清空 JSON 文件
            if os.path.exists(JSON_FILE):
                try:
                    with open(JSON_FILE, "w", encoding="utf-8") as f:
                        json.dump({}, f)  # 写入空 JSON
                    print(f"JSON Deleted: {JSON_FILE}")
                except Exception as e:
                    print(f"JSON Deleted Failed: {e}")
        else:
            os.makedirs(CACHE_DIR, exist_ok=True)
            print(f"Create .cache Directory: {CACHE_DIR}")

        # **清空所有 QLineEdit 输入框**
        for field in self.input_fields:
            if field:
                field.clear()
        print("All Inputs Cleaned")

        # **清空 QLabel 图像显示**
        self.image_label.clear()
        self.image_label.setText("Waiting for Optimisation Result")

        self.label_plot1.clear()
        self.label_plot2.clear()
        self.label_plot1.setText("No Data")
        self.label_plot2.setText("No Data")

        print("All QLabel Pic Cleaned")

        # **清空 Matplotlib 画布（如果有）**
        if hasattr(self, 'stability_plot1') and hasattr(self, 'stability_plot2'):
            for plot in [self.stability_plot1, self.stability_plot2]:
                for child in plot.children():
                    if isinstance(child, QLabel):
                        child.clear()
                        child.setText("No Data")
        print("Matplotlib Widget Cleaned")

        # **终止所有 Python 进程（除了 main.py）**
        current_pid = os.getpid()  # 获取 main.py 进程 ID
        current_script_name = os.path.basename(__file__)  # 获取当前脚本文件名 (main.py)

        for proc in psutil.process_iter(['pid', 'name', 'cmdline']):
            try:
                if proc.info['name'] == "pythonw.exe" or proc.info['name'] == "python":
                    cmdline = " ".join(proc.info['cmdline']) if proc.info['cmdline'] else ""

                    # 只终止 `Opt2D.py`、`Opt3D.py`、`Stability.py`
                    if any(script in cmdline for script in ["Opt2D.py", "Opt3D.py", "Stability.py", "Opt2D_1"]) and proc.info[
                        'pid'] != current_pid:
                        print(f"Progress Ended: {cmdline} (PID={proc.info['pid']})")
                        proc.terminate()  # 发送终止信号
                        proc.wait(timeout=5)  # 等待进程结束
            except Exception as e:
                print(f"Progress Ended Failed: {e}")

        print("All Progresses have been ended")
        # **弹出消息通知用户清理完成**
        QtWidgets.QMessageBox.information(self, "Clean Finished", "All Data, Figure and Cache Cleaned!")

    def update_pareto_image(self):
        pareto_path = os.path.join(CACHE_DIR, "pareto_front.png")
        if os.path.exists(pareto_path):
            pixmap = QPixmap(pareto_path)
            if not pixmap.isNull():
                self.image_label.setPixmap(pixmap.scaled(
                    self.image_label.width(),
                    self.image_label.height(),
                    Qt.AspectRatioMode.KeepAspectRatio,
                    Qt.TransformationMode.SmoothTransformation
                ))
                print(f"Have Updated Pareto Front: {pareto_path}")
        else:
            self.image_label.setText("Waiting for the optimisation...")

    def update_pareto_fos_vs_concrete_image(self):
        if self.radio_2D.isChecked() or self.radio_2D3.isChecked():
            pareto_path = os.path.join(CACHE_DIR, "pareto_fos_vs_concrete.png")
            if os.path.exists(pareto_path):
                pixmap = QPixmap(pareto_path)
                if not pixmap.isNull():
                    # 建议直接用SmoothTransformation
                    self.optt1_label.setPixmap(pixmap.scaled(
                        self.optt1_label.width(),
                        self.optt1_label.height(),
                        Qt.AspectRatioMode.KeepAspectRatio,
                        Qt.TransformationMode.SmoothTransformation
                    ))
                    print(f" 已更新 (c) FoS vs. Concrete Strength & Lining Thickness: {pareto_path}")

            else:
                self.optt1_label.setText("Waiting for optimisation...")
        else:
            if self.FoS_relation.isChecked():
                pareto_path = os.path.join(CACHE_DIR, "pca_parallel_coordinates_normalized_FoS.png")
                if os.path.exists(pareto_path):
                    pixmap = QPixmap(pareto_path)
                    if not pixmap.isNull():
                        # 建议直接用SmoothTransformation
                        self.optt1_label.setPixmap(pixmap.scaled(
                            self.optt1_label.width(),
                            self.optt1_label.height(),
                            Qt.AspectRatioMode.KeepAspectRatio,
                            Qt.TransformationMode.SmoothTransformation
                        ))
                        print(f" 已更新 (c) FoS vs. Concrete Strength & Lining Thickness: {pareto_path}")

            if self.umob_relation.isChecked():
                pareto_path = os.path.join(CACHE_DIR, "pca_parallel_coordinates_normalized_u_mob.png")
                if os.path.exists(pareto_path):
                    pixmap = QPixmap(pareto_path)
                    if not pixmap.isNull():
                        # 建议直接用SmoothTransformation
                        self.optt1_label.setPixmap(pixmap.scaled(
                            self.optt1_label.width(),
                            self.optt1_label.height(),
                            Qt.AspectRatioMode.KeepAspectRatio,
                            Qt.TransformationMode.SmoothTransformation
                        ))
                        print(f" 已更新 (c) FoS vs. Concrete Strength & Lining Thickness: {pareto_path}")

            if self.EC_relation.isChecked():
                pareto_path = os.path.join(CACHE_DIR, "pca_parallel_coordinates_normalized_EC.png")
                if os.path.exists(pareto_path):
                    pixmap = QPixmap(pareto_path)
                    if not pixmap.isNull():
                        # 建议直接用SmoothTransformation
                        self.optt1_label.setPixmap(pixmap.scaled(
                            self.optt1_label.width(),
                            self.optt1_label.height(),
                            Qt.AspectRatioMode.KeepAspectRatio,
                            Qt.TransformationMode.SmoothTransformation
                        ))
                        print(f" 已更新 (c) FoS vs. Concrete Strength & Lining Thickness: {pareto_path}")

    def collect_input_data(self):
        """收集所有输入框数据并存为 JSON，确保输入值合法，同时清空 JSON 和优化图片"""

        # **确保缓存目录存在**
        os.makedirs(CACHE_DIR, exist_ok=True)

        # **清空 JSON 文件**
        if os.path.exists(JSON_FILE):
            try:
                with open(JSON_FILE, "w", encoding="utf-8") as f:
                    json.dump({}, f)  # 写入空 JSON
                print(f"JSON Cleaned: {JSON_FILE}")
            except Exception as e:
                QMessageBox.critical(self, "Error", f"JSON File Delete Failed: {e}")
                return  # 终止数据收集

        # **清空 .cache 目录下的优化结果图片**
        image_files = glob.glob(os.path.join(CACHE_DIR, "*.png"))  # 获取所有 PNG 图片
        for img_file in image_files:
            try:
                os.remove(img_file)
                print(f"Figure Deleted: {img_file}")
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Figure Deleted Failed {img_file}: {e}")
                return  # 终止数据收集

        data = {}

        for field in [
            ("sigma3", self.input_sigma3),
            ("tunnelRadius", self.input_radius),
            ("c", self.input_c),
            ("phi", self.input_phi),
            ("deformationModulus", self.input_deform_mod),
            ("poissonRatio", self.input_poisson),
            ("concreteModulus", self.input_concrete_mod),
            ("concretePoisson", self.input_concrete_poisson),
            ("Xi0", self.input_xi0),
            ("liningThickness", self.input_tc),
            ("concreteUCS", self.input_sigmacc),
            ("concreteUCS_2",self.input_sigmacc_2),
            ("Xi0_2", self.input_xi0_2),
            ("liningThickness_2", self.input_tc_2),
            ("Xi0i", self.input_xi0),
            ("liningThicknessi", self.input_tc),
            ("concreteUCSi", self.input_sigmacc),
        ]:
            key, line_edit = field
            text = line_edit.text().strip()

            if not text:
                QMessageBox.warning(self, "Input Error", f"Parameter '{key}' Can't be Blank!")
                return  # 终止数据收集，等待用户输入正确值

            try:
                data[key] = float(text)
            except ValueError:
                QMessageBox.warning(self, "Input Error", f"Parameter '{key}' Need to be a Number！")
                return  # 终止数据收集

        # **确保目录存在**
        os.makedirs(CACHE_DIR, exist_ok=True)

        if self.radioButton.isChecked():
            data["algorithm"] = "NSGA2"
        elif self.radioButton_2.isChecked():
            data["algorithm"] = "SPEA2"
        elif self.radioButton_3.isChecked():
            data["algorithm"] = "MOEAD"


        # **写入 JSON 文件**
        with open(JSON_FILE, "w", encoding="utf-8") as f:
            json.dump(data, f, indent=4)

        print(f"Parameters Saved: {JSON_FILE}")



    def run_process(self, script_name):
        if hasattr(sys, '_MEIPASS'):
            script_path = os.path.join(sys._MEIPASS, script_name)
        else:
            script_path = os.path.join(BASE_DIR, script_name)

        if not os.path.exists(script_path):
            print(f" Error: Can't find {script_path}")
            return

        print(f" Running {script_name}...")
        python_path = sys.executable  #
        subprocess.Popen([python_path, script_path])

        print(f" Running {script_name}...")
        print("DEBUG: script_path =", script_path, "exists =", os.path.exists(script_path))
        python_path = os.path.join(os.path.dirname(sys.executable), "pythonw.exe")
        if not os.path.exists(python_path):
            python_path = sys.executable  # 如果 pythonw.exe 不存在，就使用当前 exe
        try:
            process = subprocess.Popen(
                [python_path, script_path],  # 关键：sys.executable 仍然是运行本 exe
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                encoding="utf-8",
                errors="replace"
            )

            print(f" Launched {script_name} (PID: {process.pid})")

            for line in process.stdout:
                print(f" {script_name} Output: {line.strip()}")

            error_output = process.stderr.read()
            if error_output:
                print(f" {script_name} Error:\n{error_output}")
                with open("error.log", "a", encoding="utf-8") as f:
                    f.write(error_output + "\n")

        except Exception as e:
            print(f" Failed to run {script_name}: {str(e)}")

    def run_optimization(self):
        """ 先运行优化程序（Opt2D/Opt2D_1/Opt3D），再运行 Stability.py """
        self.collect_input_data()  # 先收集并写入 input_data.json

        optimization_script = None
        if self.radio_2D.isChecked():
            optimization_script = "Opt2D.py"
        elif self.radio_2D3.isChecked():
            optimization_script = "Opt2D_1.py"
        elif self.radio_3D.isChecked():
            optimization_script = "Opt3D.py"

        if optimization_script:
            print(f" Running optimization: {optimization_script}")
            optimization_thread = threading.Thread(
                target=self.run_process,
                args=(optimization_script,),
                daemon=True
            )
            optimization_thread.start()
    def load_optimized_data(self):
        output_json_path = os.path.join(CACHE_DIR, "output_data.json")
        if os.path.exists(output_json_path):
            with open(output_json_path, "r", encoding="utf-8") as f:
                try:
                    new_data = json.load(f)
                except json.JSONDecodeError:
                    print(f" Error: Failed to load {output_json_path}, file may be corrupted.")
                    return self.optimized_data

            if "FoS_list" not in new_data:
                new_data["FoS_list"] = []

            # 获取当前数据长度

            print("New FoS_list len:", len(new_data["FoS_list"]))
            print("Old FoS_list len:", len(self.optimized_data.get("FoS_list", [])))

            sorted_indices = np.argsort(new_data["FoS_list"])  # 获取排序后的索引
            for key in new_data.keys():
                new_data[key] = np.array(new_data[key])[sorted_indices].tolist()  # 按照索引排序

            if len(new_data["FoS_list"]) != len(self.optimized_data.get("FoS_list", [])):
                self.slicer_bar.setMaximum(len(new_data["FoS_list"]) - 1)
                print("Slicer Bar maximum updated:", self.slicer_bar.maximum())

            self.optimized_data = new_data
            return new_data

        return self.optimized_data

    def update_stability_inputs(self):
        """ 根据 Slicer 选择的索引，更新 Stability.py 计算的参数 """
        index = self.slicer_bar.value()

        if len(self.optimized_data["FoS_list"]) == 0:
            return  # 没有数据时不做任何操作

        if index < 0:
            index = 0
        elif index >= len(self.optimized_data["FoS_list"]):
            index = len(self.optimized_data["FoS_list"]) - 1

        # 读取 Slicer 选中的参数
        sigma_cc = self.optimized_data["sigma_cc_vals"][index]
        tc = self.optimized_data["tc_vals"][index]

        # **检查是否有 xi0**
        if "xi0_vals" in self.optimized_data and len(self.optimized_data["xi0_vals"]) > index:
            xi0 = self.optimized_data["xi0_vals"][index]
            if xi0 is None or np.isnan(xi0):  # `Opt2D.py` 的 xi0 是 NaN
                xi0 = 0  # 使用默认值
        else:
            xi0 = 0  # `Opt2D.py` 没有 xi0，默认设置为 0

        # 更新 input_data.json 让 Stability.py 读取
        input_data_path = os.path.join(CACHE_DIR, "input_data.json")
        with open(input_data_path, "r", encoding="utf-8") as f:
            input_data = json.load(f)

        # 然后更新为新值（可加 round 保留小数）
        input_data.update({
            "concreteUCS": round(sigma_cc * 10, 2),
            "liningThickness": round(tc, 2),
            "Xi0": round(xi0, 2)
        })

        # 保存到 JSON 文件
        with open(input_data_path, "w", encoding="utf-8") as f:
            json.dump(input_data, f, indent=4)

        print(f"Stability Inputs Updated: sigma_cc={sigma_cc}, tc={tc}, xi0={xi0}")

    def run_stability(self):
        """运行 Stability.py 并在完成后更新界面"""

        def stability_worker():
            self.run_process("Stability.py")

            import time
            time.sleep(2)

            plot1_path = os.path.join(CACHE_DIR, "stability_plot1.png")
            plot2_path = os.path.join(CACHE_DIR, "stability_plot2.png")

            if os.path.exists(plot1_path) and os.path.exists(plot2_path):
                print(f" Stability plots found: {plot1_path}, {plot2_path}")
                self.load_stability_plots()
            else:
                print(" Stability plots NOT found!")

        stability_thread = threading.Thread(target=stability_worker, daemon=True)
        stability_thread.start()

    def load_stability_plots(self):

        plot1_path = os.path.join(CACHE_DIR, "stability_plot1.png")
        plot2_path = os.path.join(CACHE_DIR, "stability_plot2.png")

        if os.path.exists(plot1_path):
            pixmap1 = QPixmap(plot1_path)
            if not pixmap1.isNull():
                self.label_plot1.setPixmap(pixmap1.scaled(
                    self.stability_plot1.width(),
                    self.stability_plot1.height(),
                    Qt.AspectRatioMode.KeepAspectRatio,
                    Qt.TransformationMode.SmoothTransformation
                ))
                print(f" StabilityPlot1 updated: {plot1_path}")
        else:
            self.label_plot1.setText("No Stability Data")

        if os.path.exists(plot2_path):
            pixmap2 = QPixmap(plot2_path)
            if not pixmap2.isNull():
                self.label_plot2.setPixmap(pixmap2.scaled(
                    self.stability_plot2.width(),
                    self.stability_plot2.height(),
                    Qt.AspectRatioMode.KeepAspectRatio,
                    Qt.TransformationMode.SmoothTransformation
                ))
                print(f" StabilityPlot2 updated: {plot2_path}")
        else:
            self.label_plot2.setText("No Stability Data")




if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec())

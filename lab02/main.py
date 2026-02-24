import ctypes
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import tkinter as tk
from tkinter import ttk
import os
import platform

lib_name = "heat_solver.so"
if platform.system() == "Windows":
    lib_name = "heat_solver.dll"
elif platform.system() == "Darwin":
    lib_name = "heat_solver.dylib"

if not os.path.exists(lib_name):
    print("Компиляция C кода...")
    compile_cmd = f"gcc -shared -O3 -fPIC -o {lib_name} heat_solver.c"
    if os.system(compile_cmd) != 0:
        print("ОШИБКА: Не удалось скомпилировать C-код. Убедитесь, что установлен GCC.")
        exit(1)

heat_lib = ctypes.CDLL(os.path.abspath(lib_name))
solve_heat_step = heat_lib.solve_heat_step
solve_heat_step.argtypes = [
    ctypes.POINTER(ctypes.c_double), # T
    ctypes.c_int,                    # N
    ctypes.c_double,                 # h
    ctypes.c_double,                 # tau
    ctypes.c_double,                 # rho
    ctypes.c_double,                 # c_heat
    ctypes.c_double,                 # lambda_coef
    ctypes.c_int                     # steps
]

class HeatSolverGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("Моделирование теплопроводности (Метод сеток)")
        
        self.is_running = False
        self.current_time = 0.0
        self.T = None
        self.x = None

        # База данных материалов: (Плотность, Теплоемкость, Теплопроводность)
        self.materials_db = {
            "Сталь": (7800.0, 460.0, 46.0),
            "Золото": (19300.0, 129.0, 317.0),
            "Медь": (8960.0, 385.0, 401.0),
            "Алюминий": (2700.0, 900.0, 237.0),
            "Серебро": (10490.0, 235.0, 429.0),
            "Платина": (21450.0, 133.0, 71.6),
            "Вольфрам": (19300.0, 134.0, 173.0),
            "Ртуть": (13534.0, 140.0, 8.3),
            "Дерево": (700.0, 2300.0, 0.15),
            "Торий": (11700.0, 113.0, 54.0),
            "Уран": (19050.0, 116.0, 27.5)
        }

        self.setup_ui()
        self.reset_simulation()

    def setup_ui(self):
        control_frame = ttk.Frame(self.root, padding="10")
        control_frame.pack(side=tk.LEFT, fill=tk.Y)

        self.params = {
            "L (Ширина пластины, м)": tk.DoubleVar(value=1.0),
            "rho (Плотность, кг/м3)": tk.DoubleVar(value=19050.0),
            "c (Теплоемкость, Дж/кгC)": tk.DoubleVar(value=116.0),
            "lambda (Теплопровод., Вт/мC)": tk.DoubleVar(value=27.5),
            "T_initial (Нач. темп., C)": tk.DoubleVar(value=20.0),
            "T_left (Лев. граница, C)": tk.DoubleVar(value=100.0),
            "T_right (Прав. граница, C)": tk.DoubleVar(value=20.0),
            "h (Шаг по пространству, м)": tk.DoubleVar(value=0.001),
            "tau (Шаг по времени, с)": tk.DoubleVar(value=0.01),
            "Ускорение времени (х)": tk.DoubleVar(value=10.0)
        }

        row = 0
        for text, var in self.params.items():
            ttk.Label(control_frame, text=text).grid(row=row, column=0, sticky=tk.W, pady=2)
            ttk.Entry(control_frame, textvariable=var, width=10).grid(row=row, column=1, pady=2)
            row += 1

        ttk.Button(control_frame, text="Применить и Сбросить", command=self.reset_simulation).grid(row=row, column=0, columnspan=2, pady=10)
        
        btn_frame = ttk.Frame(control_frame)
        btn_frame.grid(row=row+1, column=0, columnspan=2)
        ttk.Button(btn_frame, text="Старт/Пауза", command=self.toggle_simulation).pack(side=tk.LEFT, padx=5)

        # --- ОБЛАСТЬ ВЫБОРА МАТЕРИАЛОВ ---
        mat_frame = ttk.LabelFrame(control_frame, text="Быстрый выбор материала", padding="5")
        mat_frame.grid(row=row+2, column=0, columnspan=2, pady=10, sticky="ew")
        
        m_row, m_col = 0, 0
        for name, props in self.materials_db.items():
            btn = ttk.Button(mat_frame, text=name, width=10, command=lambda p=props: self.set_material(p))
            btn.grid(row=m_row, column=m_col, padx=2, pady=2)
            m_col += 1
            if m_col > 2:  # Размещаем по 3 кнопки в ряд
                m_col = 0
                m_row += 1
        # ---------------------------------

        # --- ОБЛАСТЬ ГЕНЕРАЦИИ ТАБЛИЦЫ ---
        tbl_frame = ttk.Frame(control_frame)
        tbl_frame.grid(row=row+3, column=0, columnspan=2, pady=10)
        
        self.table_time_var = tk.DoubleVar(value=120.0)
        ttk.Label(tbl_frame, text="Время для таблицы (с):").pack(side=tk.LEFT)
        ttk.Entry(tbl_frame, textvariable=self.table_time_var, width=8).pack(side=tk.LEFT, padx=5)
        
        ttk.Button(control_frame, text="Генерировать таблицу", command=self.generate_table).grid(row=row+4, column=0, columnspan=2, pady=5)
        # ---------------------------------
        
        self.time_label = ttk.Label(control_frame, text="Время: 0.00 с", font=("Arial", 12, "bold"))
        self.time_label.grid(row=row+5, column=0, columnspan=2, pady=10)

        plot_frame = ttk.Frame(self.root)
        plot_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)

        self.fig, self.ax = plt.subplots(figsize=(6, 4))
        self.line, = self.ax.plot([], [], 'r-')
        self.ax.set_xlabel("Координата x, м")
        self.ax.set_ylabel("Температура T, °C")
        self.ax.grid(True)

        self.canvas = FigureCanvasTkAgg(self.fig, master=plot_frame)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

    def set_material(self, props):
        rho, c, lam = props
        self.params["rho (Плотность, кг/м3)"].set(rho)
        self.params["c (Теплоемкость, Дж/кгC)"].set(c)
        self.params["lambda (Теплопровод., Вт/мC)"].set(lam)
        self.reset_simulation() # Автоматически применяем новые параметры

    def reset_simulation(self):
        L = self.params["L (Ширина пластины, м)"].get()
        h = self.params["h (Шаг по пространству, м)"].get()
        
        N = int(L / h) + 1
        self.x = np.linspace(0, L, N)
        self.T = np.full(N, self.params["T_initial (Нач. темп., C)"].get(), dtype=np.float64)
        
        self.T[0] = self.params["T_left (Лев. граница, C)"].get()
        self.T[-1] = self.params["T_right (Прав. граница, C)"].get()

        self.current_time = 0.0
        self.update_plot()

    def update_plot(self):
        self.line.set_data(self.x, self.T)
        self.ax.set_xlim(0, self.params["L (Ширина пластины, м)"].get())
        self.ax.set_ylim(min(self.T) - 5, max(self.T) + 5)
        self.time_label.config(text=f"Время: {self.current_time:.4f} с")
        self.canvas.draw()

    def toggle_simulation(self):
        self.is_running = not self.is_running
        if self.is_running:
            self.run_step()

    def run_step(self):
        if not self.is_running:
            return

        tau = self.params["tau (Шаг по времени, с)"].get()
        h = self.params["h (Шаг по пространству, м)"].get()
        rho = self.params["rho (Плотность, кг/м3)"].get()
        c_heat = self.params["c (Теплоемкость, Дж/кгC)"].get()
        lam = self.params["lambda (Теплопровод., Вт/мC)"].get()
        speed_multiplier = self.params["Ускорение времени (х)"].get()
        N = len(self.T)

        model_time_per_frame = 0.05 * speed_multiplier
        steps_per_frame = max(1, int(model_time_per_frame / tau)) 
        
        T_ptr = self.T.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        solve_heat_step(T_ptr, N, h, tau, rho, c_heat, lam, steps_per_frame)

        self.current_time += tau * steps_per_frame
        self.update_plot()

        self.root.after(1, self.run_step)

    def generate_table(self):
        # Берем параметры генерации прямо из интерфейса
        target_time = self.table_time_var.get()
        L = self.params["L (Ширина пластины, м)"].get()
        rho = self.params["rho (Плотность, кг/м3)"].get()
        c_heat = self.params["c (Теплоемкость, Дж/кгC)"].get()
        lam = self.params["lambda (Теплопровод., Вт/мC)"].get()
        T_init = self.params["T_initial (Нач. темп., C)"].get()
        T_left = self.params["T_left (Лев. граница, C)"].get()
        T_right = self.params["T_right (Прав. граница, C)"].get()

        print("\n--- ГЕНЕРАЦИЯ ТАБЛИЦЫ ДЛЯ ОТЧЕТА ---")
        print(f"Рассчитываем температуру в центре пластины через {target_time} секунд модельного времени...")
        print(f"Материал: rho={rho}, c={c_heat}, lambda={lam}")
        print(f"Геометрия: L={L}м, Слева={T_left}C, Справа={T_right}C, Нач={T_init}C\n")
        
        h_values = [0.1, 0.01, 0.001, 0.0001]
        tau_values = [0.1, 0.01, 0.001, 0.0001]
        
        print(f"{'h \\ tau':<10} | {'0.1':<10} | {'0.01':<10} | {'0.001':<10} | {'0.0001':<10}")
        print("-" * 65)

        for h in h_values:
            row_str = f"{h:<10} | "
            for tau in tau_values:
                N = int(L / h) + 1
                T_array = np.full(N, T_init, dtype=np.float64)
                T_array[0] = T_left
                T_array[-1] = T_right
                
                steps = int(target_time / tau)
                if steps > 0:
                    T_ptr = T_array.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
                    solve_heat_step(T_ptr, N, h, tau, rho, c_heat, lam, steps)
                
                center_idx = N // 2
                val = T_array[center_idx]
                row_str += f"{val:<10.5f} | "
                
            print(row_str)
        print("-----------------------------------------------------------------\n")

if __name__ == "__main__":
    root = tk.Tk()
    app = HeatSolverGUI(root)
    root.mainloop()
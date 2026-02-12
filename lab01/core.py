import tkinter as tk
from tkinter import ttk
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import math


class FlightSimulatorApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Моделирование полёта с анимацией")
        
        self.animation_running = False
        self.current_after_id = None
        
        self.g = 9.81
        
        self.settings_frame = ttk.LabelFrame(root, text="Параметры")
        self.settings_frame.pack(side=tk.LEFT, fill=tk.Y, padx=10, pady=10)
        
        self.create_input("Начальная скорость (м/с):", "50", "v0")
        self.create_input("Угол бросания (град):", "45", "angle")
        self.create_input("Масса тела (кг):", "1.0", "m")
        self.create_input("Коэфф. сопротивления (k):", "0.02", "k")
        self.create_input("Шаг моделирования (dt, с):", "0.1", "dt")
        
        self.btn_calc = ttk.Button(self.settings_frame, text="Запуск с анимацией", command=self.calculate)
        self.btn_calc.pack(pady=5, fill=tk.X)
        
        self.btn_clear = ttk.Button(self.settings_frame, text="Очистить график", command=self.clear_plot)
        self.btn_clear.pack(pady=5, fill=tk.X)
        
        self.result_text = tk.Text(self.settings_frame, height=15, width=35, font=("Consolas", 9))
        self.result_text.pack(pady=10)
        
        self.figure, self.ax = plt.subplots(figsize=(6, 5))
        self.ax.set_title("Траектория полёта")
        self.ax.set_xlabel("Дальность (x), м")
        self.ax.set_ylabel("Высота (y), м")
        self.ax.grid(True)
        
        self.canvas = FigureCanvasTkAgg(self.figure, master=root)
        self.canvas.get_tk_widget().pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)

    def create_input(self, label_text, default_val, var_name):
        frame = ttk.Frame(self.settings_frame)
        frame.pack(fill=tk.X, pady=2)
        ttk.Label(frame, text=label_text).pack(side=tk.LEFT)
        entry = ttk.Entry(frame, width=10)
        entry.insert(0, default_val)
        entry.pack(side=tk.RIGHT)
        setattr(self, f"entry_{var_name}", entry)

    def clear_plot(self):
        self.animation_running = False
        if self.current_after_id:
            self.root.after_cancel(self.current_after_id)
            
        self.ax.clear()
        self.ax.set_title("Траектория полёта")
        self.ax.set_xlabel("Дальность (x), м")
        self.ax.set_ylabel("Высота (y), м")
        self.ax.grid(True)
        self.result_text.delete(1.0, tk.END)
        self.canvas.draw()

    def calculate(self):
        try:
            v0 = float(self.entry_v0.get())
            angle = float(self.entry_angle.get())
            m = float(self.entry_m.get())
            k = float(self.entry_k.get())
            dt = float(self.entry_dt.get())
        except ValueError:
            self.result_text.insert(tk.END, "Ошибка ввода!\n")
            return

        x, y = 0.0, 0.0
        alpha_rad = math.radians(angle)
        vx = v0 * math.cos(alpha_rad)
        vy = v0 * math.sin(alpha_rad)
        
        xs, ys = [x], [y]
        max_h = 0.0
        
        while y >= 0:
            v = math.sqrt(vx**2 + vy**2)
            ax_val = -(k / m) * v * vx
            ay_val = -self.g - (k / m) * v * vy
            
            x += vx * dt
            y += vy * dt
            vx += ax_val * dt
            vy += ay_val * dt
            
            if y >= 0:
                xs.append(x)
                ys.append(y)
                if y > max_h:
                    max_h = y
        
        final_v = math.sqrt(vx**2 + vy**2)
        final_x = xs[-1]
        
        res = f"--- dt = {dt} с ---\n"
        res += f"L: {final_x:.2f} м | H: {max_h:.2f} м\n"
        res += f"V_k: {final_v:.2f} м/с\n\n"
        self.result_text.insert(tk.END, res)
        self.result_text.see(tk.END)

        label_str = f"dt={dt}"
        self.run_animation(xs, ys, label_str)

    def run_animation(self, xs, ys, label):
        line, = self.ax.plot([], [], label=label, linewidth=2)
        self.ax.legend()

        current_xlim = self.ax.get_xlim()
        current_ylim = self.ax.get_ylim()
        
        needed_max_x = max(xs) * 1.1
        needed_max_y = max(ys) * 1.1

        new_xlim_right = max(current_xlim[1], needed_max_x)
        new_ylim_top = max(current_ylim[1], needed_max_y)

        if current_xlim[1] <= 1.0 and current_ylim[1] <= 1.0:
            new_xlim_right = needed_max_x
            new_ylim_top = needed_max_y

        self.ax.set_xlim(0, new_xlim_right)
        self.ax.set_ylim(0, new_ylim_top)
        self.canvas.draw()

        total_points = len(xs)
        step = max(1, total_points // 100) 
        
        self.animation_running = True
        current_idx = 0

        def animate_step():
            nonlocal current_idx
            if not self.animation_running:
                return

            current_idx += step
            
            if current_idx >= total_points:
                line.set_data(xs, ys)
                self.canvas.draw()
                self.animation_running = False
                return

            line.set_data(xs[:current_idx], ys[:current_idx])
            self.canvas.draw()
            
            self.current_after_id = self.root.after(20, animate_step)

        animate_step()


if __name__ == "__main__":
    root = tk.Tk()
    app = FlightSimulatorApp(root)
    root.mainloop()

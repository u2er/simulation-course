import tkinter as tk
from tkinter import ttk, messagebox
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg


UI_BG = "#2c3e50"
PLOT_BG = "#ffffff"
ACCENT_COLOR = "#3498db"


class WeatherModelGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("Markov Weather Model")
        self.root.geometry("1200x800")
        self.root.configure(bg=UI_BG)
        
        self.states = ['Ясно', 'Облачно', 'Пасмурно']
        self.is_running = False
        
        self.setup_layout()

    def setup_layout(self):
        # Левая панель управления
        self.side_panel = tk.Frame(self.root, bg=UI_BG, width=300, padx=20, pady=20)
        self.side_panel.pack(side=tk.LEFT, fill=tk.Y)
        self.side_panel.pack_propagate(False)

        tk.Label(self.side_panel, text="ПАРАМЕТРЫ", fg="white", bg=UI_BG, 
                 font=("Arial", 14, "bold")).pack(pady=(0, 20))

        matrix_frame = tk.Frame(self.side_panel, bg=UI_BG)
        matrix_frame.pack(fill=tk.X, pady=10)
        
        tk.Label(matrix_frame, text="Интенсивности переходов:", fg="#bdc3c7", bg=UI_BG, font=("Arial", 10)).grid(row=0, columnspan=4, pady=5)
        
        self.q_entries = {}
        defaults = {(0,1): "0.3", (0,2): "0.2", (1,0): "0.2", (1,2): "0.4", (2,0): "0.1", (2,1): "0.4"}
        
        for i in range(3):
            tk.Label(matrix_frame, text=self.states[i][:5], fg="white", bg=UI_BG, width=7).grid(row=i+1, column=0, pady=5)
            for j in range(3):
                if i == j:
                    tk.Label(matrix_frame, text="---", fg="#7f8c8d", bg=UI_BG).grid(row=i+1, column=j+1)
                else:
                    ent = tk.Entry(matrix_frame, width=5, justify='center', bg="#34495e", fg="white", insertbackground="white", borderwidth=0)
                    ent.insert(0, defaults[(i,j)])
                    ent.grid(row=i+1, column=j+1, padx=2, pady=2)
                    self.q_entries[(i,j)] = ent

        tk.Label(self.side_panel, text="Симуляция:", fg="#bdc3c7", bg=UI_BG, font=("Arial", 10)).pack(pady=(20, 5), anchor="w")
        
        self.t_max_var = tk.StringVar(value="30")
        self.delay_var = tk.StringVar(value="50")

        self.create_input(self.side_panel, "Дней всего:", self.t_max_var)
        self.create_input(self.side_panel, "Задержка (мс):", self.delay_var)

        self.start_btn = tk.Button(self.side_panel, text="ЗАПУСТИТЬ", bg="#27ae60", fg="white", 
                                   font=("Arial", 10, "bold"), command=self.start, borderwidth=0, cursor="hand2")
        self.start_btn.pack(fill=tk.X, pady=(30, 10), ipady=10)

        self.stop_btn = tk.Button(self.side_panel, text="ОСТАНОВИТЬ", bg="#c0392b", fg="white", 
                                  font=("Arial", 10, "bold"), command=self.stop, state=tk.DISABLED, borderwidth=0)
        self.stop_btn.pack(fill=tk.X, ipady=10)

        self.log = tk.Text(self.side_panel, height=8, bg="#1a252f", fg="#ecf0f1", font=("Courier", 9), borderwidth=0, padx=10, pady=10)
        self.log.pack(fill=tk.X, pady=(20, 0))

        self.main_panel = tk.Frame(self.root, bg="white")
        self.main_panel.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True, padx=10, pady=10)

        self.fig, (self.ax_sim, self.ax_stat) = plt.subplots(2, 1, figsize=(8, 8), constrained_layout=True)
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.main_panel)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        
        self.init_plots()

    def create_input(self, parent, label, var):
        frame = tk.Frame(parent, bg=UI_BG)
        frame.pack(fill=tk.X, pady=2)
        tk.Label(frame, text=label, fg="white", bg=UI_BG, font=("Arial", 9)).pack(side=tk.LEFT)
        tk.Entry(frame, textvariable=var, width=8, justify='right', bg="#34495e", fg="white", borderwidth=0).pack(side=tk.RIGHT)

    def init_plots(self):
        self.ax_sim.set_title("ДИНАМИКА СОСТОЯНИЙ ПОГОДЫ", fontsize=12, pad=10)
        self.ax_sim.set_yticks([0, 1, 2])
        self.ax_sim.set_yticklabels(self.states)
        self.ax_sim.set_xlabel("Время (дни)")
        self.ax_sim.grid(True, alpha=0.3, linestyle='--')

        self.ax_stat.set_title("СРАВНЕНИЕ РАСПРЕДЕЛЕНИЙ", fontsize=12, pad=10)
        self.ax_stat.set_ylabel("Вероятность")
        self.ax_stat.grid(True, axis='y', alpha=0.3)
        
        self.canvas.draw()

    def get_q(self):
        Q = np.zeros((3, 3))
        try:
            for (i, j), ent in self.q_entries.items():
                Q[i, j] = float(ent.get())
            for i in range(3):
                Q[i, i] = -np.sum(Q[i, :])
            return Q
        except:
            messagebox.showerror("Error", "Введите числа!")
            return None

    def start(self):
        self.Q = self.get_q()
        if self.Q is None: return
        
        A = self.Q.T.copy(); A[-1, :] = 1
        self.pi_theo = np.linalg.solve(A, [0, 0, 1])
        
        self.t_max = float(self.t_max_var.get())
        self.delay = int(self.delay_var.get())
        
        self.t, self.curr_s = 0, 0
        self.history_t, self.history_s = [0], [0]
        self.durations = np.zeros(3)
        
        self.is_running = True
        self.start_btn.config(state=tk.DISABLED, bg="#7f8c8d")
        self.stop_btn.config(state=tk.NORMAL)
        self.log.delete(1.0, tk.END)
        self.loop()

    def stop(self):
        self.is_running = False
        self.start_btn.config(state=tk.NORMAL, bg="#27ae60")
        self.stop_btn.config(state=tk.DISABLED)
        self.show_final_stats()

    def loop(self):
        if not self.is_running or self.t >= self.t_max:
            self.stop()
            return
            
        rate = -self.Q[self.curr_s, self.curr_s]
        dt = np.random.exponential(1/rate) if rate > 0 else self.t_max
        
        actual_dt = min(dt, self.t_max - self.t)
        self.durations[self.curr_s] += actual_dt
        self.t += actual_dt
        
        if self.t < self.t_max:
            probs = self.Q[self.curr_s].copy()
            probs[self.curr_s] = 0
            self.curr_s = np.random.choice([0,1,2], p=probs/rate)
            
        self.history_t.append(self.t)
        self.history_s.append(self.curr_s)
        
        self.update_plots()
        self.root.after(self.delay, self.loop)

    def update_plots(self):
        self.ax_sim.clear()
        self.ax_sim.step(self.history_t, self.history_s, where='post', color=ACCENT_COLOR, lw=2)
        self.ax_sim.fill_between(self.history_t, self.history_s, step="post", alpha=0.2, color=ACCENT_COLOR)
        self.ax_sim.set_yticks([0, 1, 2])
        self.ax_sim.set_yticklabels(self.states)
        self.ax_sim.set_xlim(0, self.t_max)
        self.ax_sim.set_ylim(-0.5, 2.5)
        self.ax_sim.set_title(f"ТЕКУЩИЙ ДЕНЬ: {self.t:.1f}", fontsize=11)
        self.ax_sim.grid(True, alpha=0.3)
        self.canvas.draw()

    def show_final_stats(self):
        pi_emp = self.durations / self.t_max
        
        self.ax_stat.clear()
        x = np.arange(3)
        w = 0.35
        
        b1 = self.ax_stat.bar(x - w/2, self.pi_theo, w, label='Теория', color='#2ecc71')
        b2 = self.ax_stat.bar(x + w/2, pi_emp, w, label='Эмпирика', color='#e67e22')
        
        self.ax_stat.set_xticks(x)
        self.ax_stat.set_xticklabels(self.states)
        self.ax_stat.legend()
        self.ax_stat.set_title("ИТОГОВОЕ РАСПРЕДЕЛЕНИЕ", fontsize=11)

        for bar in b1+b2:
            h = bar.get_height()
            self.ax_stat.text(bar.get_x()+bar.get_width()/2, h, f'{h:.2f}', ha='center', va='bottom', fontsize=9)
        
        self.canvas.draw()
        
        res = "РЕЗУЛЬТАТЫ:\n"
        for i in range(3):
            res += f"{self.states[i][:5]}: T:{self.pi_theo[i]:.2f} E:{pi_emp[i]:.2f}\n"
        self.log.insert(tk.END, res)

if __name__ == "__main__":
    root = tk.Tk()
    app = WeatherModelGUI(root)
    root.mainloop()
import customtkinter as ctk
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import tkinter.messagebox as messagebox
import math


ctk.set_appearance_mode("Dark")
ctk.set_default_color_theme("blue")

plt.rcParams.update({
    "figure.facecolor": "#212121",
    "axes.facecolor": "#212121",
    "axes.edgecolor": "#444444",
    "text.color": "#E0E0E0",
    "axes.labelcolor": "#E0E0E0",
    "xtick.color": "#E0E0E0",
    "ytick.color": "#E0E0E0",
    "grid.color": "#333333",
    "grid.linestyle": "--",
    "axes.spines.top": False,
    "axes.spines.right": False,
    "font.family": "sans-serif",
    "font.size": 10
})


class SimulationApp(ctk.CTk):
    def __init__(self):
        super().__init__()

        self.title("Имитационное моделирование случайных величин")
        self.geometry("1200x850")
        self.minsize(1000, 750)

        self.tabview = ctk.CTkTabview(self)
        self.tabview.pack(fill="both", expand=True, padx=15, pady=15)

        self.tab_discrete = self.tabview.add("Дискретная СВ")
        self.tab_normal = self.tabview.add("Нормальная СВ")

        self.setup_discrete_tab()
        self.setup_normal_tab()

    # ==========================================
    # ДИСКРЕТНАЯ СЛУЧАЙНАЯ ВЕЛИЧИНА
    # ==========================================
    def setup_discrete_tab(self):
        left_panel = ctk.CTkFrame(self.tab_discrete, width=300)
        left_panel.pack(side="left", fill="y", padx=10, pady=10)
        
        right_panel = ctk.CTkFrame(self.tab_discrete)
        right_panel.pack(side="right", fill="both", expand=True, padx=10, pady=10)

        ctk.CTkLabel(left_panel, text="Значения X и P:", font=("Helvetica", 16, "bold")).pack(pady=(10, 5))
        
        self.scroll_frame = ctk.CTkScrollableFrame(left_panel, width=250, height=250)
        self.scroll_frame.pack(fill="x", padx=10, pady=5)

        self.discrete_rows = []
        
        initial_data = [(1, 0.1), (2, 0.2), (3, 0.4), (4, 0.2), (5, 0.1)]
        for x, p in initial_data:
            self.add_discrete_row(str(x), str(p))

        btn_frame = ctk.CTkFrame(left_panel, fg_color="transparent")
        btn_frame.pack(fill="x", padx=10, pady=5)
        ctk.CTkButton(btn_frame, text="+ Добавить", width=110, command=self.add_discrete_row, fg_color="#2E7D32", hover_color="#1B5E20").pack(side="left", padx=5)
        ctk.CTkButton(btn_frame, text="- Удалить", width=110, command=self.remove_discrete_row, fg_color="#C62828", hover_color="#8E0000").pack(side="right", padx=5)

        ctk.CTkLabel(left_panel, text="Объемы выборок N (через запятую):").pack(pady=(20, 0))
        self.entry_n_discrete = ctk.CTkEntry(left_panel, width=250)
        self.entry_n_discrete.insert(0, "10, 100, 1000, 10000")
        self.entry_n_discrete.pack(pady=5)

        ctk.CTkButton(left_panel, text="Моделировать", height=40, font=("Helvetica", 14, "bold"), command=self.run_discrete_simulation).pack(pady=20, fill="x", padx=10)

        # --- Правая панель: Графики и статистика ---
        self.fig_disc = plt.Figure(figsize=(8, 5))
        self.canvas_disc = FigureCanvasTkAgg(self.fig_disc, master=right_panel)
        self.canvas_disc.get_tk_widget().pack(fill="both", expand=True)

        self.textbox_discrete = ctk.CTkTextbox(right_panel, height=150, font=("Menlo", 13))
        self.textbox_discrete.pack(fill="x", pady=(10, 0))

    def add_discrete_row(self, x_val="", p_val=""):
        row_frame = ctk.CTkFrame(self.scroll_frame, fg_color="transparent")
        row_frame.pack(fill="x", pady=2)
        
        e_x = ctk.CTkEntry(row_frame, width=100, placeholder_text="X")
        e_x.insert(0, x_val)
        e_x.pack(side="left", padx=5)
        
        e_p = ctk.CTkEntry(row_frame, width=100, placeholder_text="P")
        e_p.insert(0, p_val)
        e_p.pack(side="left", padx=5)
        
        self.discrete_rows.append((row_frame, e_x, e_p))

    def remove_discrete_row(self):
        if len(self.discrete_rows) > 2:
            frame, e_x, e_p = self.discrete_rows.pop()
            frame.destroy()
        else:
            messagebox.showwarning("Внимание", "Должно быть минимум 2 значения!")

    def run_discrete_simulation(self):
        self.textbox_discrete.delete("1.0", "end")
        self.fig_disc.clear()

        x_vals, p_vals = [], []
        try:
            for _, e_x, e_p in self.discrete_rows:
                x_vals.append(float(e_x.get().strip()))
                p_vals.append(float(e_p.get().strip()))
            
            n_list = [int(n.strip()) for n in self.entry_n_discrete.get().split(',')]

        except ValueError:
            messagebox.showerror("Ошибка", "Проверьте правильность ввода: X, P и N должны быть числами.")
            return

        x_vals = np.array(x_vals)
        p_vals = np.array(p_vals)

        sum_p = np.sum(p_vals)
        if not np.isclose(sum_p, 1.0, atol=1e-5):
            p_vals = p_vals / sum_p
            for i, (_, _, e_p) in enumerate(self.discrete_rows):
                e_p.delete(0, "end")
                e_p.insert(0, f"{p_vals[i]:.4f}")
            messagebox.showwarning("Нормировка", f"Сумма вероятностей была {sum_p:.4f}. Произведена автоматическая нормировка до 1.0.")

        theo_mean = np.sum(x_vals * p_vals)
        theo_var = np.sum((x_vals ** 2) * p_vals) - theo_mean ** 2

        cols = 2
        rows = math.ceil(len(n_list) / cols)
        axs = self.fig_disc.subplots(rows, cols)
        if rows * cols == 1: axs = [axs]
        else: axs = axs.flatten()

        alpha_level = 0.05
        bar_width = 0.35
        x_indices = np.arange(len(x_vals))

        text_output = f"ТЕОРЕТИЧЕСКИ: E = {theo_mean:.4f}, D = {theo_var:.4f}\n" + "-"*70 + "\n"

        for i, N in enumerate(n_list):
            if i >= len(axs): break
            ax = axs[i]
    
            u = np.random.rand(N)
            cdf = np.cumsum(p_vals)
            indices = np.searchsorted(cdf, u)
            sample = x_vals[indices]
            
            emp_mean = np.mean(sample)
            emp_var = np.var(sample, ddof=0)
            
            err_mean = abs(emp_mean - theo_mean) / abs(theo_mean) * 100 if theo_mean != 0 else 0
            err_var = abs(emp_var - theo_var) / abs(theo_var) * 100 if theo_var != 0 else 0

            unique, counts = np.unique(sample, return_counts=True)
            emp_freq = dict(zip(unique, counts))
            observed = np.array([emp_freq.get(x, 0) for x in x_vals])
            emp_p = observed / N
            expected = N * p_vals
            
            chi2_stat = np.sum((observed - expected)**2 / expected)
            chi2_crit = stats.chi2.ppf(1 - alpha_level, len(x_vals) - 1)
            is_accepted = "ПРИНЯТА" if chi2_stat < chi2_crit else "ОТВЕРГНУТА"

            ax.bar(x_indices - bar_width/2, p_vals, bar_width, label='Теория', color='#FF4081', alpha=0.8, edgecolor='#C51162')
            ax.bar(x_indices + bar_width/2, emp_p, bar_width, label='Эмпирика', color='#00E5FF', alpha=0.8, edgecolor='#00B8D4')
            
            ax.set_title(f"N = {N}", color="#FFFFFF", pad=10)
            ax.set_xticks(x_indices)
            ax.set_xticklabels([str(x) for x in x_vals])
            ax.grid(True, alpha=0.3)
            
            if i == 0:
                ax.legend(frameon=False)

            text_output += (f"N={N:<6} | Ср: {emp_mean:.4f} (Погр: {err_mean:>5.2f}%) | "
                            f"Дисп: {emp_var:.4f} (Погр: {err_var:>5.2f}%) | "
                            f"χ²: {chi2_stat:.2f} (Крит: {chi2_crit:.2f}) -> {is_accepted}\n")

        self.fig_disc.tight_layout()
        self.canvas_disc.draw()
        
        self.textbox_discrete.insert("end", text_output)

    # ==========================================
    # НОРМАЛЬНАЯ СЛУЧАЙНАЯ ВЕЛИЧИНА
    # ==========================================
    def setup_normal_tab(self):
        left_panel = ctk.CTkFrame(self.tab_normal, width=300)
        left_panel.pack(side="left", fill="y", padx=10, pady=10)
        
        right_panel = ctk.CTkFrame(self.tab_normal)
        right_panel.pack(side="right", fill="both", expand=True, padx=10, pady=10)

        ctk.CTkLabel(left_panel, text="Параметры распределения:", font=("Helvetica", 16, "bold")).pack(pady=(10, 15))

        ctk.CTkLabel(left_panel, text="Мат. ожидание (μ):").pack()
        self.entry_mu = ctk.CTkEntry(left_panel, width=200)
        self.entry_mu.insert(0, "0")
        self.entry_mu.pack(pady=(0, 10))

        ctk.CTkLabel(left_panel, text="Ср. кв. отклонение (σ):").pack()
        self.entry_sigma = ctk.CTkEntry(left_panel, width=200)
        self.entry_sigma.insert(0, "1")
        self.entry_sigma.pack(pady=(0, 10))

        ctk.CTkLabel(left_panel, text="Объемы выборок N (через запятую):").pack()
        self.entry_n_normal = ctk.CTkEntry(left_panel, width=200)
        self.entry_n_normal.insert(0, "10, 100, 1000, 10000")
        self.entry_n_normal.pack(pady=(0, 20))

        ctk.CTkButton(left_panel, text="Построить графики", height=40, font=("Helvetica", 14, "bold"), command=self.run_normal_simulation).pack(fill="x", padx=10)

        self.fig_norm = plt.Figure(figsize=(8, 6))
        self.canvas_norm = FigureCanvasTkAgg(self.fig_norm, master=right_panel)
        self.canvas_norm.get_tk_widget().pack(fill="both", expand=True)

    def run_normal_simulation(self):
        self.fig_norm.clear()
        
        try:
            mu = float(self.entry_mu.get())
            sigma = float(self.entry_sigma.get())
            n_list = [int(n.strip()) for n in self.entry_n_normal.get().split(',')]
        except ValueError:
            messagebox.showerror("Ошибка", "Введите корректные числовые значения.")
            return
            
        if sigma <= 0:
            messagebox.showerror("Ошибка", "σ должно быть строго больше 0.")
            return

        cols = 2
        rows = math.ceil(len(n_list) / cols)
        axs = self.fig_norm.subplots(rows, cols)
        if rows * cols == 1: axs = [axs]
        else: axs = axs.flatten()

        for i, N in enumerate(n_list):
            if i >= len(axs): break
            ax = axs[i]
            
            half_n = math.ceil(N / 2)
            u1 = u2 = np.random.rand(half_n)

            z0 = np.sqrt(-2.0 * np.log(u1)) * np.cos(2.0 * np.pi * u2)
            z1 = np.sqrt(-2.0 * np.log(u1)) * np.sin(2.0 * np.pi * u2)
 
            z = np.concatenate((z0, z1))[:N]
            sample = mu + sigma * z
            
            ax.hist(sample, bins='auto', density=True, color='#B388FF', edgecolor='#651FFF', alpha=0.7, linewidth=1.2)
            
            x_range = np.linspace(mu - 4*sigma, mu + 4*sigma, 1000)
            ax.plot(x_range, stats.norm.pdf(x_range, mu, sigma), linewidth=2.5, color='#00E676')
            
            ax.set_title(f"N = {N}", color="#FFFFFF", pad=10)
            ax.grid(True, alpha=0.2)

        self.fig_norm.tight_layout()
        self.canvas_norm.draw()


if __name__ == "__main__":
    app = SimulationApp()
    app.mainloop()
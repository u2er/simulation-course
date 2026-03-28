import time
import random as py_random 
import customtkinter as ctk
import tkinter as tk


class LinearCongruentialGenerator:
    def __init__(self, seed, a=16807, c=12345, m=(2 ** 31 - 1)):
        self.state = seed
        self.a = a
        self.c = c
        self.m = m

    def random(self):
        self.state = (self.a * self.state + self.c) % self.m
        return self.state / self.m


class MagicApp(ctk.CTk):
    def __init__(self):
        super().__init__()

        self.rng = LinearCongruentialGenerator(seed=int(time.time() * 1000))
        
        self.responses = [
            "Бесспорно", "Предрешено", "Без сомнений", "Определённо да", 
            "Пока неясно", "Спроси позже", "Даже не думай", "Никогда", "Нет"
        ]
        self.probs = [1.0 / len(self.responses)] * len(self.responses)

        self.title("Magic 8-Ball Pro")
        self.geometry("600x750")
        self.minsize(550, 650)
        
        ctk.set_appearance_mode("system")
        ctk.set_default_color_theme("blue")

        self.tabview = ctk.CTkTabview(self, corner_radius=15)
        self.tabview.pack(padx=20, pady=20, expand=True, fill="both")

        self.tabview.add("Да / Нет")
        self.tabview.add("Шар предсказаний")
        self.tabview.add("Настройки")

        self.setup_yes_no_tab()
        self.setup_magic_ball_tab()
        self.setup_settings_tab()

    def setup_yes_no_tab(self):
        tab = self.tabview.tab("Да / Нет")
        tab.grid_rowconfigure((0, 4), weight=1)
        tab.grid_columnconfigure(0, weight=1)

        ctk.CTkLabel(tab, text="Задавайте свой вопрос", font=("Helvetica", 24, "bold")).grid(row=1, column=0, pady=10)

        self.yes_no_result = ctk.CTkLabel(tab, text="❓", font=("Helvetica", 100, "bold"))
        self.yes_no_result.grid(row=2, column=0, pady=40)
        
        ctk.CTkButton(tab, text="ПОЛУЧИТЬ ОТВЕТ", font=("Helvetica", 16, "bold"), 
                      height=50, width=200, command=self.generate_yes_no).grid(row=3, column=0)

    def generate_yes_no(self):
        val = self.rng.random()
        res, color = ("ДА", "#2ecc71") if val >= 0.5 else ("НЕТ", "#e74c3c")
        self.yes_no_result.configure(text=res, text_color=color)

    def setup_magic_ball_tab(self):
        tab = self.tabview.tab("Шар предсказаний")
        tab.grid_rowconfigure((0, 4), weight=1)
        tab.grid_columnconfigure(0, weight=1)

        self.canvas = tk.Canvas(tab, width=300, height=300, highlightthickness=0)
        bg_color = self._apply_appearance_mode(ctk.ThemeManager.theme["CTkFrame"]["fg_color"])
        self.canvas.configure(bg=bg_color)
        self.canvas.grid(row=1, column=0)

        self.shake_btn = ctk.CTkButton(tab, text="ПОТРЯСТИ ШАР", font=("Helvetica", 16, "bold"), 
                                       fg_color="#5a2e98", hover_color="#451e7a",
                                       height=50, width=200, command=self.start_shake_animation)
        self.shake_btn.grid(row=2, column=0, pady=30)
        self.draw_resting_ball()

    def draw_resting_ball(self):
        self.canvas.delete("all")
        self.canvas.create_oval(20, 20, 280, 280, fill="#151515", outline="#222", width=2, tags="ball")
        self.canvas.create_oval(90, 90, 210, 210, fill="#eeeeee", outline="", tags="ball")
        self.canvas.create_text(150, 150, text="8", font=("Helvetica", 70, "bold"), fill="#151515", tags="ball")

    def draw_answer_ball(self, text):
        self.canvas.delete("all")
        self.canvas.create_oval(20, 20, 280, 280, fill="#151515", outline="#222", width=2)
        self.canvas.create_oval(60, 60, 240, 240, fill="#0a0a0a", outline="#111", width=3)
        self.canvas.create_polygon(150, 70, 65, 205, 235, 205, fill="#1c3b7a", outline="#2a52a1", width=2)
        
        wrapped_text = self.wrap_text(text, 12)
        self.canvas.create_text(150, 155, text=wrapped_text.upper(), font=("Helvetica", 10, "bold"), 
                                fill="#61a0ff", justify="center", width=110)

    def wrap_text(self, text, max_chars):
        words = text.split()
        lines = []
        current_line = []
        for w in words:
            if len(" ".join(current_line + [w])) <= max_chars:
                current_line.append(w)
            else:
                lines.append(" ".join(current_line))
                current_line = [w]
        lines.append(" ".join(current_line))
        return "\n".join(lines)

    def start_shake_animation(self):
        self.shake_btn.configure(state="disabled")
        self.draw_resting_ball()
        self.shake_count = 0
        self.animate_shake()

    def animate_shake(self):
        if self.shake_count < 15:
            dx, dy = py_random.randint(-12, 12), py_random.randint(-12, 12)
            self.canvas.move("ball", dx, dy)
            self.shake_count += 1
            self.after(30, lambda: [self.canvas.move("ball", -dx, -dy), self.animate_shake()])
        else:
            self.reveal_answer()

    def reveal_answer(self):
        val = self.rng.random()
        cumulative = 0
        response = self.responses[-1]
        for i, p in enumerate(self.probs):
            cumulative += p
            if val <= cumulative:
                response = self.responses[i]
                break
        self.draw_answer_ball(response)
        self.shake_btn.configure(state="normal")

    def setup_settings_tab(self):
        tab = self.tabview.tab("Настройки")
        
        add_frame = ctk.CTkFrame(tab, fg_color="transparent")
        add_frame.pack(fill="x", padx=20, pady=10)
        
        self.new_opt_entry = ctk.CTkEntry(add_frame, placeholder_text="Новый вариант ответа...")
        self.new_opt_entry.pack(side="left", expand=True, fill="x", padx=(0, 10))
        
        add_btn = ctk.CTkButton(add_frame, text="Добавить", width=100, command=self.add_response)
        add_btn.pack(side="right")

        self.scroll_frame = ctk.CTkScrollableFrame(tab, label_text="Список ответов и их вероятности")
        self.scroll_frame.pack(padx=20, pady=10, fill="both", expand=True)

        self.status_label = ctk.CTkLabel(tab, text="Сумма вероятностей должна быть 1.0", font=("Helvetica", 12))
        self.status_label.pack()

        apply_btn = ctk.CTkButton(tab, text="ПРИМЕНИТЬ И НОРМИРОВАТЬ", fg_color="#2ecc71", hover_color="#27ae60",
                                  height=40, command=self.apply_probabilities)
        apply_btn.pack(pady=15)

        self.refresh_settings_list()

    def refresh_settings_list(self):
        for widget in self.scroll_frame.winfo_children():
            widget.destroy()

        self.entry_widgets = []
        for i, (text, prob) in enumerate(zip(self.responses, self.probs)):
            row = ctk.CTkFrame(self.scroll_frame, fg_color="transparent")
            row.pack(fill="x", pady=2)
            
            lbl = ctk.CTkLabel(row, text=text, anchor="w", width=180)
            lbl.pack(side="left", padx=5)
            
            del_btn = ctk.CTkButton(row, text="✕", width=30, fg_color="#e74c3c", hover_color="#c0392b",
                                    command=lambda idx=i: self.delete_response(idx))
            del_btn.pack(side="right", padx=5)

            ent = ctk.CTkEntry(row, width=80)
            ent.insert(0, f"{prob:.3f}")
            ent.pack(side="right", padx=5)
            self.entry_widgets.append(ent)

    def add_response(self):
        text = self.new_opt_entry.get().strip()
        if text:
            self.responses.append(text)
            self.probs.append(0.0)
            self.new_opt_entry.delete(0, "end")
            self.refresh_settings_list()

    def delete_response(self, index):
        if len(self.responses) > 1:
            self.responses.pop(index)
            self.probs.pop(index)
            self.refresh_settings_list()

    def apply_probabilities(self):
        new_raw_probs = []
        for ent in self.entry_widgets:
            try:
                val = float(ent.get().replace(",", "."))
                new_raw_probs.append(max(0, val))
            except Exception:
                new_raw_probs.append(0.0)
        
        total = sum(new_raw_probs)
        if total == 0:
            new_raw_probs = [1.0] * len(new_raw_probs)
            total = sum(new_raw_probs)
        
        self.probs = [v / total for v in new_raw_probs]
        self.refresh_settings_list()
        
        if abs(total - 1.0) > 0.01:
            self.status_label.configure(text=f"Нормировано (было {total:.2f})", text_color="#e67e22")
        else:
            self.status_label.configure(text="Применено успешно!", text_color="#2ecc71")


if __name__ == "__main__":
    app = MagicApp()
    app.mainloop()

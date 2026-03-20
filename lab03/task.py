import random
import math
import tkinter as tk
from tkinter import messagebox


class Cell:
    wind_x = 1.0  
    wind_y = 0.0

    def __init__(self, grow_p, fire_p, cloud_p, dead_p, rain_p, rain_end_p, scat_p):
        self.update_params(grow_p, fire_p, cloud_p, dead_p, rain_p, rain_end_p, scat_p)
        self.status = "tree"       
        self.substatus = "nothing" 
        self.neighbors = {}
        
        self.next_status = self.status
        self.next_substatus = self.substatus
        
        self.rect_ids = [None, None, None]
        self.current_colors = [None, None, None]

    def update_params(self, grow_p, fire_p, cloud_p, dead_p, rain_p, rain_end_p, scat_p):
        self.grow_p = grow_p
        self.fire_p = fire_p
        self.cloud_p = cloud_p
        self.dead_p = dead_p
        self.rain_p = rain_p
        self.rain_end_p = rain_end_p
        self.scat_p = scat_p

    def set_neighbors(self, tl=None, tc=None, tr=None, ml=None, mr=None, bl=None, bc=None, br=None):
        self.neighbors = dict(tl=tl, tc=tc, tr=tr, ml=ml, mr=mr, bl=bl, bc=bc, br=br)

    def tick(self):
        fire_influence = 0.0
        cloud_influence = 0.0
        rain_influence = 0.0

        dirs = {
            'tl': (-1, -1), 'tc': (0, -1), 'tr': (1, -1),
            'ml': (-1,  0),                'mr': (1,  0),
            'bl': (-1,  1), 'bc': (0,  1), 'br': (1,  1)
        }

        for key, neighbor in self.neighbors.items():
            if neighbor is None:
                continue
            
            nx, ny = dirs[key]
            dist = math.hypot(nx, ny)
            dir_x, dir_y = -nx / dist, -ny / dist
            
            dot_product = (Cell.wind_x * dir_x) + (Cell.wind_y * dir_y)
            base_weight = 1.0 / dist
            weight = max(0.0, base_weight + dot_product)

            if neighbor.status == "fire":
                fire_influence += weight
            if neighbor.substatus == "cloud":
                cloud_influence += weight
            if neighbor.substatus == "rain":
                rain_influence += weight

        self.next_substatus = self.substatus
        
        if self.substatus == "nothing":
            prob = self.cloud_p + (cloud_influence * 0.2)
            if random.random() < min(1.0, prob):
                self.next_substatus = "cloud"
        elif self.substatus == "cloud":
            prob_rain = self.rain_p + (rain_influence * 0.3)
            if random.random() < min(1.0, prob_rain):
                self.next_substatus = "rain"
            elif random.random() < self.scat_p:
                self.next_substatus = "nothing"
        elif self.substatus == "rain":
            if random.random() < self.rain_end_p:
                self.next_substatus = "nothing"

        self.next_status = self.status
        
        if self.status in ["dead", "empty"]:
            if random.random() < self.grow_p:
                self.next_status = "tree"
        elif self.status == "tree":
            if self.next_substatus != "rain":
                prob_fire = self.fire_p + (fire_influence * 0.25)
                if random.random() < min(1.0, prob_fire):
                    self.next_status = "fire"
        elif self.status == "fire":
            if self.next_substatus == "rain":
                self.next_status = "dead"
            elif random.random() < self.dead_p:
                self.next_status = "dead"

    def apply_tick(self):
        self.status = self.next_status
        self.substatus = self.next_substatus


class SimulationApp:
    COLORS_COMBINED = {
        ("tree", "nothing"): "#2ca02c",   
        ("tree", "cloud"): "#98df8a",     
        ("tree", "rain"): "#1f77b4",      
        ("fire", "nothing"): "#d62728",   
        ("fire", "cloud"): "#ff9896",     
        ("fire", "rain"): "#9467bd",      
        ("dead", "nothing"): "#555555",   
        ("dead", "cloud"): "#7f7f7f",     
        ("dead", "rain"): "#8c564b",      
        ("empty", "nothing"): "#555555",  
        ("empty", "cloud"): "#7f7f7f",
        ("empty", "rain"): "#8c564b",
    }
    
    COLORS_GROUND = {
        "tree": "#2ca02c",  
        "fire": "#d62728",  
        "dead": "#333333",  
        "empty": "#333333"
    }

    COLORS_WEATHER = {
        "nothing": "#111111", 
        "cloud": "#cccccc",   
        "rain": "#1f77b4"     
    }

    DEFAULT_PARAMS = {
        'grid_size': 60,
        'grow_p': 0.002,       
        'fire_p': 0.00005,     
        'dead_p': 0.1,         
        'cloud_p': 0.0005,     
        'scat_p': 0.05,        
        'rain_p': 0.01,        
        'rain_end_p': 0.1,     
        'wind_x': 1.0,         
        'wind_y': 0.0
    }

    def __init__(self, root):
        self.root = root
        self.root.title("Симуляция лесного пожара")
        
        self.is_running = False
        self.grid = []
        
        self.canvas_size = 300 
        
        self.entries = {}
        self.active_params = self.DEFAULT_PARAMS.copy()
        
        self.setup_ui()
        self.init_grid()
        self.root.after(50, self.loop)

    def setup_ui(self):
        main_frame = tk.Frame(self.root)
        main_frame.pack(padx=10, pady=10)

        f1 = tk.Frame(main_frame)
        f1.grid(row=0, column=0, padx=5, pady=5)
        tk.Label(f1, text="Общая карта", font=("Arial", 10, "bold")).pack()
        self.canvas_combined = tk.Canvas(f1, width=self.canvas_size, height=self.canvas_size, bg="#555555")
        self.canvas_combined.pack()

        f2 = tk.Frame(main_frame)
        f2.grid(row=0, column=1, padx=5, pady=5)
        tk.Label(f2, text="Земля (Лес и Пожары)", font=("Arial", 10, "bold")).pack()
        self.canvas_ground = tk.Canvas(f2, width=self.canvas_size, height=self.canvas_size, bg="#333333")
        self.canvas_ground.pack()

        f3 = tk.Frame(main_frame)
        f3.grid(row=1, column=0, padx=5, pady=5)
        tk.Label(f3, text="Погода (Тучи и Дождь)", font=("Arial", 10, "bold")).pack()
        self.canvas_weather = tk.Canvas(f3, width=self.canvas_size, height=self.canvas_size, bg="#111111")
        self.canvas_weather.pack()

        for cv in (self.canvas_combined, self.canvas_ground, self.canvas_weather):
            cv.bind("<Button-1>", self.paint_fire)
            cv.bind("<B1-Motion>", self.paint_fire)
            
            cv.bind("<Button-3>", self.paint_cloud)
            cv.bind("<B3-Motion>", self.paint_cloud)
            
            cv.bind("<Button-2>", self.paint_cloud)
            cv.bind("<B2-Motion>", self.paint_cloud)
            
            cv.bind("<Control-Button-1>", self.paint_cloud)
            cv.bind("<Control-B1-Motion>", self.paint_cloud)

        self.control_frame = tk.Frame(main_frame, width=self.canvas_size)
        self.control_frame.grid(row=1, column=1, padx=5, pady=5, sticky="nsew")

        self.btn_play = tk.Button(self.control_frame, text="▶ Старт", bg="lightgreen", font=("Arial", 11, "bold"), command=self.toggle_play)
        self.btn_play.pack(fill=tk.X, pady=(0, 5))
        
        params_container = tk.Frame(self.control_frame)
        params_container.pack(fill=tk.X)
        
        col1 = tk.Frame(params_container)
        col1.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=(0, 2))
        col2 = tk.Frame(params_container)
        col2.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True, padx=(2, 0))

        frame_grid = tk.LabelFrame(col1, text="Поле")
        frame_grid.pack(fill=tk.X, pady=2)
        self.add_entry(frame_grid, "Размер:", "grid_size", self.DEFAULT_PARAMS['grid_size'])
        tk.Button(frame_grid, text="Сброс", command=self.init_grid).grid(row=1, column=0, columnspan=2, pady=2, sticky="ew")

        frame_prob = tk.LabelFrame(col1, text="Лес и огонь")
        frame_prob.pack(fill=tk.X, pady=2)
        self.add_entry(frame_prob, "Рост:", "grow_p", self.DEFAULT_PARAMS['grow_p'])
        self.add_entry(frame_prob, "Пожар:", "fire_p", self.DEFAULT_PARAMS['fire_p'])
        self.add_entry(frame_prob, "Выгорание:", "dead_p", self.DEFAULT_PARAMS['dead_p'])

        frame_weather = tk.LabelFrame(col2, text="Погода")
        frame_weather.pack(fill=tk.X, pady=2)
        self.add_entry(frame_weather, "Туча:", "cloud_p", self.DEFAULT_PARAMS['cloud_p'])
        self.add_entry(frame_weather, "Рассеив.:", "scat_p", self.DEFAULT_PARAMS['scat_p'])
        self.add_entry(frame_weather, "Дождь:", "rain_p", self.DEFAULT_PARAMS['rain_p'])
        self.add_entry(frame_weather, "Конец:", "rain_end_p", self.DEFAULT_PARAMS['rain_end_p'])

        frame_wind = tk.LabelFrame(col2, text="Ветер")
        frame_wind.pack(fill=tk.X, pady=2)
        self.add_entry(frame_wind, "Ось X:", "wind_x", self.DEFAULT_PARAMS['wind_x'])
        self.add_entry(frame_wind, "Ось Y:", "wind_y", self.DEFAULT_PARAMS['wind_y'])

        btn_apply = tk.Button(self.control_frame, text="Применить параметры", bg="lightblue", font=("Arial", 10, "bold"), command=self.apply_parameters)
        btn_apply.pack(fill=tk.X, pady=5)
        
        tk.Label(self.control_frame, text="ЛКМ - Поджечь\nПКМ (или клик 2 пальцами) - Туча", fg="gray").pack()

    def add_entry(self, parent, text, var_name, default_val):
        row = parent.grid_size()[1]
        tk.Label(parent, text=text).grid(row=row, column=0, sticky="w", padx=2, pady=1)
        entry = tk.Entry(parent, width=7, justify="right")
        entry.insert(0, str(default_val))
        entry.grid(row=row, column=1, sticky="e", padx=2, pady=1)
        self.entries[var_name] = entry

    def apply_parameters(self):
        try:
            new_params = {}
            for key, entry in self.entries.items():
                val = entry.get().strip()
                if key == 'grid_size':
                    new_params[key] = int(val)
                else:
                    new_params[key] = float(val)
            
            self.active_params = new_params
            Cell.wind_x = self.active_params['wind_x']
            Cell.wind_y = self.active_params['wind_y']

        except ValueError:
            messagebox.showerror("Ошибка", "Вводите только числа. Разделитель - точка.")

    def toggle_play(self):
        self.is_running = not self.is_running
        if self.is_running:
            self.btn_play.config(text="⏸ Пауза", bg="salmon")
        else:
            self.btn_play.config(text="▶ Старт", bg="lightgreen")

    def init_grid(self):
        try:
            size = int(self.entries['grid_size'].get())
            self.active_params['grid_size'] = size
        except ValueError:
            size = self.active_params['grid_size']

        for cv in (self.canvas_combined, self.canvas_ground, self.canvas_weather):
            cv.delete("all")
            
        self.grid = []
        self.cell_w = self.canvas_size / size

        for y in range(size):
            row = []
            for x in range(size):
                cell = Cell(
                    self.active_params['grow_p'], self.active_params['fire_p'],
                    self.active_params['cloud_p'], self.active_params['dead_p'],
                    self.active_params['rain_p'], self.active_params['rain_end_p'],
                    self.active_params['scat_p']
                )
                x1, y1 = x * self.cell_w, y * self.cell_w
                x2, y2 = x1 + self.cell_w, y1 + self.cell_w
                
                c0 = self.COLORS_COMBINED.get((cell.status, cell.substatus), "black")
                c1 = self.COLORS_GROUND.get(cell.status, "black")
                c2 = self.COLORS_WEATHER.get(cell.substatus, "black")
                
                cell.rect_ids[0] = self.canvas_combined.create_rectangle(x1, y1, x2, y2, fill=c0, outline="")
                cell.rect_ids[1] = self.canvas_ground.create_rectangle(x1, y1, x2, y2, fill=c1, outline="")
                cell.rect_ids[2] = self.canvas_weather.create_rectangle(x1, y1, x2, y2, fill=c2, outline="")
                
                cell.current_colors = [c0, c1, c2]
                row.append(cell)

            self.grid.append(row)

        for y in range(size):
            for x in range(size):
                neighbors = {}
                for dy, key_y in zip([-1, 0, 1], ['t', 'm', 'b']):
                    for dx, key_x in zip([-1, 0, 1], ['l', 'c', 'r']):
                        if dx == 0 and dy == 0:
                            continue

                        nx, ny = x + dx, y + dy

                        if 0 <= nx < size and 0 <= ny < size:
                            neighbors[f"{key_y}{key_x}"] = self.grid[ny][nx]
                        else:
                            neighbors[f"{key_y}{key_x}"] = None

                self.grid[y][x].set_neighbors(**neighbors)

    def paint_fire(self, event):
        self._paint_cell(event, status="fire")

    def paint_cloud(self, event):
        self._paint_cell(event, substatus="cloud")

    def _paint_cell(self, event, status=None, substatus=None):
        x, y = int(event.x // self.cell_w), int(event.y // self.cell_w)
        size = self.active_params['grid_size']
        
        if 0 <= x < size and 0 <= y < size:
            cell = self.grid[y][x]

            if status:
                cell.status = status

            if substatus:
                cell.substatus = substatus
            
            c0 = self.COLORS_COMBINED.get((cell.status, cell.substatus), "black")
            c1 = self.COLORS_GROUND.get(cell.status, "black")
            c2 = self.COLORS_WEATHER.get(cell.substatus, "black")

            if cell.current_colors[0] != c0:
                self.canvas_combined.itemconfig(cell.rect_ids[0], fill=c0)
                cell.current_colors[0] = c0

            if cell.current_colors[1] != c1:
                self.canvas_ground.itemconfig(cell.rect_ids[1], fill=c1)
                cell.current_colors[1] = c1

            if cell.current_colors[2] != c2:
                self.canvas_weather.itemconfig(cell.rect_ids[2], fill=c2)
                cell.current_colors[2] = c2

    def loop(self):
        if self.is_running:
            params = (
                self.active_params['grow_p'], self.active_params['fire_p'],
                self.active_params['cloud_p'], self.active_params['dead_p'],
                self.active_params['rain_p'], self.active_params['rain_end_p'],
                self.active_params['scat_p']
            )

            for row in self.grid:
                for cell in row:
                    cell.update_params(*params)
                    cell.tick()

            for row in self.grid:
                for cell in row:
                    cell.apply_tick()
                    
            self.draw_cells()

        self.root.after(50, self.loop)

    def draw_cells(self):
        for row in self.grid:
            for cell in row:
                c0 = self.COLORS_COMBINED.get((cell.status, cell.substatus), "black")
                if cell.current_colors[0] != c0:
                    self.canvas_combined.itemconfig(cell.rect_ids[0], fill=c0)
                    cell.current_colors[0] = c0

                c1 = self.COLORS_GROUND.get(cell.status, "black")
                if cell.current_colors[1] != c1:
                    self.canvas_ground.itemconfig(cell.rect_ids[1], fill=c1)
                    cell.current_colors[1] = c1

                c2 = self.COLORS_WEATHER.get(cell.substatus, "black")
                if cell.current_colors[2] != c2:
                    self.canvas_weather.itemconfig(cell.rect_ids[2], fill=c2)
                    cell.current_colors[2] = c2


if __name__ == "__main__":
    root = tk.Tk()
    app = SimulationApp(root)
    root.mainloop()

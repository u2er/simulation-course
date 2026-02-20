import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.colors import ListedColormap
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import tkinter as tk
from tkinter import ttk
from typing import Tuple


class ForestFireAutomaton:
    """
    Cellular automaton for simulating forest fire spread.
    
    States:
    0 - Empty/Burnt ground
    1 - Healthy tree
    2 - Burning tree
    3 - Burnt tree (cooling)
    """
    
    EMPTY = 0
    HEALTHY = 1
    BURNING = 2
    BURNT = 3
    
    def __init__(self, width: int, height: int, tree_density: float = 0.6,
                 ignition_prob: float = 0.001, wind_direction: Tuple[int, int] = (1, 0)):
        """
        Initialize the forest fire automaton.
        
        Args:
            width: Width of the grid
            height: Height of the grid
            tree_density: Initial proportion of trees (0-1)
            ignition_prob: Probability of spontaneous ignition per step
            wind_direction: Direction of wind (dx, dy) affecting fire spread
        """
        self.width = width
        self.height = height
        self.ignition_prob = ignition_prob
        self.wind_direction = np.array(wind_direction) / (np.linalg.norm(wind_direction) + 1e-6)
        
        # Initialize grid
        self.grid = np.zeros((height, width), dtype=int)
        self.grid[np.random.random((height, width)) < tree_density] = self.HEALTHY
        
        # Burning duration for each cell
        self.burn_duration = np.zeros((height, width), dtype=int)
        
        # Initialize a fire
        center_x, center_y = width // 2, height // 2
        self.grid[center_y, center_x] = self.BURNING
        self.burn_duration[center_y, center_x] = 3
        
        self.generation = 0
        
    def _get_neighbors(self, x: int, y: int, radius: int = 1) -> list:
        """Get valid neighbors within radius."""
        neighbors = []
        for dx in range(-radius, radius + 1):
            for dy in range(-radius, radius + 1):
                if dx == 0 and dy == 0:
                    continue
                nx, ny = (x + dx) % self.width, (y + dy) % self.height
                neighbors.append((nx, ny))
        return neighbors
    
    def _wind_spread_probability(self, from_x: int, from_y: int, to_x: int, to_y: int) -> float:
        """
        Rule 1: Wind-influenced fire spread.
        Fire spreads more easily in the direction of wind.
        """
        dx = to_x - from_x
        dy = to_y - from_y
        distance = np.sqrt(dx**2 + dy**2)
        
        if distance == 0:
            return 0.0
        
        direction = np.array([dx, dy]) / distance
        wind_alignment = np.dot(direction, self.wind_direction)
        
        # Base spread probability increases with wind alignment
        base_prob = 0.6
        return base_prob + 0.3 * max(0, wind_alignment)
    
    def _age_based_flammability(self, burn_time: int) -> float:
        """
        Rule 2: Age-based flammability.
        Fresh trees burn easier, older burnt trees stop spreading fire.
        """
        if burn_time <= 0:
            return 1.0
        if burn_time <= 2:
            return 0.8
        elif burn_time <= 4:
            return 0.5
        else:
            return 0.2
    
    def _distance_decay(self, distance: int) -> float:
        """
        Rule 3: Distance-based decay.
        Fire spreads less effectively over longer distances.
        """
        if distance == 1:
            return 1.0
        elif distance <= 2:
            return 0.6
        else:
            return 0.3
    
    def _spontaneous_ignition(self) -> float:
        """
        Rule 4: Temperature-based spontaneous ignition.
        Hotter regions (more burning neighbors) increase ignition probability.
        """
        return self.ignition_prob
    
    def _humidity_effect(self, burning_neighbors: int) -> float:
        """
        Rule 5: Humidity-based burning duration.
        More burning neighbors increase heat, reducing burn duration.
        """
        base_duration = 3
        if burning_neighbors > 2:
            return max(1, base_duration - 1)
        return base_duration
    
    def step(self):
        """Execute one step of the cellular automaton."""
        new_grid = self.grid.copy()
        new_burn_duration = self.burn_duration.copy()
        
        for y in range(self.height):
            for x in range(self.width):
                current_state = self.grid[y, x]
                
                if current_state == self.HEALTHY:
                    # Check if any neighbors are burning
                    neighbors = self._get_neighbors(x, y, radius=2)
                    burning_neighbors = []
                    
                    for nx, ny in neighbors:
                        if self.grid[ny, nx] == self.BURNING:
                            burning_neighbors.append((nx, ny))
                    
                    if burning_neighbors:
                        # Fire spread calculation
                        spread_prob = 0.0
                        for bx, by in burning_neighbors:
                            dx = bx - x
                            dy = by - y
                            distance = int(np.ceil(np.sqrt(dx**2 + dy**2)))
                            
                            wind_prob = self._wind_spread_probability(bx, by, x, y)
                            distance_prob = self._distance_decay(distance)
                            flammability = self._age_based_flammability(0)
                            
                            spread_prob += wind_prob * distance_prob * flammability
                        
                        if np.random.random() < min(spread_prob, 0.95):
                            new_grid[y, x] = self.BURNING
                            new_burn_duration[y, x] = int(self._humidity_effect(len(burning_neighbors)))
                    
                    # Spontaneous ignition
                    elif np.random.random() < self._spontaneous_ignition():
                        new_grid[y, x] = self.BURNING
                        new_burn_duration[y, x] = 3
                
                elif current_state == self.BURNING:
                    # Update burn duration
                    new_burn_duration[y, x] -= 1
                    
                    if new_burn_duration[y, x] <= 0:
                        new_grid[y, x] = self.BURNT
                        new_burn_duration[y, x] = 0
                
                elif current_state == self.BURNT:
                    # Slow regrowth - trees eventually regrow
                    if np.random.random() < 0.002:  # Very slow regrowth
                        new_grid[y, x] = self.HEALTHY
        
        self.grid = new_grid
        self.burn_duration = new_burn_duration
        self.generation += 1
    
    def get_state(self) -> np.ndarray:
        """Get current grid state."""
        return self.grid.copy()
    
    def reset(self, tree_density: float):
        """Reset the simulation with new parameters."""
        self.grid = np.zeros((self.height, self.width), dtype=int)
        self.grid[np.random.random((self.height, self.width)) < tree_density] = self.HEALTHY
        self.burn_duration = np.zeros((self.height, self.width), dtype=int)
        
        center_x, center_y = self.width // 2, self.height // 2
        self.grid[center_y, center_x] = self.BURNING
        self.burn_duration[center_y, center_x] = 3
        self.generation = 0


class SimulationGUI:
    """GUI for controlling the forest fire simulation."""
    
    def __init__(self, root):
        self.root = root
        self.root.title("Forest Fire Cellular Automaton Simulator")
        self.root.geometry("1200x700")
        
        # Simulation parameters
        self.width = 100
        self.height = 100
        self.tree_density = 0.7
        self.ignition_prob = 0.0005
        self.wind_x = 1.0
        self.wind_y = 0.3
        
        self.automaton = None
        self.is_running = False
        self.animation = None
        
        self._create_widgets()
        self._create_simulation()
    
    def _create_widgets(self):
        """Create GUI widgets."""
        # Main frame
        main_frame = ttk.Frame(self.root)
        main_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        
        # Left panel - Controls
        left_frame = ttk.Frame(main_frame)
        left_frame.pack(side=tk.LEFT, fill=tk.BOTH, padx=5)
        
        # Title
        title = ttk.Label(left_frame, text="Simulation Controls", font=("Arial", 14, "bold"))
        title.pack(pady=10)
        
        # Grid Size Frame
        grid_frame = ttk.LabelFrame(left_frame, text="Grid Size")
        grid_frame.pack(fill=tk.X, padx=5, pady=5)
        
        ttk.Label(grid_frame, text="Width:").grid(row=0, column=0, sticky=tk.W, padx=5, pady=5)
        self.width_entry = ttk.Entry(grid_frame, width=10)
        self.width_entry.insert(0, str(self.width))
        self.width_entry.grid(row=0, column=1, padx=5, pady=5)
        
        ttk.Label(grid_frame, text="Height:").grid(row=1, column=0, sticky=tk.W, padx=5, pady=5)
        self.height_entry = ttk.Entry(grid_frame, width=10)
        self.height_entry.insert(0, str(self.height))
        self.height_entry.grid(row=1, column=1, padx=5, pady=5)
        
        # Forest Parameters Frame
        forest_frame = ttk.LabelFrame(left_frame, text="Forest Parameters")
        forest_frame.pack(fill=tk.X, padx=5, pady=5)
        
        ttk.Label(forest_frame, text="Tree Density:").grid(row=0, column=0, sticky=tk.W, padx=5, pady=5)
        self.tree_density_entry = ttk.Entry(forest_frame, width=10)
        self.tree_density_entry.insert(0, str(self.tree_density))
        self.tree_density_entry.grid(row=0, column=1, padx=5, pady=5)
        
        ttk.Label(forest_frame, text="(0.0 - 1.0)").grid(row=0, column=2, sticky=tk.W, padx=5)
        
        ttk.Label(forest_frame, text="Ignition Prob:").grid(row=1, column=0, sticky=tk.W, padx=5, pady=5)
        self.ignition_entry = ttk.Entry(forest_frame, width=10)
        self.ignition_entry.insert(0, str(self.ignition_prob))
        self.ignition_entry.grid(row=1, column=1, padx=5, pady=5)
        
        ttk.Label(forest_frame, text="(0.0 - 1.0)").grid(row=1, column=2, sticky=tk.W, padx=5)
        
        # Wind Parameters Frame
        wind_frame = ttk.LabelFrame(left_frame, text="Wind Direction")
        wind_frame.pack(fill=tk.X, padx=5, pady=5)
        
        ttk.Label(wind_frame, text="Wind X:").grid(row=0, column=0, sticky=tk.W, padx=5, pady=5)
        self.wind_x_entry = ttk.Entry(wind_frame, width=10)
        self.wind_x_entry.insert(0, str(self.wind_x))
        self.wind_x_entry.grid(row=0, column=1, padx=5, pady=5)
        
        ttk.Label(wind_frame, text="Wind Y:").grid(row=1, column=0, sticky=tk.W, padx=5, pady=5)
        self.wind_y_entry = ttk.Entry(wind_frame, width=10)
        self.wind_y_entry.insert(0, str(self.wind_y))
        self.wind_y_entry.grid(row=1, column=1, padx=5, pady=5)
        
        # Control Buttons Frame
        button_frame = ttk.LabelFrame(left_frame, text="Controls")
        button_frame.pack(fill=tk.X, padx=5, pady=5)
        
        self.start_button = ttk.Button(button_frame, text="Start", command=self._toggle_simulation)
        self.start_button.pack(fill=tk.X, padx=5, pady=5)
        
        self.reset_button = ttk.Button(button_frame, text="Reset", command=self._reset_simulation)
        self.reset_button.pack(fill=tk.X, padx=5, pady=5)
        
        self.step_button = ttk.Button(button_frame, text="Step", command=self._step_simulation)
        self.step_button.pack(fill=tk.X, padx=5, pady=5)
        
        # Speed Control Frame
        speed_frame = ttk.LabelFrame(left_frame, text="Animation Speed")
        speed_frame.pack(fill=tk.X, padx=5, pady=5)
        
        ttk.Label(speed_frame, text="Delay (ms):").grid(row=0, column=0, sticky=tk.W, padx=5, pady=5)
        self.speed_entry = ttk.Entry(speed_frame, width=10)
        self.speed_entry.insert(0, "100")
        self.speed_entry.grid(row=0, column=1, padx=5, pady=5)
        
        # Statistics Frame
        stats_frame = ttk.LabelFrame(left_frame, text="Statistics")
        stats_frame.pack(fill=tk.X, padx=5, pady=5)
        
        self.generation_label = ttk.Label(stats_frame, text="Generation: 0")
        self.generation_label.pack(anchor=tk.W, padx=5, pady=5)
        
        self.trees_label = ttk.Label(stats_frame, text="Healthy Trees: 0")
        self.trees_label.pack(anchor=tk.W, padx=5, pady=5)
        
        self.burning_label = ttk.Label(stats_frame, text="Burning Trees: 0")
        self.burning_label.pack(anchor=tk.W, padx=5, pady=5)
        
        self.burnt_label = ttk.Label(stats_frame, text="Burnt Ground: 0")
        self.burnt_label.pack(anchor=tk.W, padx=5, pady=5)
        
        # Right panel - Visualization
        right_frame = ttk.Frame(main_frame)
        right_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True, padx=5)
        
        self.canvas_frame = ttk.Frame(right_frame)
        self.canvas_frame.pack(fill=tk.BOTH, expand=True)
    
    def _create_simulation(self):
        """Create the initial simulation."""
        try:
            self.width = int(self.width_entry.get())
            self.height = int(self.height_entry.get())
            self.tree_density = float(self.tree_density_entry.get())
            self.ignition_prob = float(self.ignition_entry.get())
            self.wind_x = float(self.wind_x_entry.get())
            self.wind_y = float(self.wind_y_entry.get())
        except ValueError:
            return
        
        self.automaton = ForestFireAutomaton(
            width=self.width,
            height=self.height,
            tree_density=self.tree_density,
            ignition_prob=self.ignition_prob,
            wind_direction=(self.wind_x, self.wind_y)
        )
        
        self._update_canvas()
    
    def _update_canvas(self):
        """Update the matplotlib canvas."""
        # Clear previous canvas
        for widget in self.canvas_frame.winfo_children():
            widget.destroy()
        
        # Create new figure
        fig = Figure(figsize=(6, 6), dpi=100)
        ax = fig.add_subplot(111)
        
        colors = ['white', 'green', 'red', 'gray']
        cmap = ListedColormap(colors)
        
        self.im = ax.imshow(self.automaton.get_state(), cmap=cmap, vmin=0, vmax=3)
        ax.set_title("Forest Fire Simulation")
        
        # Embed in tkinter
        canvas = FigureCanvasTkAgg(fig, master=self.canvas_frame)
        canvas.draw()
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        
        self.figure = fig
        self.ax = ax
        self.canvas = canvas
        
        self._update_statistics()
    
    def _update_statistics(self):
        """Update statistics labels."""
        if self.automaton is None:
            return
        
        grid = self.automaton.grid
        healthy = np.sum(grid == ForestFireAutomaton.HEALTHY)
        burning = np.sum(grid == ForestFireAutomaton.BURNING)
        burnt = np.sum(grid == ForestFireAutomaton.BURNT)
        
        self.generation_label.config(text=f"Generation: {self.automaton.generation}")
        self.trees_label.config(text=f"Healthy Trees: {healthy}")
        self.burning_label.config(text=f"Burning Trees: {burning}")
        self.burnt_label.config(text=f"Burnt Ground: {burnt}")
    
    def _toggle_simulation(self):
        """Toggle simulation running state."""
        if self.automaton is None:
            self._create_simulation()
        
        self.is_running = not self.is_running
        self.start_button.config(text="Pause" if self.is_running else "Start")
        
        if self.is_running:
            self._animate()
    
    def _animate(self):
        """Animate the simulation."""
        if not self.is_running or self.automaton is None:
            return
        
        try:
            delay = int(self.speed_entry.get())
        except ValueError:
            delay = 100
        
        self.automaton.step()
        self.im.set_data(self.automaton.get_state())
        self.canvas.draw_idle()
        self._update_statistics()
        
        self.root.after(delay, self._animate)
    
    def _step_simulation(self):
        """Execute one step of the simulation."""
        if self.automaton is None:
            self._create_simulation()
        
        self.automaton.step()
        self.im.set_data(self.automaton.get_state())
        self.canvas.draw_idle()
        self._update_statistics()
    
    def _reset_simulation(self):
        """Reset the simulation with current parameters."""
        self.is_running = False
        self.start_button.config(text="Start")
        
        try:
            self.width = int(self.width_entry.get())
            self.height = int(self.height_entry.get())
            self.tree_density = float(self.tree_density_entry.get())
            self.ignition_prob = float(self.ignition_entry.get())
            self.wind_x = float(self.wind_x_entry.get())
            self.wind_y = float(self.wind_y_entry.get())
        except ValueError:
            return
        
        self.automaton = ForestFireAutomaton(
            width=self.width,
            height=self.height,
            tree_density=self.tree_density,
            ignition_prob=self.ignition_prob,
            wind_direction=(self.wind_x, self.wind_y)
        )
        
        self._update_canvas()


def main():
    root = tk.Tk()
    gui = SimulationGUI(root)
    root.mainloop()


if __name__ == '__main__':
    main()
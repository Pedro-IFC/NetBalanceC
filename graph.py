import tkinter as tk
from tkinter import ttk
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import csv
import numpy as np
import random

# Constantes
N = 20
B_VALUE = 100.0

def generate_tester(min_matrix, max_matrix):
    """
    Gera uma matriz de teste com base nos valores mínimos e máximos.
    Esta função simula a lógica do código C original.
    """
    tester = [[0.0] * N for _ in range(N)]
    for i in range(N):
        for j in range(i, N): # Itera na triangular superior
            min_val = int(min_matrix[i][j])
            max_val = int(max_matrix[i][j])
            v = 0.0
            if max_val > 0 and max_val >= min_val:
                steps = (max_val - min_val) // 10 + 1
                v = float(min_val + 10 * random.randint(0, steps - 1))
            tester[i][j] = v
            tester[j][i] = v # Garante simetria
    return tester

def load_history_data(filepath="history_advanced.csv"):
    generations, global_fitness, generational_fitness, genes_history = [], [], [], []
    try:
        with open(filepath, newline='', encoding="utf-8") as csvfile:
            reader = csv.reader(csvfile)
            next(reader, None)
            for row in reader:
                generations.append(int(row[0]))
                global_fitness.append(float(row[1]))
                generational_fitness.append(float(row[2]))
                genes_history.append([int(g) for g in row[3:]])
    except FileNotFoundError:
        print(f"Erro: O arquivo '{filepath}' não foi encontrado.")
    return generations, global_fitness, generational_fitness, genes_history

def load_tester_matrices(filepath="tester.csv"):
    min_matrix = [[0.0] * N for _ in range(N)]
    max_matrix = [[0.0] * N for _ in range(N)]
    try:
        with open(filepath, newline='', encoding="utf-8") as csvfile:
            reader = csv.reader(csvfile)
            next(reader, None)
            for row in reader:
                i, j, min_val, max_val = map(float, row)
                min_matrix[int(i)][int(j)] = min_val
                min_matrix[int(j)][int(i)] = min_val
                max_matrix[int(i)][int(j)] = max_val
                max_matrix[int(j)][int(i)] = max_val
    except FileNotFoundError:
        print(f"Erro: O arquivo '{filepath}' não foi encontrado.")
        return None, None
    return min_matrix, max_matrix

class GeneticAlgoAnalyzer(tk.Tk):
    def __init__(self, generations, global_fitness, generational_fitness, genes_history, min_matrix, max_matrix):
        super().__init__()
        self.title("Análise de Algoritmo Genético")
        self.state('zoomed')

        self.generations = generations
        self.global_fitness = global_fitness
        self.generational_fitness = generational_fitness
        self.genes_history = genes_history
        self.min_matrix = min_matrix
        self.max_matrix = max_matrix
        self.b_vector = [B_VALUE] * N
        self.initial_solution_x = self._calculate_solution(0)
        self.active_nodes = [tk.BooleanVar(value=True) for _ in range(N)]

        main_frame = ttk.Frame(self)
        main_frame.pack(fill=tk.BOTH, expand=1)
        self.canvas = tk.Canvas(main_frame)
        self.canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=1)
        scrollbar = ttk.Scrollbar(main_frame, orient=tk.VERTICAL, command=self.canvas.yview)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        self.canvas.configure(yscrollcommand=scrollbar.set)
        
        self.scrollable_frame = ttk.Frame(self.canvas)
        self.canvas_window_id = self.canvas.create_window((0, 0), window=self.scrollable_frame, anchor="nw")

        # MODIFICAÇÃO: Vincula a função _on_canvas_configure ao evento de redimensionamento do canvas
        # Isso garante que a área rolável sempre ocupe 100% da largura
        self.canvas.bind('<Configure>', self._on_canvas_configure)
        self.canvas.bind_all("<MouseWheel>", self._on_mousewheel)
        
        self._configure_grid()
        self._create_widgets()
        self._initial_display_update()

    # MODIFICAÇÃO: Nova função para lidar com o redimensionamento do canvas
    def _on_canvas_configure(self, event):
        canvas_width = event.width
        # Ajusta a largura do frame interno para ser igual à largura do canvas
        self.canvas.itemconfig(self.canvas_window_id, width=canvas_width)
        # Atualiza a região de rolagem
        self.canvas.configure(scrollregion=self.canvas.bbox("all"))

    def _on_mousewheel(self, event):
        self.canvas.yview_scroll(int(-1*(event.delta/120)), "units")

    def _configure_grid(self):
        self.scrollable_frame.grid_columnconfigure(0, weight=1)

    def _create_widgets(self):
        self._create_fitness_panel()
        self._create_genes_panel()
        self._create_solution_panel()
        self._create_simulation_panel()

    def _create_fitness_panel(self):
        frame = ttk.LabelFrame(self.scrollable_frame, text="Evolução do Fitness", padding="10")
        frame.grid(row=0, column=0, sticky="nsew", padx=10, pady=5)
        fig, ax = plt.subplots(figsize=(10, 5)) 
        ax.plot(self.generations, self.global_fitness, label="Melhor Fitness Global", color='blue')
        ax.plot(self.generations, self.generational_fitness, label="Fitness Geracional", color='red', alpha=0.7)
        ax.set_xlabel("Geração")
        ax.set_ylabel("Fitness")
        ax.legend()
        ax.grid(True)
        plt.tight_layout()
        canvas = FigureCanvasTkAgg(fig, master=frame)
        canvas.draw()
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

    def _create_genes_panel(self):
        frame = ttk.LabelFrame(self.scrollable_frame, text="Visualização dos Genes e Controles", padding="10")
        frame.grid(row=1, column=0, sticky="nsew", padx=10, pady=5)

        controls_frame = ttk.Frame(frame)
        controls_frame.pack(pady=5, fill=tk.X)
        
        self.slider = tk.Scale(controls_frame, from_=0, to=len(self.generations)-1, orient=tk.HORIZONTAL,
                               command=self.on_slider_change, label="Selecione a Geração")
        self.slider.pack(pady=5, fill=tk.X, expand=True)
        
        entry_frame = ttk.Frame(controls_frame)
        entry_frame.pack(pady=5)
        ttk.Label(entry_frame, text="Ou digite o Nº:").pack(side=tk.LEFT)
        self.entry_gen = ttk.Entry(entry_frame, width=10)
        self.entry_gen.pack(side=tk.LEFT, padx=5)
        ttk.Button(entry_frame, text="Atualizar", command=self.on_manual_entry).pack(side=tk.LEFT)
        
        # MODIFICAÇÃO: justify alterado para tk.CENTER para centralizar o texto dos genes
        self.genes_text = tk.Label(frame, text="", justify=tk.CENTER, font=("Courier", 11))
        self.genes_text.pack(pady=20, fill=tk.BOTH, expand=True)

    def _create_solution_panel(self):
        frame = ttk.LabelFrame(self.scrollable_frame, text="Comparação da Solução do Vetor X", padding="10")
        frame.grid(row=2, column=0, sticky="nsew", padx=10, pady=5)
        frame.grid_rowconfigure(0, weight=1)
        frame.grid_columnconfigure(0, weight=1)

        columns = ('node', 'gen0', 'gen_curr', 'abs_diff', 'perc_diff')
        self.solution_tree = ttk.Treeview(frame, columns=columns, show='headings', height=N)
        
        self.solution_tree.heading('node', text='Nó (i)')
        self.solution_tree.heading('gen0', text='Solução Geração 0')
        self.solution_tree.heading('gen_curr', text='Solução Geração Atual')
        self.solution_tree.heading('abs_diff', text='Dif. Absoluta')
        self.solution_tree.heading('perc_diff', text='Dif. %')
        
        # MODIFICAÇÃO: anchor de todas as colunas alterado para 'center' para centralizar os dados
        self.solution_tree.column('node', width=80, anchor='center')
        self.solution_tree.column('gen0', width=150, anchor='center')
        self.solution_tree.column('gen_curr', width=160, anchor='center')
        self.solution_tree.column('abs_diff', width=140, anchor='center')
        self.solution_tree.column('perc_diff', width=120, anchor='center')

        tree_scrollbar = ttk.Scrollbar(frame, orient="vertical", command=self.solution_tree.yview)
        self.solution_tree.configure(yscrollcommand=tree_scrollbar.set)
        
        self.solution_tree.grid(row=0, column=0, sticky='nsew')
        tree_scrollbar.grid(row=0, column=1, sticky='ns')

    def _create_simulation_panel(self):
        frame = ttk.LabelFrame(self.scrollable_frame, text="Simulação de Desligamento de Nós", padding="10")
        frame.grid(row=3, column=0, sticky="nsew", padx=10, pady=5)

        # MODIFICAÇÃO: Adicionado um frame wrapper para centralizar o frame dos checkboxes
        center_wrapper_frame = ttk.Frame(frame)
        center_wrapper_frame.pack(pady=5, fill="x")

        # Frame interno que conterá os checkboxes e será centralizado
        checkbox_frame = ttk.Frame(center_wrapper_frame)
        checkbox_frame.pack() # pack() sem argumentos centraliza o widget
        
        cols = 10 
        for i in range(N):
            var = self.active_nodes[i]
            # Os checkboxes são adicionados ao frame interno
            cb = ttk.Checkbutton(checkbox_frame, text=str(i), variable=var, command=self.update_simulation_fitness)
            cb.grid(row=i // cols, column=i % cols, padx=5, pady=2, sticky="w")
        
        self.simulation_fitness_label = ttk.Label(frame, text="Fitness Simulado: -", font=("Arial", 14, "bold"))
        self.simulation_fitness_label.pack(pady=15)

    def _initial_display_update(self):
        if self.generations:
            initial_idx = len(self.generations) - 1
            self.slider.set(initial_idx)
            self.entry_gen.insert(0, str(self.generations[initial_idx]))
            self.update_displays(initial_idx)
    
    def calculate_fitness_py(self, A_matrix, b_vector):
        try:
            sol = np.linalg.solve(A_matrix, b_vector)
            total_abs = np.sum(np.abs(sol))
            if total_abs == 0.0:
                return 0.0
            return 10.0 / total_abs
        except np.linalg.LinAlgError:
            return 0.0

    def update_simulation_fitness(self):
        gen_index = self.slider.get()
        if not self.min_matrix or gen_index >= len(self.genes_history):
            self.simulation_fitness_label.config(text="Fitness Simulado: Erro")
            return

        genes_flat = self.genes_history[gen_index]
        positions = np.array([genes_flat[i:i+N] for i in range(0, len(genes_flat), N)])
        
        random.seed(0) 
        tester = np.array(generate_tester(self.min_matrix, self.max_matrix))

        A_sim = positions * tester
        b_sim = np.array(self.b_vector, dtype=float)
        for i in range(N):
            if not self.active_nodes[i].get():
                A_sim[i, :] = 0.0
                A_sim[:, i] = 0.0
                b_sim[i] = 0.0
                A_sim[i][i] = 1
        
        fitness_score = self.calculate_fitness_py(A_sim, b_sim)
        self.simulation_fitness_label.config(text=f"Fitness Simulado: {fitness_score:.4f}")

    def _calculate_solution(self, generation_index):
        if not self.min_matrix or generation_index >= len(self.genes_history):
            return None
        genes_flat = self.genes_history[generation_index]
        positions = np.array([genes_flat[i:i+N] for i in range(0, len(genes_flat), N)])
        random.seed(0) 
        tester = np.array(generate_tester(self.min_matrix, self.max_matrix))
        A = positions * tester
        try:
            return np.linalg.solve(A, np.array(self.b_vector))
        except np.linalg.LinAlgError:
            return None

    def format_genes(self, genes_list):
        return "\n".join(" ".join(map(str, genes_list[i:i+N])) for i in range(0, len(genes_list), N))

    def update_displays(self, generation_index):
        self.update_genes_display(generation_index)
        self.update_solution_display(generation_index)
        self.update_simulation_panel()
        self.scrollable_frame.update_idletasks()
        self.canvas.configure(scrollregion=self.canvas.bbox("all"))

    def update_genes_display(self, generation_index):
        genes_str = self.format_genes(self.genes_history[generation_index])
        self.genes_text.config(text=genes_str)

    def update_simulation_panel(self):
        for var in self.active_nodes:
            var.set(True)
        self.update_simulation_fitness()

    def update_solution_display(self, generation_index):
        for item in self.solution_tree.get_children():
            self.solution_tree.delete(item)
        
        current_solution_x = self._calculate_solution(generation_index)
        
        if self.initial_solution_x is None or current_solution_x is None:
            msg = "Solução inicial não pôde ser calculada." if self.initial_solution_x is None else "Solução atual não pôde ser calculada (matriz singular)."
            self.solution_tree.insert("", "end", values=("Erro", msg, "", "", ""))
            return

        for i in range(N):
            initial_val = self.initial_solution_x[i]
            current_val = current_solution_x[i]
            abs_diff = abs(current_val - initial_val)
            
            if initial_val != 0:
                perc_diff = ((current_val - initial_val) / abs(initial_val)) * 100
                perc_str = f"{perc_diff:+.2f}%"
            else:
                perc_str = "N/A"

            self.solution_tree.insert("", "end", values=(
                f"x[{i}]", f"{initial_val:9.4f}", f"{current_val:9.4f}",
                f"{abs_diff:9.4f}", perc_str
            ))

    def on_slider_change(self, val):
        generation_idx = int(val)
        self.update_displays(generation_idx)
        self.entry_gen.delete(0, tk.END)
        self.entry_gen.insert(0, str(self.generations[generation_idx]))

    def on_manual_entry(self):
        try:
            gen_num = int(self.entry_gen.get())
            if gen_num in self.generations:
                gen_idx = self.generations.index(gen_num)
                self.slider.set(gen_idx)
        except ValueError:
            print("Entrada manual inválida.")

def main():
    generations, global_fitness, generational_fitness, genes_history = load_history_data()
    min_matrix, max_matrix = load_tester_matrices()
    if not generations:
        print("Não foi possível carregar os dados do histórico. Encerrando.")
        return
    app = GeneticAlgoAnalyzer(generations, global_fitness, generational_fitness, genes_history, min_matrix, max_matrix)
    app.mainloop()

if __name__ == '__main__':
    main()
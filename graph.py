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
            next(reader, None)  # Pular cabeçalho
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
                min_matrix[int(j)][int(i)] = min_val # Garante simetria
                max_matrix[int(i)][int(j)] = max_val
                max_matrix[int(j)][int(i)] = max_val # Garante simetria
    except FileNotFoundError:
        print(f"Erro: O arquivo '{filepath}' não foi encontrado.")
        return None, None
    return min_matrix, max_matrix

class GeneticAlgoAnalyzer(tk.Tk):
    def __init__(self, generations, global_fitness, generational_fitness, genes_history, min_matrix, max_matrix):
        super().__init__()
        self.title("Análise de Algoritmo Genético")
        self.geometry("1400x900")

        self.generations = generations
        self.global_fitness = global_fitness
        self.generational_fitness = generational_fitness
        self.genes_history = genes_history
        self.min_matrix = min_matrix
        self.max_matrix = max_matrix
        self.b_vector = [B_VALUE] * N

        # Calcula a solução inicial (Geração 0) uma única vez
        self.initial_solution_x = self._calculate_solution(0)

        # --- Estrutura de Rolagem ---
        main_frame = ttk.Frame(self)
        main_frame.pack(fill=tk.BOTH, expand=1)
        self.canvas = tk.Canvas(main_frame)
        self.canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=1)
        scrollbar = ttk.Scrollbar(main_frame, orient=tk.VERTICAL, command=self.canvas.yview)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        self.canvas.configure(yscrollcommand=scrollbar.set)
        self.canvas.bind('<Configure>', lambda e: self.canvas.configure(scrollregion=self.canvas.bbox("all")))
        self.canvas.bind_all("<MouseWheel>", self._on_mousewheel)
        self.scrollable_frame = ttk.Frame(self.canvas)
        self.canvas.create_window((0, 0), window=self.scrollable_frame, anchor="nw")
        
        # --- Conteúdo da Aplicação ---
        self._configure_grid()
        self._create_widgets()
        self._initial_display_update()

    def _on_mousewheel(self, event):
        self.canvas.yview_scroll(int(-1*(event.delta/120)), "units")

    def _configure_grid(self):
        self.scrollable_frame.grid_columnconfigure(0, weight=1, minsize=700)
        self.scrollable_frame.grid_columnconfigure(1, weight=1, minsize=700)

    def _create_widgets(self):
        self._create_fitness_panel()
        self._create_genes_panel()
        self._create_solution_panel()

    def _create_fitness_panel(self):
        frame = ttk.Frame(self.scrollable_frame, padding="10")
        frame.grid(row=0, column=0, sticky="nsew")
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.plot(self.generations, self.global_fitness, label="Melhor Fitness Global", color='blue')
        ax.plot(self.generations, self.generational_fitness, label="Fitness Geracional", color='red', alpha=0.7)
        ax.set_xlabel("Geração")
        ax.set_ylabel("Fitness")
        ax.set_title("Evolução do Fitness")
        ax.legend()
        ax.grid(True)
        plt.tight_layout()
        canvas = FigureCanvasTkAgg(fig, master=frame)
        canvas.draw()
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

    def _create_genes_panel(self):
        frame = ttk.Frame(self.scrollable_frame, padding="10")
        frame.grid(row=0, column=1, rowspan=2, sticky="nsew")
        ttk.Label(frame, text="Visualização dos Genes", font=("Arial", 16)).pack(pady=10)
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
        self.genes_text = tk.Label(frame, text="", justify=tk.LEFT, font=("Courier", 11))
        self.genes_text.pack(pady=20, fill=tk.BOTH, expand=True)

    def _create_solution_panel(self):
        frame = ttk.Frame(self.scrollable_frame, padding="10")
        frame.grid(row=1, column=0, sticky="nsew")
        frame.grid_rowconfigure(1, weight=1)
        frame.grid_columnconfigure(0, weight=1)

        ttk.Label(frame, text="Comparação da Solução do Vetor X", font=("Arial", 16)).grid(row=0, column=0, pady=10, sticky="ew")

        columns = ('node', 'gen0', 'gen_curr', 'abs_diff', 'perc_diff')
        self.solution_tree = ttk.Treeview(frame, columns=columns, show='headings')
        
        self.solution_tree.heading('node', text='Nó (i)')
        self.solution_tree.heading('gen0', text='Solução Geração 0')
        self.solution_tree.heading('gen_curr', text='Solução Geração Atual')
        self.solution_tree.heading('abs_diff', text='Dif. Absoluta')
        self.solution_tree.heading('perc_diff', text='Dif. %')
        
        self.solution_tree.column('node', width=80, anchor='center')
        self.solution_tree.column('gen0', width=150, anchor='e')
        self.solution_tree.column('gen_curr', width=160, anchor='e')
        self.solution_tree.column('abs_diff', width=140, anchor='e')
        self.solution_tree.column('perc_diff', width=120, anchor='e')

        tree_scrollbar = ttk.Scrollbar(frame, orient="vertical", command=self.solution_tree.yview)
        self.solution_tree.configure(yscrollcommand=tree_scrollbar.set)
        
        self.solution_tree.grid(row=1, column=0, sticky='nsew')
        tree_scrollbar.grid(row=1, column=1, sticky='ns')

    def _initial_display_update(self):
        if self.generations:
            initial_idx = len(self.generations) - 1
            self.slider.set(initial_idx)
            self.entry_gen.insert(0, str(self.generations[initial_idx]))
            self.update_displays(initial_idx)
    
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
        self.scrollable_frame.update_idletasks()
        self.canvas.configure(scrollregion=self.canvas.bbox("all"))

    def update_genes_display(self, generation_index):
        genes_str = self.format_genes(self.genes_history[generation_index])
        self.genes_text.config(text=genes_str)

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
                f"x[{i}]",
                f"{initial_val:9.4f}",
                f"{current_val:9.4f}",
                f"{abs_diff:9.4f}",
                perc_str
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
                # O slider já chama on_slider_change, que chama update_displays
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
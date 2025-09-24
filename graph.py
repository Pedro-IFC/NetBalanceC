import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import tkinter as tk
import csv
import math
import numpy as np
import networkx as nx

def media_movel(dados, janela=10):
    """Calcula a média móvel de uma lista de dados."""
    if len(dados) < janela:
        return dados
    return [sum(dados[i:i+janela]) / janela for i in range(len(dados) - janela + 1)]

# --- Carregar os dados ---
geracoes = []
fitness_global = []
fitness_geracional = []
genes_history = [] 

try:
    with open("history_advanced.csv", newline='', encoding="utf-8") as csvfile:
        reader = csv.reader(csvfile)
        header = next(reader, None) 
        for row in reader:
            geracoes.append(int(row[0]))
            fitness_global.append(float(row[1]))
            fitness_geracional.append(float(row[2]))
            genes = [int(g) for g in row[3:]]
            genes_history.append(genes)
except FileNotFoundError:
    print("Erro: O arquivo 'history_advanced.csv' não foi encontrado.")
    exit()

janela_media_movel = 1
fitness_geracional_suave = media_movel(fitness_geracional, janela_media_movel)
geracoes_suave = geracoes[:len(fitness_geracional_suave)]

# --- Calcular a dimensão da matriz ---
if genes_history:
    num_genes = len(genes_history[0])
    matrix_dim = int(math.sqrt(num_genes))
    if matrix_dim * matrix_dim != num_genes:
        print(f"Atenção: O número de genes ({num_genes}) não é um quadrado perfeito. A visualização da matriz pode não ser ideal.")
else:
    matrix_dim = 0

# --- Funções para a interface ---

# Bloco de código para o grafo
# ----------------------------------------------------
def draw_graph(genes_list):
    """
    Desenha um grafo a partir da lista de genes.
    A lista de genes é tratada como uma matriz de adjacência.
    """
    for widget in frame_graph.winfo_children():
        widget.destroy()

    if not genes_list or matrix_dim == 0:
        label = tk.Label(frame_graph, text="Dados de grafo indisponíveis.")
        label.pack(expand=True)
        return
        
    adj_matrix = np.array(genes_list).reshape(matrix_dim, matrix_dim)
    G = nx.from_numpy_array(adj_matrix, create_using=nx.DiGraph)

    node_colors = plt.cm.get_cmap('tab10', matrix_dim)
    node_color_list = [node_colors(i) for i in range(matrix_dim)]
    
    pos = nx.circular_layout(G)
    
    fig_graph, ax_graph = plt.subplots(figsize=(6, 4), dpi=100)
    
    def on_hover(event):
        if event.inaxes != ax_graph:
            return

        ax_graph.clear()
        
        nx.draw_networkx_nodes(G, pos, ax=ax_graph, node_color=node_color_list, node_size=700)
        nx.draw_networkx_labels(G, pos, ax=ax_graph)

        highlight_node = None
        min_dist = float('inf')
        for node, (x, y) in pos.items():
            dist = math.sqrt((x - event.xdata)**2 + (y - event.ydata)**2)
            if dist < min_dist:
                min_dist = dist
                highlight_node = node
        
        if min_dist < 0.1:
            out_edges = list(G.out_edges(highlight_node))
            if out_edges:
                nx.draw_networkx_edges(G, pos, edgelist=out_edges, ax=ax_graph, 
                                       edge_color='red', width=2, arrows=True)
            
            in_edges = list(G.in_edges(highlight_node))
            if in_edges:
                nx.draw_networkx_edges(G, pos, edgelist=in_edges, ax=ax_graph, 
                                       edge_color='green', width=2, arrows=True)

            nx.draw_networkx_nodes(G, pos, ax=ax_graph, nodelist=[highlight_node],
                                   node_color='yellow', node_size=1000)

        fig_graph.canvas.draw_idle()

    fig_graph.canvas.mpl_connect('motion_notify_event', on_hover)

    ax_graph.set_title(f'Grafo da Geração {geracoes[slider.get()]}')
    nx.draw(G, pos, ax=ax_graph, with_labels=True, node_color=node_color_list,
            edge_color='gray', node_size=700)
    plt.tight_layout()

    canvas_graph = FigureCanvasTkAgg(fig_graph, master=frame_graph)
    canvas_widget_graph = canvas_graph.get_tk_widget()
    canvas_widget_graph.pack(fill=tk.BOTH, expand=True)
    canvas_graph.draw()
# ----------------------------------------------------

def update_genes_display(generation_index):
    """
    Atualiza a exibição dos genes em formato de matriz e desenha o grafo.
    """
    if generation_index < len(geracoes):
        genes = genes_history[generation_index]
        matrix_str = ""
        for i in range(0, len(genes), matrix_dim):
            row = genes[i:i+matrix_dim]
            matrix_str += ' '.join(map(str, row)) + '\n'
        
        genes_text.config(text=f"Melhor Indivíduo da Geração {geracoes[generation_index]}:\n\n{matrix_str}")
        draw_graph(genes)
    else:
        genes_text.config(text="Geração não encontrada.")
        for widget in frame_graph.winfo_children():
            widget.destroy()

def on_slider_change(val):
    """Lida com a mudança da barra de rolagem."""
    gen_index = int(float(val))
    update_genes_display(gen_index)
    entry_gen.delete(0, tk.END)
    entry_gen.insert(0, str(geracoes[gen_index]))

def on_manual_entry():
    """Lida com a entrada manual do número da geração."""
    try:
        gen_number = int(entry_gen.get())
        if gen_number in geracoes:
            gen_index = geracoes.index(gen_number)
            slider.set(gen_index)  # Isso irá acionar on_slider_change automaticamente
        else:
            tk.messagebox.showwarning("Geração Inválida", "O número de geração inserido não existe.")
    except ValueError:
        tk.messagebox.showerror("Erro de Entrada", "Por favor, insira um número inteiro válido.")

# --- Configuração da interface Tkinter ---
root = tk.Tk()
root.title("Visualizador de Histórico Genético")
root.attributes('-fullscreen', True)

root.bind('<Escape>', lambda event: root.attributes('-fullscreen', False))

root.grid_columnconfigure(0, weight=1)
root.grid_columnconfigure(1, weight=1)
root.grid_rowconfigure(0, weight=1)
root.grid_rowconfigure(1, weight=1)

# Painel Superior Esquerdo - Gráfico de Fitness
frame_top_left = tk.Frame(root, bg="white")
frame_top_left.grid(row=0, column=0, sticky="nsew", padx=10, pady=10)

fig, ax = plt.subplots(figsize=(6, 4), dpi=100)
ax.set_title('Histórico de Fitness')
ax.set_xlabel('Geração')
ax.set_ylabel('Fitness')
ax.grid(True)
ax.plot(geracoes, fitness_global, color='dodgerblue', linewidth=2.5, marker='o', markersize=4,
        label='Melhor Fitness Global Acumulado')
ax.plot(geracoes_suave, fitness_geracional_suave, color='tomato', linewidth=2,
        label=f'Melhor Fitness Geracional (Média Móvel {janela_media_movel})')
ax.legend()
plt.tight_layout()

canvas = FigureCanvasTkAgg(fig, master=frame_top_left)
canvas_widget = canvas.get_tk_widget()
canvas_widget.pack(fill=tk.BOTH, expand=True)

# Painel Superior Direito - Seleção de Geração e Genes
frame_genes = tk.Frame(root)
frame_genes.grid(row=0, column=1, sticky="nsew", padx=10, pady=10)

title_label = tk.Label(frame_genes, text="Selecione uma Geração", font=("Arial", 16))
title_label.pack(pady=10)

if geracoes:
    slider = tk.Scale(frame_genes, from_=0, to=len(geracoes)-1, orient=tk.HORIZONTAL,
                      command=on_slider_change, length=400,
                      label="Nº da Geração")
    slider.pack(pady=5)
    
    # Campo de entrada manual e botão
    frame_entry = tk.Frame(frame_genes)
    frame_entry.pack(pady=5)
    tk.Label(frame_entry, text="Ou digite o Nº:").pack(side=tk.LEFT)
    entry_gen = tk.Entry(frame_entry, width=10)
    entry_gen.pack(side=tk.LEFT, padx=5)
    btn_update = tk.Button(frame_entry, text="Atualizar", command=on_manual_entry)
    btn_update.pack(side=tk.LEFT)

    slider.set(len(geracoes)-1)
    entry_gen.insert(0, str(geracoes[len(geracoes)-1]))

else:
    no_data_label = tk.Label(frame_genes, text="Nenhum dado de geração encontrado.", font=("Arial", 12))
    no_data_label.pack(pady=10)
    slider = None # Para evitar erros se não houver dados
    entry_gen = None

genes_text = tk.Label(frame_genes, text="", justify=tk.LEFT, font=("Courier", 10))
genes_text.pack(pady=20, fill=tk.BOTH, expand=True)

# Painel Inferior Esquerdo
frame_bottom_left = tk.Frame(root, bg="lightgray")
frame_bottom_left.grid(row=1, column=0, sticky="nsew", padx=10, pady=10)
tk.Label(frame_bottom_left, text="Painel Inferior Esquerdo", bg="lightgray").pack(expand=True)

# Painel Inferior Direito - Grafo (local ajustado)
frame_graph = tk.Frame(root, bg="white")
frame_graph.grid(row=1, column=1, sticky="nsew", padx=10, pady=10)

if geracoes:
    update_genes_display(len(geracoes) - 1)

root.mainloop()
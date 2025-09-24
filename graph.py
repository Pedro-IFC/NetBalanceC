import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import csv
import numpy as np # Adicionado para manipulação de matrizes

def media_movel(dados, janela=10):
    """Calcula a média móvel de uma lista de dados."""
    if len(dados) < janela:
        return dados
    return [sum(dados[i:i+janela]) / janela for i in range(len(dados) - janela + 1)]

# --- Leitura dos Dados do Novo Arquivo ---
geracoes = []
fitness_global = []
fitness_geracional = []
genes_history = [] # Lista para guardar os genes de cada geração

# MODIFICAÇÃO: Abrir o novo arquivo com o histórico avançado
with open("history_advanced.csv", newline='', encoding="utf-8") as csvfile:
    reader = csv.reader(csvfile)
    
    # MODIFICAÇÃO: Pular a linha do cabeçalho
    next(reader, None) 
    
    for row in reader:
        # MODIFICAÇÃO: Ajustar os índices das colunas
        geracoes.append(int(row[0]))
        fitness_global.append(float(row[1]))
        fitness_geracional.append(float(row[2]))
        
        # Leitura dos genes (convertendo para inteiros)
        # Os genes começam na coluna 3 e vão até o final
        genes = [int(g) for g in row[3:]]
        genes_history.append(genes)

# --- Geração do Gráfico de Fitness (sem alterações na lógica) ---
janela_media_movel = 20
fitness_geracional_suave = media_movel(fitness_geracional, janela_media_movel)
geracoes_suave = geracoes[:len(fitness_geracional_suave)]

plt.style.use('seaborn-v0_8-whitegrid') # Estilo mais moderno
plt.figure(figsize=(12, 7))

plt.plot(geracoes, fitness_global, color='dodgerblue', linewidth=2.5, marker='o', markersize=4,
         label='Melhor Fitness Global Acumulado')
plt.plot(geracoes_suave, fitness_geracional_suave, color='tomato', linewidth=2,
         label=f'Melhor Fitness Geracional (Média Móvel {janela_media_movel})', alpha=0.9)

plt.title('Evolução do Fitness ao Longo das Gerações', fontsize=16, fontweight='bold', pad=20)
plt.xlabel('Geração', fontsize=12)
plt.ylabel('Fitness', fontsize=12)
plt.gca().xaxis.set_major_locator(MaxNLocator(integer=True))
plt.legend(fontsize=11)
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.tight_layout()
plt.show()
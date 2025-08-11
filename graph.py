import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import csv

def media_movel(dados, janela=10):
    return [sum(dados[i:i+janela])/janela for i in range(len(dados)-janela+1)]

geracoes = []
fitness_global = []
fitness_geracional = []

with open("history.csv", newline='') as csvfile:
    reader = csv.reader(csvfile)
    next(reader)
    for row in reader:
        geracoes.append(int(row[0]))
        fitness_global.append(float(row[1]))
        fitness_geracional.append(float(row[2]))

janela = 20
fitness_geracional_suave = media_movel(fitness_geracional, janela)
geracoes_suave = geracoes[:len(fitness_geracional_suave)]

plt.figure(figsize=(12, 7))
plt.plot(geracoes, fitness_global, color='blue', linewidth=2, marker='o', markersize=4, label='Melhor Fitness Global')
plt.plot(geracoes_suave, fitness_geracional_suave, color='red', linewidth=1.5, label=f'Melhor Fitness Geracional (média {janela})', alpha=0.8)

plt.title('Evolução do Fitness ao Longo das Gerações', fontsize=14, fontweight='bold')
plt.xlabel('Geração', fontsize=12)
plt.ylabel('Fitness', fontsize=12)
plt.gca().xaxis.set_major_locator(MaxNLocator(integer=True))
plt.legend(fontsize=10)
plt.grid(True, linestyle='--', alpha=0.6)
plt.tight_layout()
plt.show()

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import csv

geracoes = []
fitness = []

# Leitura dos dados do arquivo CSV
with open("history.csv", newline='') as csvfile:
    reader = csv.reader(csvfile)
    for row in reader:
        geracoes.append(int(row[0]))
        fitness.append(float(row[1]))

# Plot scatter simples
plt.figure(figsize=(10, 6))
plt.scatter(geracoes, fitness, color='blue', label='AG Padrão')
plt.title('Evolução do Fitness (5x5)')
plt.xlabel('Geração')
plt.ylabel('Melhor Fitness')
plt.gca().xaxis.set_major_locator(MaxNLocator(integer=True))
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

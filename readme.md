# Projeto em NetBalance

Este projeto utiliza bibliotecas científicas do Python para análises numéricas e visualizações gráficas, rodando em ambiente interativo com JupyterLab.

## 📦 Requisitos

- Python 3.8 ou superior
- pip (gerenciador de pacotes)
- virtualenv (recomendado)

## 📁 Instalação

Clone o repositório (ou navegue até a pasta do projeto) e siga os passos:

### 1. Compile o projeto com os dados desejados:

```bash
gcc -o genetic_solver genetic_solver.c
```

### 2. Execute o algoritmo genético de otimização:

```bash
./genetic_solver 
```

### 3. Crie e ative um ambiente virtual (opcional, mas recomendado):

```bash
python -m venv venv
source venv/bin/activate  # Linux/macOS
venv\Scripts\activate     # Windows
```

### 4. Instale as dependências:
```bash
pip install -r requirements.txt
```

### 5. Inicie o ambiente gráfico:
```bash
python ./graph.py
```

O navegador será aberto automaticamente com a interface do JupyterLab, onde você poderá rodar os notebooks .ipynb do projeto.

## 🛠️ Execução
Abra o arquivo .ipynb principal no JupyterLab e execute célula por célula.

## 🔒 Licença
Este projeto é livre para uso educacional. Para outros usos, entre em contato.
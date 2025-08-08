/*
 * Algoritmo Genético para resolver sistema linear e otimizar fitness em C
 * Adaptado de implementação em Python, sem uso de objetos.
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define N 20
#define POP_SIZE 20
#define GEN 2000
#define MU_TAX_BASE 0.05
#define TOURNAMENT_SIZE 2

int initial_positions[N][N] = {
    {1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1}
};

int min_matrix[N][N] = {
    {80,60,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {60,80,60,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,60,80,60,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,60,80,60,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,60,80,60,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,60,80,60,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,60,80,60,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,60,80,60,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,60,80,60,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,60,80,60,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,60,80,60,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,60,80,60,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,60,80,60,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,60,80,60,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,60,80,60,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,60,80,60,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,60,80,60,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,60,80,60,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,60,80,60},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,60,80}
};

int max_matrix[N][N] = {
    {100,80,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {80,100,80,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,80,100,80,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,80,100,80,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,80,100,80,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,80,100,80,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,80,100,80,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,80,100,80,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,80,100,80,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,80,100,80,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,80,100,80,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,80,100,80,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,80,100,80,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,80,100,80,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,80,100,80,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,80,100,80,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,80,100,80,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,80,100,80,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,80,100,80},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,80,100}
};

double b_vector[N] = {
    100,100,100,100,100,100,100,100,100,100,
    100,100,100,100,100,100,100,100,100,100
};


double history[GEN];

/* Função para arredondar ao mínimo míltiplo de 10 */
double round_to_nearest_5(int value) {
    return round(value/10.0)*10.0;
}

/* Gera matriz tester com valores aleatórios dentro do pool */
void generate_tester(double tester[N][N]) {
    int i, j;
    int pool_size = 0;
    int *pool = malloc(N*N*10 * sizeof(int));

    /* Preenche pool de valores */
    for (i = 0; i < N; ++i) {
        for (j = i; j < N; ++j) {
            if (max_matrix[i][j] > 0) {
                double minv = round_to_nearest_5(min_matrix[i][j]);
                double maxv = round_to_nearest_5(max_matrix[i][j]);
                for (int v = (int)minv; v <= (int)maxv; v += 10) {
                    pool[pool_size++] = v;
                }
            }
        }
    }

    /* Atribui valor aleatório de pool respeitando limites */
    for (i = 0; i < N; ++i) {
        for (j = i; j < N; ++j) {
            if (max_matrix[i][j] > 0) {
                int val;
                do {
                    val = pool[rand() % pool_size];
                } while (val < min_matrix[i][j] || val > max_matrix[i][j]);
                tester[i][j] = tester[j][i] = val;
            } else {
                tester[i][j] = tester[j][i] = 0.0;
            }
        }
    }
    free(pool);
}

/* Solução de sistema linear por eliminação de Gauss */
int solve_linear(double A[N][N], double b[N], double x[N]) {
    int i, j, k;
    double ratio;
    double aug[N][N+1];

    /* Monta matriz aumentada */
    for (i = 0; i < N; ++i) {
        for (j = 0; j < N; ++j) aug[i][j] = A[i][j];
        aug[i][N] = b[i];
    }

    /* Eliminação */
    for (i = 0; i < N; ++i) {
        if (fabs(aug[i][i]) < 1e-12) return 0; /* singular */
        for (k = i+1; k < N; ++k) {
            ratio = aug[k][i] / aug[i][i];
            for (j = i; j <= N; ++j) {
                aug[k][j] -= ratio * aug[i][j];
            }
        }
    }
    /* Retrossubstituição */
    for (i = N-1; i >= 0; --i) {
        x[i] = aug[i][N];
        for (j = i+1; j < N; ++j) {
            x[i] -= aug[i][j] * x[j];
        }
        x[i] /= aug[i][i];
    }
    return 1;
}

/* Calcula fitness de um indivíduo */
double fitness(int positions[N][N]) {
    float final_point = 0.0;
    for (size_t f = 0; f < 20; f++)
    {
        int i, j;
        double tester[N][N];
        generate_tester(tester);
        double A[N][N];
        for (i = 0; i < N; ++i)
            for (j = 0; j < N; ++j)
                A[i][j] = positions[i][j] * tester[i][j];
        double sol[N];
        if (!solve_linear(A, b_vector, sol)) return 0.0;
        float total_abs = 0.0;
        for (i = 0; i < N; ++i) total_abs += fabs(sol[i]);
        if (total_abs == 0.0) return 0.0;
        final_point+= 10*N / total_abs;
    }
    return final_point;
}

/* Cópias, mutação, cruzamento e randomização */
void copy_positions(int src[N][N], int dst[N][N]) {
    for (int i=0;i<N;++i) for (int j=0;j<N;++j) dst[i][j] = src[i][j];
}

void randomize(int dst[N][N]) {
    int i, j;
    for (i=0; i<N; ++i)
        for (j=0; j<N; ++j)
            dst[i][j] = (i==j ? 1 : 0);
    for (i=0; i<N; ++i) for (j=i+1;j<N;++j) {
        int bit = rand() % 2;
        dst[i][j] = dst[j][i] = bit;
    }
}

void mutate(int src[N][N], int dst[N][N], double mu) {
    copy_positions(src, dst);
    for (int i=0; i<N; ++i) for (int j=i+1; j<N; ++j) {
        if ((double)rand()/RAND_MAX < mu) {
            dst[i][j] = dst[j][i] = 1 - dst[i][j];
        }
    }
}

void cross(int p1[N][N], int p2[N][N], int dst[N][N]) {
    copy_positions(p1, dst);
    for (int i=0; i<N; ++i) for (int j=i+1; j<N; ++j) {
        if ((double)rand()/RAND_MAX < 0.3) {
            dst[i][j] = dst[j][i] = p2[i][j];
        }
    }
}

/* Seleção por torneio: retorna índice do melhor pai */
int select_parent(int pop, double fitnesses[POP_SIZE]) {
    int best_idx = rand() % pop;
    double best_fit = fitnesses[best_idx];
    for (int i=1; i<TOURNAMENT_SIZE; ++i) {
        int idx = rand() % pop;
        if (fitnesses[idx] > best_fit) {
            best_fit = fitnesses[idx];
            best_idx = idx;
        }
    }
    return best_idx;
}

int main() {
    srand(time(NULL));
    int population[POP_SIZE][N][N];
    int new_pop[POP_SIZE][N][N];
    double fitnesses[POP_SIZE];
    double mu = MU_TAX_BASE;
    int gens_no_improve = 0;

    float history[GEN];

    /* Inicializa população */
    copy_positions(initial_positions, population[0]);
    for (int i=1; i<POP_SIZE; ++i) randomize(population[i]);

    /* Fitness inicial e solução do sistema tester = max_matrix */
    printf("Inicial:\n");
    for (int i=0;i<N;++i) {
        for (int j=0;j<N;++j) printf("%d ", initial_positions[i][j]);
        printf("\n");
    }
    double A0[N][N], sol0[N];
    for (int i=0; i<N; ++i) for (int j=0; j<N; ++j)
        A0[i][j] = initial_positions[i][j] * max_matrix[i][j];
    if (solve_linear(A0, b_vector, sol0)) {
        printf("Solucao tester = max_matrix:\n");
        for (int i=0;i<N;++i) printf("== %f\n", sol0[i]);
    } else printf("Erro ao resolver sistema inicial\n");
    double fit0 = fitness(initial_positions);
    printf("Fitness inicial: %f\n", fit0);

    /* Loop principal */
    int best_idx = 0;
    double best_fit = -INFINITY;
    int best_positions[N][N];

    for (int gen=0; gen<GEN; ++gen) {
        /* Avalia fitness */
        for (int i=0; i<POP_SIZE; ++i) {
            fitnesses[i] = fitness(population[i]);
        }
        /* Encontra melhor atual */
        int curr_best = 0;
        for (int i=1; i<POP_SIZE; ++i)
            if (fitnesses[i] > fitnesses[curr_best]) curr_best = i;
        double curr_fit = fitnesses[curr_best];

        if (curr_fit > best_fit) {
            best_fit = curr_fit;
            best_idx = curr_best;
            copy_positions(population[best_idx], best_positions);
            gens_no_improve = 0;
            mu = MU_TAX_BASE;
        } else {
            gens_no_improve++;
            if (gens_no_improve % 50 == 0) mu += 0.025;
        }
        history[gen] = best_fit;
        history[gen] = curr_fit;

        /* Gera nova população */
        int cnt = 0, attempts = 0;
        while (cnt < POP_SIZE && attempts < 10*POP_SIZE) {
            int p1 = select_parent(POP_SIZE, fitnesses);
            int p2 = select_parent(POP_SIZE, fitnesses);
            int child[N][N];
            cross(population[p1], population[p2], child);
            int mutated[N][N];
            mutate(child, mutated, mu);
            double f1 = fitness(population[p1]);
            double f2 = fitness(population[p2]);
            double fc = fitness(mutated);
            if (fc > f1 && fc > f2) {
                copy_positions(mutated, new_pop[cnt++]);
            }
            attempts++;
        }
        while (cnt < POP_SIZE) {
            copy_positions(best_positions, new_pop[cnt++]);
        }
        /* Swap population */
        for (int i=0;i<POP_SIZE;++i)
            copy_positions(new_pop[i], population[i]);
    }

    printf("\nMelhor individuo encontrado (fitness = %f):\n", best_fit);
    for (int i=0;i<N;++i) {
        for (int j=0;j<N;++j) printf("%d ", best_positions[i][j]);
        printf("\n");
    }
    FILE *f = fopen("history.csv", "w");
    if (!f) {
        perror("Erro ao abrir arquivo");
        return 1;
    }

    for (int gen = 0; gen < GEN; gen++) {
        fprintf(f, "%d,%.6f\n", gen, history[gen]);
    }
    fclose(f);

    return 0;
}
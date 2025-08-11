#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define N 20
#define POP_SIZE 50
#define GEN 6000
#define MU_TAX_BASE 0.05
#define TOURNAMENT_SIZE 50

#define EVAL_MATRICES 50 
#define EVAL_LOOPS 50

#define REGEN_INTERVAL 20

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

double round_to_nearest_5(int value) {
    return round(value/10.0)*10.0;
}

void generate_tester(double tester[N][N]) {
    int i, j;
    int pool_size = 0;
    int *pool = malloc(N*N*10 * sizeof(int));

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

int solve_linear(double A[N][N], double b[N], double x[N]) {
    int i, j, k;
    double ratio;
    double aug[N][N+1];

    for (i = 0; i < N; ++i) {
        for (j = 0; j < N; ++j) aug[i][j] = A[i][j];
        aug[i][N] = b[i];
    }

    for (i = 0; i < N; ++i) {
        if (fabs(aug[i][i]) < 1e-12) return 0;
        for (k = i+1; k < N; ++k) {
            ratio = aug[k][i] / aug[i][i];
            for (j = i; j <= N; ++j) {
                aug[k][j] -= ratio * aug[i][j];
            }
        }
    }
    for (i = N-1; i >= 0; --i) {
        x[i] = aug[i][N];
        for (j = i+1; j < N; ++j) {
            x[i] -= aug[i][j] * x[j];
        }
        x[i] /= aug[i][i];
    }
    return 1;
}
double fitness(int positions[N][N], double testers[][N][N], int num_testers, int loops) {
    long double final_point = 0.0;
    for (size_t f = 0; f < (size_t)loops; f++)
    {
        int i, j;
        int idx = f % num_testers; /* consulta rotativa à tabela pré-gerada */
        /* uso direto testers[idx][i][j] */
        double A[N][N];
        for (i = 0; i < N; ++i)
            for (j = 0; j < N; ++j)
                A[i][j] = positions[i][j] * testers[idx][i][j];
        double sol[N];
        if (!solve_linear(A, b_vector, sol)) continue;
        long double total_abs = 0.0;
        for (i = 0; i < N; ++i) total_abs += fabsl(sol[i]);
        if (total_abs == 0.0) continue;
        final_point += 1000 - total_abs;
    }
    return (double)(final_point/100.0);
}

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

int select_parent(int pop, double fitnesses[POP_SIZE], int exclude) {
    int best_idx = rand() % pop;
    double best_fit = fitnesses[best_idx];
    for (int i=1; i<TOURNAMENT_SIZE; ++i) {
        int idx = rand() % pop;
        if (fitnesses[idx] > best_fit && idx != exclude) {
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

    float history_global[GEN];
    float history_generation[GEN];

    double evaluation_matrices[EVAL_MATRICES][N][N];
    for (int t = 0; t < EVAL_MATRICES; ++t) {
        generate_tester(evaluation_matrices[t]);
    }

    copy_positions(initial_positions, population[0]);
    for (int i=1; i<POP_SIZE; ++i) randomize(population[i]);

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
    double fit0 = fitness(initial_positions, evaluation_matrices, EVAL_MATRICES, EVAL_LOOPS);
    printf("Fitness inicial: %f\n", fit0);

    int best_idx = 0;
    double best_fit_global = -INFINITY;
    int best_positions[N][N];

    clock_t start_time = clock();

    for (int gen=0; gen<GEN; ++gen) {

        if (gen > 0 && (gen % REGEN_INTERVAL) == 0) {
            for (int t = 0; t < EVAL_MATRICES; ++t) {
                generate_tester(evaluation_matrices[t]);
            }
            printf("[geracao %d] Regeneradas %d evaluation_matrices (interval = %d)\n",
                   gen, EVAL_MATRICES, REGEN_INTERVAL);
        }

        for (int i=0; i<POP_SIZE; ++i) {
            fitnesses[i] = fitness(population[i], evaluation_matrices, EVAL_MATRICES, EVAL_LOOPS);
        }
        int curr_best_idx = 0;
        for (int i=1; i<POP_SIZE; ++i)
            if (fitnesses[i] > fitnesses[curr_best_idx]) curr_best_idx = i;
        double best_fit_generation = fitnesses[curr_best_idx];

        if (best_fit_generation > best_fit_global) {
            best_fit_global = best_fit_generation;
            best_idx = curr_best_idx;
            copy_positions(population[best_idx], best_positions);
            gens_no_improve = 0;
            mu = MU_TAX_BASE;
        } else {
            gens_no_improve++;
            if (gens_no_improve % 50 == 0) mu += 0.025;
        }

        history_global[gen] = best_fit_global;
        history_generation[gen] = best_fit_generation;

        int cnt = 0, attempts = 0;
        while (cnt < POP_SIZE && attempts < 10*POP_SIZE) {
            int p1 = select_parent(POP_SIZE, fitnesses, POP_SIZE+1);
            int p2 = select_parent(POP_SIZE, fitnesses, p1);
            int child[N][N];
            cross(population[p1], population[p2], child);
            int mutated[N][N];
            mutate(child, mutated, mu);
            double f1 = fitness(population[p1], evaluation_matrices, EVAL_MATRICES, EVAL_LOOPS);
            double f2 = fitness(population[p2], evaluation_matrices, EVAL_MATRICES, EVAL_LOOPS);
            double fc = fitness(mutated, evaluation_matrices, EVAL_MATRICES, EVAL_LOOPS);
            if (fc > f1 && fc > f2) {
                copy_positions(mutated, new_pop[cnt++]);
            }
            attempts++;
        }
        while (cnt < POP_SIZE) {
            copy_positions(best_positions, new_pop[cnt++]);
        }
        for (int i=0;i<POP_SIZE;++i)
            copy_positions(new_pop[i], population[i]);
    }

    clock_t end_time = clock();
    double elapsed_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;

    printf("\nMelhor individuo encontrado (fitness = %f):\n", best_fit_global);
    for (int i=0;i<N;++i) {
        for (int j=0;j<N;++j) printf("%d ", best_positions[i][j]);
        printf("\n");
    }

    printf("Tempo total de execução: %.3f segundos\n", elapsed_time);

    FILE *f = fopen("history.csv", "w");
    if (!f) {
        perror("Erro ao abrir arquivo");
        return 1;
    }

    fprintf(f, "geracao,melhor_fitness_global,melhor_fitness_geracional\n");
    for (int gen = 0; gen < GEN; gen++) {
        fprintf(f, "%d,%.6f,%.6f\n", gen, history_global[gen], history_generation[gen]);
    }
    fclose(f);

    return 0;
}
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#define N 20
#define POP_SIZE 50
#define GEN 1000
#define MU_TAX_BASE 0.05
#define TOURNAMENT_SIZE 10

#define EVAL_MATRICES 10
#define EVAL_LOOPS 10
#define REGEN_INTERVAL 500

typedef int mati[N][N];
typedef double matd[N][N];

typedef struct {
    double global_best_fit;
    double generation_best_fit;
    mati genes;
} HistoryEntry;

static int initial_positions[N][N] = {
    {1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1}
};

static int min_matrix[N][N] = {
    {100,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30},
    {30,100,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30},
    {30,30,100,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30},
    {30,30,30,100,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30},
    {30,30,30,30,100,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30},
    {30,30,30,30,30,100,30,30,30,30,30,30,30,30,30,30,30,30,30,30},
    {30,30,30,30,30,30,100,30,30,30,30,30,30,30,30,30,30,30,30,30},
    {30,30,30,30,30,30,30,100,30,30,30,30,30,30,30,30,30,30,30,30},
    {30,30,30,30,30,30,30,30,100,30,30,30,30,30,30,30,30,30,30,30},
    {30,30,30,30,30,30,30,30,30,100,30,30,30,30,30,30,30,30,30,30},
    {30,30,30,30,30,30,30,30,30,30,100,30,30,30,30,30,30,30,30,30},
    {30,30,30,30,30,30,30,30,30,30,30,100,30,30,30,30,30,30,30,30},
    {30,30,30,30,30,30,30,30,30,30,30,30,100,30,30,30,30,30,30,30},
    {30,30,30,30,30,30,30,30,30,30,30,30,30,100,30,30,30,30,30,30},
    {30,30,30,30,30,30,30,30,30,30,30,30,30,30,100,30,30,30,30,30},
    {30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,100,30,30,30,30},
    {30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,100,30,30,30},
    {30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,100,30,30},
    {30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,100,30},
    {30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,100}
};

static int max_matrix[N][N] = {
    {100,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80},
    {80,100,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80},
    {80,80,100,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80},
    {80,80,80,100,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80},
    {80,80,80,80,100,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80},
    {80,80,80,80,80,100,80,80,80,80,80,80,80,80,80,80,80,80,80,80},
    {80,80,80,80,80,80,100,80,80,80,80,80,80,80,80,80,80,80,80,80},
    {80,80,80,80,80,80,80,100,80,80,80,80,80,80,80,80,80,80,80,80},
    {80,80,80,80,80,80,80,80,100,80,80,80,80,80,80,80,80,80,80,80},
    {80,80,80,80,80,80,80,80,80,100,80,80,80,80,80,80,80,80,80,80},
    {80,80,80,80,80,80,80,80,80,80,100,80,80,80,80,80,80,80,80,80},
    {80,80,80,80,80,80,80,80,80,80,80,100,80,80,80,80,80,80,80,80},
    {80,80,80,80,80,80,80,80,80,80,80,80,100,80,80,80,80,80,80,80},
    {80,80,80,80,80,80,80,80,80,80,80,80,80,100,80,80,80,80,80,80},
    {80,80,80,80,80,80,80,80,80,80,80,80,80,80,100,80,80,80,80,80},
    {80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,100,80,80,80,80},
    {80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,100,80,80,80},
    {80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,100,80,80},
    {80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,100,80},
    {80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,100}
};

static double b_vector[N] = {100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100};

static HistoryEntry history[GEN];

void save_tester_config() {
    FILE *f = fopen("tester.csv", "w");
    if (!f) {
        perror("Erro ao abrir tester.csv");
        return;
    }

    fprintf(f, "i,j,min_value,max_value\n");

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            fprintf(f, "%d,%d,%d,%d\n", i, j, min_matrix[i][j], max_matrix[i][j]);
        }
    }
    fclose(f);
}


int solve_linear_fallback(const double A_in[N][N], const double b_in[N], double x_out[N]) {
    double aug[N][N+1];
    for (int i=0;i<N;++i) {
        for (int j=0;j<N;++j) aug[i][j]=A_in[i][j];
        aug[i][N]=b_in[i];
    }
    for (int i=0;i<N;++i) {
        int pivot=i; double maxv=fabs(aug[i][i]);
        for (int r=i+1;r<N;++r){ double av=fabs(aug[r][i]); if (av>maxv){ maxv=av; pivot=r; } }
        if (maxv < 1e-12) return 0;
        if (pivot != i) for (int c=i;c<=N;++c){ double t=aug[i][c]; aug[i][c]=aug[pivot][c]; aug[pivot][c]=t; }
        double diag = aug[i][i];
        for (int k=i+1;k<N;++k){
            double f = aug[k][i] / diag;
            if (f==0.0) continue;
            for (int j=i;j<=N;++j) aug[k][j] -= f * aug[i][j];
        }
    }
    for (int i=N-1;i>=0;--i){
        double s = aug[i][N];
        for (int j=i+1;j<N;++j) s -= aug[i][j] * x_out[j];
        x_out[i] = s / aug[i][i];
    }
    return 1;
}

int solve_linear_lapack_with_cache(const double A_in[N][N], const double b_in[N], double x_out[N]) {
    return solve_linear_fallback(A_in, b_in, x_out);
}

static inline double drand_01(void){ return rand() / (double)RAND_MAX; }
static inline void copy_positions(const int src[N][N], int dst[N][N]) { memcpy(dst, src, N*N*sizeof(int)); }

void generate_tester(double tester[N][N]) {
    for (int i=0;i<N;++i) for (int j=0;j<=i;++j) {
        int minv = min_matrix[i][j];
        int maxv = max_matrix[i][j];
        int v = 0;
        if (maxv > 0 && maxv >= minv) {
            int steps = (maxv - minv) / 10 + 1;
            v = minv + 10 * (rand() % steps);
        }
        tester[i][j] = tester[j][i] = (double)v;
    }
}

double fitness(const int positions[N][N], double testers[][N][N], int num_testers, int loops) {
    double final_point = 0.0;
    for (int f=0; f<loops; ++f) {
        int idx = f % num_testers;
        double A[N][N];
        for (int i=0;i<N;++i) for (int j=0;j<N;++j) A[i][j] = positions[i][j] ? testers[idx][i][j] : 0.0;
        double sol[N] = {0};
        if (!solve_linear_lapack_with_cache(A, b_vector, sol)) continue;
        double total_abs = 0.0;
        for (int i=0;i<N;++i) total_abs += fabs(sol[i]);
        if (total_abs == 0.0) continue;
        final_point += 1.0 / total_abs;
    }
    return final_point*10;
}

void randomize(int p[N][N]) {
    for (int i=0;i<N;++i) for (int j=0;j<=i;++j) {
        int v = rand() & 1;
        if(max_matrix[i][j] == 0){
            v = 0;
        }else if(min_matrix[i][j] == max_matrix[i][j]){
            v = 1;
        }
        p[i][j] = p[j][i] = v;
    }
}

void mutate(const int src[N][N], int dst[N][N], double mu) {
    copy_positions(src, dst);
    for (int i=0; i<N; ++i) {
        for (int j=i+1; j<N; ++j) {
            if(max_matrix[i][j] == 0){
                dst[i][j] = dst[j][i] = 0;
            }else if(min_matrix[i][j] == max_matrix[i][j]){
                dst[i][j] = dst[j][i] = 1;
            }else if (drand_01() < mu) {
                dst[i][j] = dst[j][i] = 1 - dst[i][j];
            }
        }
    }
}

void cross(const int p1[N][N], const int p2[N][N], int dst[N][N]) {
    copy_positions(p1, dst);
    for (int i=0; i<N; ++i) {
        for (int j=i+1; j<N; ++j) {
            if (drand_01() < 0.3) {
                dst[i][j] = dst[j][i] = p2[i][j];
            }
        }
    }
}

int select_parent(int pop, const double fitnesses[POP_SIZE], int exclude) {
    int best_idx = rand() % pop;
    double best_fit = fitnesses[best_idx];
    int iterations = (TOURNAMENT_SIZE < pop) ? TOURNAMENT_SIZE : pop;
    for (int k=1;k<iterations;++k) {
        int idx = rand() % pop;
        if (idx == exclude) continue;
        if (fitnesses[idx] > best_fit) { best_fit = fitnesses[idx]; best_idx = idx; }
    }
    return best_idx;
}

int main(void) {
    srand((unsigned)time(NULL));

    static int population[POP_SIZE][N][N];
    static int new_pop[POP_SIZE][N][N];
    static double fitnesses[POP_SIZE];
    static double evaluation_matrices[EVAL_MATRICES][N][N];

    double mu = MU_TAX_BASE;
    int gens_no_improve = 0;

    for (int t=0;t<EVAL_MATRICES;++t) generate_tester(evaluation_matrices[t]);
    
    // MODIFICAÇÃO: Salva a configuração inicial dos testers
    save_tester_config();
    printf("Configuração inicial dos testers salva em 'tester.csv'\n");

    copy_positions(initial_positions, population[0]);
    for (int i=1;i<POP_SIZE;++i) randomize(population[i]);

    double fit0 = fitness(initial_positions, evaluation_matrices, EVAL_MATRICES, EVAL_LOOPS);
    printf("Fitness inicial: %f\n", fit0);

    double best_fit_global = fit0;
    int best_positions[N][N];
    copy_positions(initial_positions, best_positions);

    clock_t start_time = clock();

    for (int gen=0; gen<GEN; ++gen) {
        if (gen > 0 && (gen % REGEN_INTERVAL) == 0) {
            for (int t=0;t<EVAL_MATRICES;++t) generate_tester(evaluation_matrices[t]);
            printf("[geracao %d] Regeneradas %d evaluation_matrices\n", gen, EVAL_MATRICES);
            
            // MODIFICAÇÃO: Salva a configuração dos testers novamente
            save_tester_config();
            printf("[geracao %d] Configuração dos testers salva novamente em 'tester.csv'\n", gen);
        }

        #pragma omp parallel for if(POP_SIZE>1)
        for (int i=0;i<POP_SIZE;++i) fitnesses[i] = fitness(population[i], evaluation_matrices, EVAL_MATRICES, EVAL_LOOPS);

        int curr_best_idx = 0;
        for (int i=1;i<POP_SIZE;++i) if (fitnesses[i] > fitnesses[curr_best_idx]) curr_best_idx = i;
        double best_fit_generation = fitnesses[curr_best_idx];

        if (best_fit_generation > best_fit_global) {
            best_fit_global = best_fit_generation;
            copy_positions(population[curr_best_idx], best_positions);
            gens_no_improve = 0;
            mu = MU_TAX_BASE;
        } else {
            gens_no_improve++;
            if (gens_no_improve % 50 == 0) mu += 0.025;
        }

        history[gen].global_best_fit = best_fit_global;
        history[gen].generation_best_fit = best_fit_generation;
        copy_positions(population[curr_best_idx], history[gen].genes);

        int cnt = 0, attempts = 0;
        while (cnt < POP_SIZE && attempts < 10 * POP_SIZE) {
            int p1 = select_parent(POP_SIZE, fitnesses, POP_SIZE + 1);
            int p2 = select_parent(POP_SIZE, fitnesses, p1);

            int child[N][N];
            cross(population[p1], population[p2], child);

            int mutated[N][N];
            mutate(child, mutated, mu);

            double f1 = fitnesses[p1];
            double f2 = fitnesses[p2];
            double fc = fitness(mutated, evaluation_matrices, EVAL_MATRICES, EVAL_LOOPS);

            if (fc > f1 && fc > f2) {
                copy_positions(mutated, new_pop[cnt++]);
            }
            attempts++;
        }

        while (cnt < POP_SIZE) {
            copy_positions(best_positions, new_pop[cnt++]);
        }

        for (int i=0;i<POP_SIZE;++i) copy_positions(new_pop[i], population[i]);

        if ((gen % 100) == 0) {
            printf("Geração %d, melhor geracional = %f, melhor global = %f\n", gen, best_fit_generation, best_fit_global);
        }
    }

    clock_t end_time = clock();
    double elapsed_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;

    printf("\nMelhor fitness global: %f\n", best_fit_global);
    printf("Melhor indivíduo encontrado:\n");
    for (int i=0;i<N;++i){
        for (int j=0;j<N;++j) printf("%d ", best_positions[i][j]);
        printf("\n");
    }
    printf("Tempo total de execução: %.3f segundos\n", elapsed_time);

    FILE *f = fopen("history_advanced.csv","w");
    if (!f) {
        perror("Erro ao abrir arquivo");
        return 1;
    }
    
    // Escreve o cabeçalho
    fprintf(f, "Generation,GlobalBestFitness,GenerationBestFitness");
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            fprintf(f, ",Gene_%d_%d", i, j);
        }
    }
    fprintf(f, "\n");

    // Escreve o estado inicial (Geração 0)
    fprintf(f, "0,%f,%f", fit0, fit0);
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            fprintf(f, ",%d", initial_positions[i][j]);
        }
    }
    fprintf(f, "\n");

    // Escreve os dados de cada geração da simulação (a partir da Geração 1)
    for (int gen = 0; gen < GEN; ++gen) {
        fprintf(f, "%d,%f,%f", gen + 1, history[gen].global_best_fit, history[gen].generation_best_fit);
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                fprintf(f, ",%d", history[gen].genes[i][j]);
            }
        }
        fprintf(f, "\n");
    }
    fclose(f);
    printf("\nHistórico avançado salvo em 'history_advanced.csv'\n");

    return 0;
}
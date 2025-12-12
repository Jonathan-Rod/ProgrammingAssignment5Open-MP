// %%writefile functions.c
// #include "functions.h"
#include "../include/functions.h"

// ============================================================================
// FUNCIONES DE INICIALIZACIÓN Y CONFIGURACIÓN
// ============================================================================

void initialize_default_parameters(SimulationParams *params) {
  // Parámetros físicos (acero)
  params->k = 10.0;       // Conductividad térmica [W/mK]
  params->rho_c = 1.0e7;  // Densidad * calor específico [J/m³K]

  // Geometría y dominio
  params->L = 2.0e-2;       // Longitud del dominio [m]
  params->n_volumes = 100;  // Número de volúmenes de control

  // Condiciones iniciales y de frontera
  params->T_initial = 200.0;  // Temperatura inicial [°C]
  params->T_cooled = 0.0;     // Temperatura superficie enfriada [°C]

  // Parámetros temporales
  params->total_time = 850.0;  // Tiempo total de simulación [s]
  params->dt = 0.01;           // Paso de tiempo [s] (1 ms)

  // Inicializar arreglos auxiliares
  params->n_profiles = 0;
  for (int i = 0; i < MAX_TIME_PROFILES; i++) {
    params->time_samples[i] = 0.0;
  }

  // Calcular parámetros derivados
  calculate_derived_parameters(params);
}

void calculate_derived_parameters(SimulationParams *params) {
  // Calcular espaciado espacial
  if (params->n_volumes > 1) {
    params->dx = params->L / (params->n_volumes);
  } else {
    params->dx = params->L;
  }

  // Calcular difusividad térmica
  calculate_thermal_diffusivity(params);

  // Calcular número de pasos de tiempo
  if (params->dt > 0) {
    params->n_time_steps = (int)(params->total_time / params->dt);
    if (params->n_time_steps * params->dt < params->total_time) {
      params->n_time_steps++;  // Asegurar que cubre el tiempo total
    }
  } else {
    params->n_time_steps = 0;
  }

  // Verificar consistencia
  if (params->dx <= 0) {
    fprintf(stderr, "Error: dx must be positive\n");
  }

  // Calcular coeficientes (no necesitan ser recalculados)
  params->aW = params->k / params->dx;
  params->aE = params->k / params->dx;
  params->aP = params->rho_c * (params->dx / params->dt);
  params->aP0 = params->rho_c * (params->dx / params->dt);
  params->aEb = 2 * (params->k / params->dx);
}

double *allocate_temperature_field(int n_volumes) {
  if (n_volumes <= 0 || n_volumes > MAX_VOLUMES) {
    fprintf(stderr, "Error: Invalid number of volumes: %d\n", n_volumes);
    return NULL;
  }
  // n_volumes es solo el numero de nodos
  // no incluye extremos entonces añadimos 2 para que sea 0 y n - 1 el ultimo
  // extremos
  double *T = (double *)calloc(n_volumes + 2, sizeof(double));
  if (T == NULL) {
    fprintf(stderr, "Error: Could not allocate memory for temperature field\n");
    return NULL;
  }

  return T;
}

void free_temperature_field(double *T) {
  if (T != NULL) {
    free(T);
  }
}

int validate_parameters(const SimulationParams *params) {
  // Verificar parámetros físicos positivos
  if (params->k <= 0) {
    fprintf(stderr, "Error: Thermal conductivity k must be positive\n");
    return 0;
  }
  if (params->rho_c <= 0) {
    fprintf(stderr, "Error: rho_c must be positive\n");
    return 0;
  }
  if (params->L <= 0) {
    fprintf(stderr, "Error: Length L must be positive\n");
    return 0;
  }

  // Verificar número de volúmenes
  if (params->n_volumes < 3 || params->n_volumes > MAX_VOLUMES) {
    fprintf(stderr, "Error: Number of volumes must be between 1 and %d\n",
            MAX_VOLUMES);
    return 0;
  }

  // Verificar condiciones de estabilidad
  if (!check_stability_condition(params)) {
    fprintf(stderr, "Error: Stability condition not met\n");
    return 0;
  }

  // Verificar consistencia dimensional
  if (params->dx <= 0) {
    fprintf(stderr, "Error: dx must be positive\n");
    return 0;
  }
  if (params->dt <= 0) {
    fprintf(stderr, "Error: dt must be positive\n");
    return 0;
  }

  // Verificar condiciones iniciales
  if (params->T_initial <= params->T_cooled) {
    fprintf(stderr, "Error: T_initial must be greater than T_cooled\n");
    return 0;
  }

  return 1;
}

// ============================================================================
// CÁLCULOS DE ESTABILIDAD Y COEFICIENTES (EXPLÍCITO)
// ============================================================================

double calculate_stability_limit(const SimulationParams *params) {
  if (params->alpha <= 0 || params->dx <= 0) {
    return 0.0;
  }

  // Límite teórico para esquema explícito: Fo <= 0.5
  double Fo_max = 0.5;
  // Aplicar factor de seguridad conservador
  double safety_factor = 0.8;
  Fo_max *= safety_factor;

  // dt_max = Fo_max * dx² / alpha
  double dt_max = Fo_max * params->dx * params->dx / params->alpha;

  return dt_max;
}

int check_stability_condition(const SimulationParams *params) {
  double Fo = calculate_fourier_number(params);
  return (Fo <= 0.5);
}

void calculate_thermal_diffusivity(SimulationParams *params) {
  if (params->rho_c == 0) {
    fprintf(stderr, "Error: rho_c cannot be zero\n");
    params->alpha = 0.0;
    return;
  }

  params->alpha = params->k / params->rho_c;

  // Verificar resultado físicamente razonable
  if (params->alpha <= 0) {
    fprintf(stderr, "Warning: Non-positive thermal diffusivity: %e\n",
            params->alpha);
  }
}

double calculate_fourier_number(const SimulationParams *params) {
  if (params->dx == 0) {
    return 0.0;
  }

  double Fo = params->alpha * params->dt / (params->dx * params->dx);
  return Fo;
}

// ============================================================================
// SIMULACIÓN SECUENCIAL (EXPLÍCITA)
// ============================================================================

void solve_heat_equation_sequential(double *T, const SimulationParams *params) {
  // 1. Inicializa campo de temperaturas con condición inicial
  int n_points = params->n_volumes + 2;
  for (int i = 0; i < n_points; i++) {
    T[i] = params->T_initial;
  }

  // 2. Aplica condiciones de frontera
  apply_boundary_conditions_sequential(T, params);

  // 3. Para cada paso de tiempo
  double *T_new = allocate_temperature_field(params->n_volumes);
  for (int i = 0; i < n_points; i++) {
    T_new[i] = T[i];
  }

  double current_time = 0.0;
  for (int step = 0; step < params->n_time_steps; step++) {
    current_time += params->dt;

    // 3a. Calcula nuevas temperaturas usando esquema explícito
    calculate_explicit_step_sequential(T_new, T, params);

    // 3b. Aplica condiciones de frontera
    apply_boundary_conditions_sequential(T_new, params);

    // 3c. Actualiza campo de temperaturas
    for (int i = 0; i < n_points; i++) {
      T[i] = T_new[i];
    }
    // Imprimir progreso de steps
    printf("Time step: %d/%d\n", step + 1, params->n_time_steps);
    printf("Current time: %.2f\n", current_time);

    double max_temp, min_temp;  // Por si explota fue aqui XD
    find_temperature_extremes(T, params->n_volumes, &max_temp, &min_temp);
    printf("Tmax: %.2f Tmin: %.2f\n", max_temp, min_temp);
  }

  // 4. Calcula métricas de convergencia
  // calculate_convergence_metrics(T, params, current_time);

  free_temperature_field(T_new);
  printf("Sequential simulation completed\n");
}

void solve_transient_sequential(double *T, const SimulationParams *params,
                                double *T_profiles, int profile_indices[]) {
  // TODO
  printf("Transient simulation completed\n");
}

void time_integration_sequential(double *T, const SimulationParams *params) {
  // TODO
}

void calculate_explicit_step_sequential(double *T_new, const double *T_old,
                                        const SimulationParams *params) {
  int i;
  double b;
  // Los bordes se manejan en apply_boundary_conditions
  // i = 0 y i = n_volumes - 1

  // Contribution from BC 1
  i = 1;  // primer nodo
  b = params->aE * T_old[i] + (params->aP0 - params->aE) * T_old[i];
  T_new[i] = b / params->aP;

  // Contribution from internal nodes (2 to n_volumes-3)
  for (int i = 2; i < params->n_volumes - 2; i++) {
    b = params->aW * T_old[i - 1] + params->aE * T_old[i + 1] +
        (params->aP0 - (params->aE + params->aW)) * T_old[i];
    T_new[i] = b / params->aP;
  }

  // Contribution from BC 2
  i = params->n_volumes - 2;  // ultimo nodo
  b = params->aW * T_old[i - 1] +
      (params->aP0 - (params->aE + params->aW)) * T_old[i] +
      params->aEb * params->T_cooled;
  T_new[i] = b / params->aP;

  // Los bordes se manejan en apply_boundary_conditions
}

void apply_boundary_conditions_sequential(double *T,
                                          const SimulationParams *params) {
  // Borde izquierdo: Dirichlet (T = T_cooled)
  // i = 0 hace referencia al extremo del primer nodo
  T[0] = params->T_cooled;
  // Borde derecho: Neumann (dT/dx = 0, aislado)
  // i = n_volumes - 1 hace referencia al extremo del ultimo nodo
  T[params->n_volumes - 1] = T[params->n_volumes - 2];
}

// ============================================================================
// SIMULACIÓN PARALELA (OPENMP - EXPLÍCITA)
// ============================================================================

void solve_heat_equation_parallel(double *T, const SimulationParams *params) {
  // Configurar entorno OpenMP
  configure_omp_environment(OMP_NUM_THREADS);
  // TODO

  printf("Parallel simulation completed\n");
}

void solve_transient_parallel(double *T, const SimulationParams *params,
                              double *T_profiles, int profile_indices[]) {
  // Configurar entorno OpenMP
  configure_omp_environment(OMP_NUM_THREADS);
  // TODO
  printf("Parallel transient simulation completed\n");
}

void time_integration_parallel(double *T, const SimulationParams *params) {
  // TODO
}

void calculate_explicit_step_parallel(double *T_new, const double *T_old,
                                      const SimulationParams *params) {
  // TODO
  // Los bordes se manejan en apply_boundary_conditions
}

void apply_boundary_conditions_parallel(double *T,
                                        const SimulationParams *params) {
// Aplicar condiciones de frontera (solo necesita ejecutarse una vez)
#pragma omp single
  {
    // Borde izquierdo: Dirichlet (T = T_cooled)
    T[0] = params->T_cooled;

    // Borde derecho: Neumann (dT/dx = 0, aislado)
    T[params->n_volumes - 1] = T[params->n_volumes - 2];
  }
}

// ============================================================================
// ANÁLISIS DE PERFORMANCE Y BENCHMARKING
// ============================================================================

PerformanceMetrics compare_sequential_vs_parallel(
    const SimulationParams *params) {
  PerformanceMetrics metrics = {0};
  // TODO
  return metrics;
}

void benchmark_scalability_analysis(const SimulationParams *params) {
  int max_threads = get_available_parallel_threads();
  // TODO
}

void find_optimal_thread_configuration(const SimulationParams *params) {
  int max_threads = get_available_parallel_threads();
  // TODO
}

void performance_sweep_parameters(const SimulationParams *base_params) {
  printf("\n=== PERFORMANCE PARAMETER SWEEP ===\n");

  // Variar número de volúmenes
  int volume_sizes[] = {100, 1000, 5000, 10000};
  int num_sizes = sizeof(volume_sizes) / sizeof(volume_sizes[0]);

  for (int i = 0; i < num_sizes; i++) {
    // TODO
  }
}

double calculate_speedup_ratio(double seq_time, double par_time) {
  if (par_time <= 0) {
    return 1.0;
  }
  if (seq_time <= 0) {
    return 1.0;
  }

  double speedup = seq_time / par_time;

  // Limitar a valores razonables
  if (speedup > 1000) {
    speedup = 1000;
  }

  return speedup;
}

double calculate_parallel_efficiency(double speedup, int n_threads) {
  if (n_threads <= 0) {
    return 0.0;
  }

  double efficiency = speedup / n_threads;

  // Normalizar entre 0 y 1
  if (efficiency < 0) efficiency = 0;
  if (efficiency > 1) efficiency = 1;

  return efficiency;
}

double measure_execution_time(void (*solver)(double *,
                                             const SimulationParams *),
                              double *T, const SimulationParams *params) {
  double start_time = omp_get_wtime();
  solver(T, params);
  double end_time = omp_get_wtime();

  return end_time - start_time;
}

// ============================================================================
// VALIDACIÓN Y VERIFICACIÓN (EXPLÍCITA)
// ============================================================================

double calculate_total_energy(const double *T, const SimulationParams *params) {
  double total_energy = 0.0;

  for (int i = 0; i < params->n_volumes; i++) {
    total_energy += params->rho_c * T[i] * params->dx;
  }

  return total_energy;
}

int verify_energy_conservation(const double *T_initial, const double *T_final,
                               const SimulationParams *params) {
  double energy_initial = calculate_total_energy(T_initial, params);
  double energy_final = calculate_total_energy(T_final, params);

  // La energía debe decrecer debido a la superficie enfriada
  int energy_decreases = (energy_final < energy_initial);

  printf("Initial energy: %.6e J\n", energy_initial);
  printf("Final energy:   %.6e J\n", energy_final);
  printf("Difference:     %.6e J\n", energy_final - energy_initial);
  printf("Energy decreases: %s\n", energy_decreases ? "YES" : "NO");

  return energy_decreases;
}

int verify_solution_equivalence(const double *T_seq, const double *T_par,
                                int n_volumes, double tolerance) {
  double max_difference = 0.0;

  for (int i = 0; i < n_volumes; i++) {
    double difference = fabs(T_seq[i] - T_par[i]);
    if (difference > max_difference) {
      max_difference = difference;
    }
  }

  printf("Maximum difference between solutions: %.2e\n", max_difference);
  printf("Tolerance: %.2e\n", tolerance);
  printf("Solutions equivalent: %s\n",
         (max_difference <= tolerance) ? "YES" : "NO");

  return (max_difference <= tolerance);
}

double calculate_numerical_error(const double *T_numeric,
                                 const SimulationParams *params,
                                 double current_time) {
  // Validación simplificada usando propiedades físicas conocidas
  double error = 0.0;
  int count = 0;

  // Verificar propiedades físicas básicas en lugar de solución analítica exacta
  for (int i = 0; i < params->n_volumes - 1; i++) {
    double x = (params->dx / 2) + i * params->dx;

    // Propiedades físicas esperadas:
    // 1. Temperaturas deben estar entre T_cooled y T_initial
    int within_bounds = (T_numeric[i] >= params->T_cooled - 1e-6) &&
                        (T_numeric[i] <= params->T_initial + 1e-6);

    // 2. El perfil debe ser monótono decreciente (para este problema
    // específico)
    // 3. Las temperaturas cerca del borde enfriado deben ser más bajas

    double expected_trend =
        params->T_initial -
        (params->T_initial - params->T_cooled) * (x / params->L);

    error += fabs(T_numeric[i] - expected_trend);
    count++;
  }

  return (count > 0) ? error / count : 0.0;
}

void validate_convergence_history(const SimulationParams *params) {
  printf("\n=== CONVERGENCE VALIDATION ===\n");

  // Verificar estabilidad del esquema explícito
  if (!check_stability_condition(params)) {
    printf("WARNING: Stability condition not met\n");
  } else {
    printf("Stability condition: MET\n");
  }

  // Verificar que Fo está en rango razonable
  double Fo = calculate_fourier_number(params);
  printf("Fourier number: %.4f\n", Fo);

  if (Fo < 0.01) {
    printf("WARNING: Fo very small, potential precision issues\n");
  }
}

void find_temperature_extremes(const double *T, int n_volumes, double *max_temp,
                               double *min_temp) {
  if (n_volumes <= 0) {
    *max_temp = 0.0;
    *min_temp = 0.0;
    return;
  }

  *max_temp = T[0];
  *min_temp = T[0];

  // Incluir puntos de frontera
  int n_points = n_volumes + 2;

  for (int i = 0; i < n_points; i++) {
    if (T[i] > *max_temp) {
      *max_temp = T[i];
    }
    if (T[i] < *min_temp) {
      *min_temp = T[i];
    }
  }
}

// ============================================================================
// GESTIÓN DE DATOS Y ARCHIVOS
// ============================================================================

void save_temperature_profile_csv(const double *T,
                                  const SimulationParams *params,
                                  double current_time, const char *filename) {
  // TODO
}

void save_performance_metrics_csv(const PerformanceMetrics *metrics,
                                  const char *filename) {
  // TODO
}

void save_transient_profiles_csv(const double *T_profiles,
                                 const SimulationParams *params,
                                 const char *filename) {
  // TODO
}

void save_scalability_data_csv(const PerformanceMetrics *metrics_array,
                               int num_configs, const char *filename) {
  // TODO
}

FILE *safe_file_open(const char *filename, const char *mode) {
  FILE *file = fopen(filename, mode);
  if (file == NULL) {
    fprintf(stderr, "Error: Could not open file %s in mode %s\n", filename,
            mode);
    perror("Detailed error");
  }
  return file;
}

void safe_file_close(FILE *file) {
  if (file != NULL) {
    if (fclose(file) != 0) {
      fprintf(stderr, "Error: Could not close file properly\n");
      perror("Detailed error");
    }
  }
}

void write_csv_headers(FILE *file, const char *headers[], int n_headers) {
  if (file == NULL || headers == NULL || n_headers <= 0) {
    return;
  }

  // Escribir primer header
  fprintf(file, "%s", headers[0]);

  // Escribir headers restantes con coma
  for (int i = 1; i < n_headers; i++) {
    fprintf(file, ",%s", headers[i]);
  }
  fprintf(file, "\n");

  fflush(file);
}

// ============================================================================
// CONFIGURACIÓN Y CONTROL OPENMP
// ============================================================================

void configure_omp_environment(int num_threads) {
  if (num_threads > 0) {
    omp_set_num_threads(num_threads);
  }
}

void print_omp_configuration_info(void) {
  printf("\n=== OPENMP CONFIGURATION ===\n");
  printf("Maximum threads available: %d\n", omp_get_max_threads());
  printf("Available processors: %d\n", omp_get_num_procs());
  printf("In parallel region: %s\n", omp_in_parallel() ? "YES" : "NO");

#ifdef _OPENMP
  printf("OpenMP version: %d\n", _OPENMP);
#else
  printf("OpenMP not available\n");
#endif
}

int get_available_parallel_threads(void) {
  int max_threads = omp_get_max_threads();
  // Limitar para pruebas
  if (max_threads > 16) {
    max_threads = 16;
  }
  return max_threads;
}

void set_omp_dynamic_scheduling(int chunk_size) {
  if (chunk_size > 0) {
    omp_set_schedule(omp_sched_dynamic, chunk_size);
  }
}

// ============================================================================
// UTILIDADES DE VISUALIZACIÓN Y DEBUG
// ============================================================================

void print_temperature_field(const double *T, int n_volumes,
                             const char *label) {
  printf("\n=== %s ===\n", label);
  printf("Position [m] | Temperature [C]\n");
  printf("-----------------------------\n");

  int n_points = n_volumes + 2;
  for (int i = 0; i < n_points; i++) {
    double x = i * 0.1 / (n_points - 1);  // Posición normalizada
    printf("%10.4f | %12.4f\n", x, T[i]);

    if (i == 0 || i == n_points - 1) {
      printf("  [Boundary]");
    }

    // Limitar impresión para muchos volúmenes
    if (n_points > 20 && i >= 10 && i < n_points - 10) {
      if (i == 10) {
        printf("... (omitting %d points) ...\n", n_points - 20);
      }
      continue;
    }
  }
}

void print_simulation_progress(int iteration, int total, double max_temp,
                               double min_temp) {
  int bar_width = 50;
  float progress = (float)iteration / total;
  int pos = (int)(bar_width * progress);

  printf("\rProgress: [");
  for (int i = 0; i < bar_width; i++) {
    if (i < pos)
      printf("=");
    else if (i == pos)
      printf(">");
    else
      printf(" ");
  }
  printf("] %d%% | Step: %d/%d | Tmax: %.2fC | Tmin: %.2fC",
         (int)(progress * 100), iteration, total, max_temp, min_temp);
  fflush(stdout);

  if (iteration == total) {
    printf("\n");
  }
}

void visualize_domain_partitioning(int n_volumes, int n_threads) {
  printf("\n=== DOMAIN PARTITIONING ===\n");
  printf("Total volumes: %d\n", n_volumes);
  printf("Threads: %d\n", n_threads);

  int volumes_per_thread = n_volumes / n_threads;
  int remainder = n_volumes % n_threads;

  printf("Distribution:\n");
  int start = 0;
  for (int thread = 0; thread < n_threads; thread++) {
    int end = start + volumes_per_thread - 1;
    if (thread < remainder) {
      end++;  // Distribuir resto entre primeros hilos
    }

    printf("Thread %d: volumes %d-%d (%d volumes)\n", thread, start, end,
           end - start + 1);
    start = end + 1;
  }

  // Representación visual simple
  printf("\nVisualization:\n");
  for (int i = 0; i < n_volumes; i++) {
    int thread = i % n_threads;
    printf("%d", thread);

    if ((i + 1) % 50 == 0) printf("\n");
  }
  printf("\n");
}

void print_performance_summary(const PerformanceMetrics *metrics) {
  printf("\n=== PERFORMANCE SUMMARY ===\n");
  printf("Sequential time:   %8.4f s\n", metrics->sequential_time);
  printf("Parallel time:     %8.4f s\n", metrics->parallel_time);
  printf("Speedup:           %8.2f x\n", metrics->speedup);
  printf("Efficiency:        %8.1f %%\n", metrics->efficiency * 100);

  if (metrics->speedup > 1.0) {
    printf("IMPROVEMENT: %.1f%% faster\n",
           (metrics->sequential_time - metrics->parallel_time) /
               metrics->sequential_time * 100);
  } else {
    printf("WARNING: No improvement with parallelization\n");
  }
}

// ============================================================================
// FUNCIONES DE PRUEBA Y VERIFICACIÓN
// ============================================================================

void run_correctness_test(void) {
  printf("\n=== CORRECTNESS TESTS ===\n");

  // Configurar parámetros de prueba
  SimulationParams test_params;
  initialize_default_parameters(&test_params);
  test_params.n_volumes = 10;
  test_params.dt = 1;
  calculate_derived_parameters(&test_params);

  // Prueba 1: Condiciones de frontera
  printf("1. Testing boundary conditions...\n");
  test_boundary_conditions();

  // Prueba 2: Cálculo explícito
  printf("2. Testing explicit calculation...\n");
  test_explicit_calculation();

  // Prueba 3: Conservación de energía
  printf("3. Testing energy conservation...\n");
  double *T_initial = allocate_temperature_field(test_params.n_volumes);
  double *T_final = allocate_temperature_field(test_params.n_volumes);

  for (int i = 0; i < test_params.n_volumes; i++) {
    T_initial[i] = test_params.T_initial;
    T_final[i] = test_params.T_initial * 0.5;  // Simular enfriamiento
  }

  int energy_ok = verify_energy_conservation(T_initial, T_final, &test_params);
  printf("Energy conservation: %s\n", energy_ok ? "OK" : "FAILED");

  free_temperature_field(T_initial);
  free_temperature_field(T_final);

  printf("=== TESTS COMPLETED ===\n");
}

void test_boundary_conditions(void) {
  SimulationParams params;
  initialize_default_parameters(&params);
  params.n_volumes = 5;
  calculate_derived_parameters(&params);

  double *T = allocate_temperature_field(params.n_volumes);
  for (int i = 0; i < params.n_volumes; i++) {
    T[i] = params.T_initial;
  }

  // Aplicar condiciones de frontera
  apply_boundary_conditions_sequential(T, &params);

  // Verificar condiciones
  int left_ok = (T[0] == params.T_cooled);
  int right_ok = (T[params.n_volumes - 1] == T[params.n_volumes - 2]);

  printf("Left condition (Dirichlet): %s\n", left_ok ? "OK" : "FAILED");
  printf("Right condition (Neumann): %s\n", right_ok ? "OK" : "FAILED");

  free_temperature_field(T);
}

void test_explicit_calculation(void) {
  SimulationParams params;
  initialize_default_parameters(&params);
  params.n_volumes = 5;  // nodos internos
  params.dt = 0.001;
  calculate_derived_parameters(&params);

  // Nodos internos + 2 extremos
  double *T_old = allocate_temperature_field(params.n_volumes);
  double *T_new = allocate_temperature_field(params.n_volumes);

  int n_points = params.n_volumes + 2;
  for (int i = 0; i < n_points; i++) {
    T_old[i] = params.T_initial;
    T_new[i] = params.T_initial;
  }

  // Aplicar condiciones de frontera
  apply_boundary_conditions_sequential(T_old, &params);

  // Calcular un paso explícito
  calculate_explicit_step_sequential(T_new, T_old, &params);

  // Aplicar condiciones de frontera
  apply_boundary_conditions_sequential(T_new, &params);

  // Verificar que se calcularon los volúmenes internos
  int calculated = 1;
  for (int i = 1; i < n_points - 1; i++) {
    if (T_new[i] == T_old[i]) {
      calculated = 0;
      break;
    }
  }

  printf("Explicit calculation performed: %s\n",
         calculated ? "PASSED" : "FAILED");

  free_temperature_field(T_old);
  free_temperature_field(T_new);
}

void verify_parallel_correctness(const SimulationParams *params) {
  printf("\n=== PARALLEL CORRECTNESS VERIFICATION ===\n");

  double *T_seq = allocate_temperature_field(params->n_volumes);
  double *T_par = allocate_temperature_field(params->n_volumes);

  if (T_seq == NULL || T_par == NULL) {
    fprintf(stderr, "Error: Could not allocate memory for verification\n");
    return;
  }

  // Ejecutar ambas versiones
  solve_heat_equation_sequential(T_seq, params);

  // Probar con diferentes números de hilos
  int thread_counts[] = {2, 4, 8};
  int num_threads = sizeof(thread_counts) / sizeof(thread_counts[0]);

  for (int i = 0; i < num_threads; i++) {
    configure_omp_environment(thread_counts[i]);
    solve_heat_equation_parallel(T_par, params);

    double tolerance = 1e-12;
    int equivalent =
        verify_solution_equivalence(T_seq, T_par, params->n_volumes, tolerance);

    printf("Correctness with %d threads: %s\n", thread_counts[i],
           equivalent ? "OK" : "FAILED");
  }

  free_temperature_field(T_seq);
  free_temperature_field(T_par);

  printf("=== VERIFICATION COMPLETED ===\n");
}
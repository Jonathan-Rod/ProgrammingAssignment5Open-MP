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
  params->n_profiles = 10;
  double profile_step = params->total_time / (params->n_profiles - 1);
  for (int i = 0; i < params->n_profiles; i++) {
    params->time_samples[i] = profile_step * i;
  }
  params->T_profiles =
      allocate_temperature_profiles(params->n_profiles, params->n_volumes);

  // Calcular parámetros derivados
  calculate_derived_parameters(params);
}

void calculate_derived_parameters(SimulationParams *params) {
  // Calcular espaciado espacial
  if (params->n_volumes > 0) {
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

double **allocate_temperature_profiles(int n_profiles, int n_volumes) {
  double **T_profiles = malloc(n_profiles * sizeof(double *));
  for (int i = 0; i < n_profiles; i++) {
    T_profiles[i] = allocate_temperature_field(n_volumes);
  }
  return T_profiles;
}

void free_temperature_field(double *T) {
  if (T != NULL) {
    free(T);
  }
}

void free_temperature_profiles(double **TT, int profiles) {
  for (int i = 0; i < profiles; i++) {
    free_temperature_field(TT[i]);
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
  int n_points = params->n_volumes + 2;

  // 1. Inicializar temperaturas
  for (int i = 0; i < n_points; i++) {
    T[i] = params->T_initial;
  }

  apply_boundary_conditions_sequential(T, params);

  // 2. Arreglo temporal
  double *T_temp = allocate_temperature_field(params->n_volumes);

  // 3. Bucle temporal
  for (int step = 0; step < params->n_time_steps; step++) {
    // Calcular nuevo paso basado en T
    calculate_explicit_step_sequential(T_temp, T, params);
    apply_boundary_conditions_sequential(T_temp, params);

    // Actualizar T
    for (int i = 0; i < n_points; i++) {
      T[i] = T_temp[i];
    }
  }
  free_temperature_field(T_temp);
}

void solve_transient_sequential(double *T, SimulationParams *params) {
  int n_points = params->n_volumes + 2;
  int original_n_steps = params->n_time_steps;

  // Inicializar temperaturas
  for (int i = 0; i < n_points; i++) {
    T[i] = params->T_initial;
  }
  apply_boundary_conditions_sequential(T, params);

  for (int i = 0; i < params->n_profiles; i++) {
    double target_time = params->time_samples[i];

    // Calcular pasos para este perfil
    int steps_needed = (int)(target_time / params->dt + 0.5);
    if (steps_needed > original_n_steps) {
      steps_needed = original_n_steps;
    }

    // Ejecutar simulación hasta este tiempo
    params->n_time_steps = steps_needed;
    solve_heat_equation_sequential(T, params);

    // Guardar perfil
    for (int j = 0; j < n_points; j++) {
      params->T_profiles[i][j] = T[j];
    }

    printf("Calculando perfil %d de %d en tiempo %.2f\n", i + 1,
           params->n_profiles, target_time);
  }

  // Restaurar valor original
  params->n_time_steps = original_n_steps;
}

void calculate_explicit_step_sequential(double *T_new, const double *T_old,
                                        const SimulationParams *params) {
  int i;
  double b;
  int n_points = params->n_volumes + 2;

  // Primer nodo interno (i=1) - Neumann izquierdo: T[0] = T[1]
  i = 1;
  b = params->aW * T_old[1] + params->aE * T_old[i + 1] +
      (params->aP0 - (params->aW + params->aE)) * T_old[i];
  T_new[i] = b / params->aP;

  // Nodos internos centrales (i=2 a n_volumes - 1)
  for (i = 2; i <= n_points - 3; i++) {
    b = params->aW * T_old[i - 1] + params->aE * T_old[i + 1] +
        (params->aP0 - (params->aE + params->aW)) * T_old[i];
    T_new[i] = b / params->aP;
  }

  // Último nodo interno (i=n_volumes+1)
  i = n_points - 2;
  b = params->aW * T_old[i - 1] + params->aEb * params->T_cooled +
      (params->aP0 - params->aW) * T_old[i];
  T_new[i - 1] = b / params->aP;
}

void apply_boundary_conditions_sequential(double *T,
                                          const SimulationParams *params) {
  int n_points = params->n_volumes + 2;

  // Borde izquierdo: Neumann (T[0] = T[1])
  T[0] = T[1];

  // Borde derecho: Dirichlet (T[n_points-1] = T_cooled)
  T[n_points - 1] = params->T_cooled;
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

void solve_transient_parallel(double *T, const SimulationParams *params) {
  // Configurar entorno OpenMP
  configure_omp_environment(OMP_NUM_THREADS);
  // TODO
  printf("Parallel transient simulation completed\n");
}

void calculate_explicit_step_parallel(double *T_new, const double *T_old,
                                      const SimulationParams *params) {
  // TODO
  // Los bordes se manejan en apply_boundary_conditions
}

void apply_boundary_conditions_parallel(double *T,
                                        const SimulationParams *params) {
  // TODO
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

void performance_sweep_parameters(const SimulationParams *base_params) {
  printf("\n=== PERFORMANCE PARAMETER SWEEP ===\n");

  // Variar número de profiles
  int profiles_numbers[] = {100, 1000, 5000, 10000};
  int profiles_size = sizeof(profiles_numbers) / sizeof(profiles_numbers[0]);

  for (int i = 0; i < profiles_size; i++) {
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

int verify_solution_equivalence(const double *T_seq, const double *T_par,
                                int n_volumes, double tolerance) {
  double max_difference = 0.0;

  for (int i = 0; i < n_volumes; i++) {
    double difference = fabs(T_seq[i] - T_par[i]);
    if (difference > max_difference) {
      max_difference = difference;
    }
  }
  printf("Maximum difference between solutions: %.2f\n", max_difference);
  printf("Tolerance: %.2f\n", tolerance);
  printf("Solutions equivalent: %s\n",
         max_difference <= tolerance ? "true" : "false");
  return (max_difference <= tolerance);
}

void validate_convergence(const SimulationParams *params) {
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

// ============================================================================
// GESTIÓN DE DATOS Y ARCHIVOS
// ============================================================================

void save_temperature_profile_csv(const double *T,
                                  const SimulationParams *params,
                                  double current_time, const char *filename) {
  if (T == NULL || params == NULL || filename == NULL) {
    fprintf(stderr, "Error: Invalid parameters\n");
    return;
  }

  int n_points = params->n_volumes + 2;

  FILE *file = safe_file_open(filename, "w");
  if (file == NULL) return;

  fprintf(file, "x (m),T (C)"); // compatible with csv files (otherwise wont open in colab)
  for (int i = 0; i < n_points; i++) {
    double x;
    if (i == 0) {
      x = 0.0;  // Borde izquierdo
    } else if (i == n_points - 1) {
      x = params->L;  // Borde derecho
    } else {
      x = (i - 0.5) * params->dx;  // Centro del volumen
    }

    fprintf(file, "%.6f,%.6f\n", x, T[i]);
  }

  safe_file_close(file);
}

void save_performance_metrics_csv(const PerformanceMetrics *metrics,
                                  const char *filename) {
  // TODO
}

void save_transient_profiles_csv(const SimulationParams *params,
                                 const char *filename) {
  if (params == NULL || filename == NULL) {
    fprintf(stderr, "Error: Invalid parameters\n");
    return;
  }
  char filename_temp[256];

  for (int i = 0; i < params->n_profiles; i++) {
    sprintf(filename_temp, "%s_profile_%d_at_%.2fs.csv", filename, i + 1,
            params->time_samples[i]);
    save_temperature_profile_csv(params->T_profiles[i], params,
                                 params->time_samples[i], filename_temp);
  }
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
  printf("\n--- TABLE: %s ---\n", label);
  printf("Position [m] | Temperature [C]\n");
  printf("-----------------------------\n");

  int n_points = n_volumes + 2;
  for (int i = 0; i < n_points; i++) {
    double x = i * 0.1 / (n_points - 1);  // Posición normalizada
    printf("%10.4f | %12.4f", x, T[i]);

    if (i == 0) {
      printf("  [Left Boundary]\n");
    } else if (i == n_points - 1) {
      printf("  [Right Boundary]\n");
    } else {
      printf("  [Inner Nodes]\n");
    }
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

  SimulationParams params;
  initialize_default_parameters(&params);

  printf("1. Testing convergence checker...\n");
  // validate convergence (fail)
  params.n_volumes = 5;
  params.dt = 9; // unstable
  calculate_derived_parameters(&params);
  int not_converge = !check_stability_condition(&params); // should not converge 
  printf("Test covergence failure: %s\n", not_converge ? "PASS" : "FAILED");
  // validate convergence (pass)
  params.dt = 1; // stable
  calculate_derived_parameters(&params);
  int converge = check_stability_condition(&params); // should converge
  printf("Test covergence pass: %s\n", converge ? "PASS" : "FAILED");

  printf("2. Testing parallel correctness...\n");
  // verify_parallel_correctness(&params);

  printf("=== TESTS COMPLETED ===\n");
}

void verify_parallel_correctness(const SimulationParams *params) {
  printf("\n=== PARALLEL CORRECTNESS VERIFICATION ===\n");
  SimulationParams params_seq = *params;
  SimulationParams params_par = *params;
  double *T_seq = allocate_temperature_field(params_seq.n_volumes);
  double *T_par = allocate_temperature_field(params_par.n_volumes);

  if (T_seq == NULL || T_par == NULL) {
    fprintf(stderr, "Error: Could not allocate memory for verification\n");
    return;
  }

  // Ejecutar ambas versiones
  solve_transient_sequential(T_seq, &params_seq);

  // Probar con diferentes números de hilos
  int num_threads = OMP_NUM_THREADS;
  configure_omp_environment(num_threads);
  solve_transient_parallel(T_par, &params_par);

  double tolerance = 1e-12;
  int equivalent = 1;
  for (int i = 0; i < params->n_profiles; i++) {
    equivalent = verify_solution_equivalence(params_seq.T_profiles[i],
                                             params_par.T_profiles[i],
                                             params->n_volumes, tolerance);
    if (!equivalent) break;
  }
  printf("Correctness: %s\n", equivalent ? "OK" : "FAILED");
  free_temperature_field(T_seq);
  free_temperature_field(T_par);
  free_temperature_profiles(params_seq.T_profiles, params_seq.n_profiles);
  free_temperature_profiles(params_par.T_profiles, params_par.n_profiles);

  printf("=== VERIFICATION COMPLETED ===\n");
}
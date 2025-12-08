// %%writefile main.c
// #include "functions.h"
#include "../include/functions.h"

void print_configuration(const SimulationParams *params) {
  printf("=== SIMULATION PARAMETERS ===\n");
  printf("Physical parameters:\n");
  printf("  Thermal conductivity (k): %.1f W/mK\n", params->k);
  printf("  rho*c: %.1e J/m³K\n", params->rho_c);
  printf("  Thermal diffusivity (alpha): %.2e m²/s\n", params->alpha);
  printf("\nGeometry:\n");
  printf("  Length (L): %.2e m\n", params->L);
  printf("  Number of volumes: %d\n", params->n_volumes);
  printf("  Spatial step (dx): %.2e m\n", params->dx);
  printf("\nTemperatures:\n");
  printf("  Initial temperature: %.1f °C\n", params->T_initial);
  printf("  Cooled surface temperature: %.1f °C\n", params->T_cooled);
  printf("\nTime parameters:\n");
  printf("  Total time: %.1f s\n", params->total_time);
  printf("  Time step (dt): %.3e s\n", params->dt);
  printf("  Number of time steps: %d\n", params->n_time_steps);
  printf("  Node Coefficients: %lf, %lf, %lf, %lf, %lf\n", params->aW, 
        params->aE, 
        params->aP, 
        params->aP0,
        params->aEb);
}

int main() {
  printf("=== HEAT EQUATION SOLVER - 1D EXPLICIT METHOD ===\n");
  printf("Using default parameters\n\n");

  // Inicializar parámetros con valores por defecto
  SimulationParams params;
  initialize_default_parameters(&params);

  // Imprimir configuración
  print_configuration(&params);

  // Validar parámetros
  if (!validate_parameters(&params)) {
    fprintf(stderr, "Error: Invalid parameters\n");
    return 1;
  }

  // Ejecutar tests de corrección
  printf("\n");
  run_correctness_test();

  // Comparar versiones secuencial y paralela
  printf("\n");
  PerformanceMetrics metrics = compare_sequential_vs_parallel(&params);

  // Análisis de escalabilidad
  printf("\n");
  benchmark_scalability_analysis(&params);

  // Ejemplo de simulación transitoria
  printf("\n=== TRANSIENT SIMULATION EXAMPLE ===\n");
  double *T = allocate_temperature_field(params.n_volumes);
  double *T_profiles =
      malloc(params.n_volumes * 5 * sizeof(double));  // 5 perfiles

  if (T == NULL || T_profiles == NULL) {
    fprintf(stderr, "Error: Could not allocate memory\n");
    return 1;
  }

  // Definir tiempos para guardar perfiles (20%, 40%, 60%, 80%, 100%)
  int profile_indices[5];
  for (int i = 0; i < 5; i++) {
    profile_indices[i] = (i + 1) * params.n_time_steps / 5;
  }

  // Ejecutar simulación transitoria secuencial
  solve_transient_sequential(T, &params, T_profiles, profile_indices);

  // Guardar último perfil
  printf("\nSaving results...\n");
  save_temperature_profile_csv(T, &params, params.total_time,
                               "final_profile.csv");

  // Guardar métricas
  save_performance_metrics_csv(&metrics, "performance_metrics.csv");

  // Calcular y mostrar temperaturas finales
  double max_temp, min_temp;
  find_temperature_extremes(T, params.n_volumes, &max_temp, &min_temp);
  printf("\nFinal temperature distribution:\n");
  printf("  Maximum temperature: %.2f °C\n", max_temp);
  printf("  Minimum temperature: %.2f °C\n", min_temp);
  printf("  Temperature drop: %.2f °C\n", params.T_initial - min_temp);

  // Calcular energía total
  double final_energy = calculate_total_energy(T, &params);
  printf("  Final thermal energy: %.6e J\n", final_energy);

  // Liberar memoria
  free_temperature_field(T);
  free(T_profiles);
  return 0;
}
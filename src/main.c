#include "../include/functions.h"

void print_parameters(SimulationParams *params) {
  printf("\n---Parameters---\n");

  printf("\nPhysical Parameters:\n");
  printf("    k: %.3f W/mK\n", params->k);
  printf("rho_c: %.3f J/m^3K\n", params->rho_c);
  printf("alpha: %.3f m^2/s\n", params->alpha);

  printf("\nGeometry Parameters:\n");
  printf("        L: %.3f m\n", params->L);
  printf("       dx: %.6f m\n", params->dx);
  printf("n_volumes: %d\n", params->n_volumes);

  printf("\nThermal Conditions:\n");
  printf("T_initial: %.1f C\n", params->T_initial);
  printf(" T_cooled: %.1f C\n", params->T_cooled);

  printf("\nTemporal Parameters:\n");
  printf("          dt: %.3f s\n", params->dt);
  printf("  total_time: %.1f s\n", params->total_time);
  printf("n_time_steps: %d\n", params->n_time_steps);

  printf("\nNumerical Coefficients:\n");
  printf(" aW: %.3f\n", params->aW);
  printf(" aE: %.3f\n", params->aE);
  printf(" aP: %.3f\n", params->aP);
  printf("aP0: %.3f\n", params->aP0);
  printf("aEb: %.3f\n", params->aEb);

  printf("\nProfile Parameters:\n");
  printf("   n_profiles: %d\n", params->n_profiles);
  printf("profile_times:\n");
  for (int i = 0; i < params->n_profiles; i++) {
    printf("            %d: %.1f s\n", i + 1, params->time_samples[i]);
  }
  printf("---End Parameters---\n\n");
}

int main() {
  SimulationParams params;
  initialize_default_parameters(&params);
  print_parameters(&params);

  run_correctness_test();

  double *T = allocate_temperature_field(params.n_volumes);

  // 1. Resolver ecuación de calor
  solve_heat_equation_sequential(T, &params);
  save_temperature_profile_csv(T, &params, params.total_time,
                               "../data/sequential_solve_heat_equation.csv");

  print_temperature_field(T, params.n_volumes,
                          "../data/sequential_solve_heat_equation.csv");

  // 2. Resolver transitoria
  validate_convergence(&params);
  solve_transient_sequential(T, &params);
  save_transient_profiles_csv(&params, "../data/transient_seq/sequential_solve_transient");
  // liberar T_profiles si fue asignado
  if (params.T_profiles != NULL) {
    free_temperature_profiles(params.T_profiles, params.n_profiles);
  }

  solve_transient_parallel(&params, 2);
  save_transient_profiles_csv(&params, "../data/transient_par/parallel_solve_transient");
  // liberar T_profiles si fue asignado
  if (params.T_profiles != NULL) {
    free_temperature_profiles(params.T_profiles, params.n_profiles);
  }

  // Liberar memoria
  free_temperature_field(T);


  return 0;

  // Liberar memoria
  free_temperature_field(T);

  return 0;
}
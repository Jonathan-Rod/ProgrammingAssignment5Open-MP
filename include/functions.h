// %%writefile functions.h
#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// ============================================================================
// CONSTANTES Y CONFIGURACIÓN
// ============================================================================

#define MAX_VOLUMES 1000
#define MAX_FILENAME 256
#define CSV_DELIMITER ","
#define MAX_TIME_PROFILES 10000
#define STABILITY_LIMIT 0.5  // Límite de Fourier number para esquema explícito

// Configuración OpenMP por defecto
#ifndef OMP_NUM_THREADS
#define OMP_NUM_THREADS 4
#endif

// ============================================================================
//  ESTRUCTURAS PRINCIPALES
// ============================================================================

/**
 * @brief Estructura que contiene todos los parámetros de simulación
 *
 * Esta estructura almacena los parámetros físicos, geométricos y numéricos
 * necesarios para resolver la ecuación de calor transitoria 1D con método
 * explícito.
 */
typedef struct {
  // Parámetros físicos
  double k;      // Conductividad térmica [W/mK]
  double rho_c;  // Densidad * calor específico [J/m³K]
  double alpha;  // Difusividad térmica [m²/s] (k/rho_c)

  // Geometría y dominio
  double L;       // Longitud del dominio [m]
  double dx;      // Espaciado espacial [m]
  int n_volumes;  // Número de volúmenes de control

  // Condiciones iniciales y de frontera
  double T_initial;  // Temperatura inicial [°C]
  double T_cooled;   // Temperatura superficie enfriada [°C]

  // Parámetros temporales
  double dt;          // Paso de tiempo [s]
  double total_time;  // Tiempo total de simulación [s]
  int n_time_steps;   // Número total de pasos de tiempo

  // Node Coefficients
  double aW;   // West point's coefficient
  double aE;   // East point's coefficient
  double aP;   // Point's coefficient
  double aP0;  // Point's coefficient (Previous time step)
  double aEb;  // East point's coefficient (at Boundary)

  // Para análisis de perfiles transitorios
  int n_profiles;                          // Número de perfiles a guardar
  double** T_profiles;                     // Array de Arrays
  double time_samples[MAX_TIME_PROFILES];  // Tiempos para guardar perfiles
} SimulationParams;                        // DONE

/**
 * @brief Estructura para almacenar métricas de performance computacional
 *
 * Contiene tiempos de ejecución, speedup, eficiencia y otras métricas
 * para comparar implementaciones secuenciales y paralelas.
 */
typedef struct {
  double sequential_time;  // Tiempo de ejecución secuencial [s]
  double parallel_time;    // Tiempo de ejecución paralelo [s]
  double speedup;          // Speedup = sequential_time/parallel_time
  double efficiency;       // Eficiencia = speedup/n_threads
  int optimal_threads;     // Número óptimo de hilos encontrado
  double max_temperature;  // Temperatura máxima en el dominio
  double min_temperature;  // Temperatura mínima en el dominio
} PerformanceMetrics;      // DONE

// ============================================================================
// FUNCIONES DE INICIALIZACIÓN Y CONFIGURACIÓN
// ============================================================================

/**
 * @brief Inicializa los parámetros de simulación con valores por defecto
 *
 * Paso a paso:
 * 1. Establece parámetros físicos por defecto (k, rho_c)
 * 2. Configura geometría básica (L, n_volumes)
 * 3. Define condiciones iniciales y de frontera
 * 4. Establece parámetros temporales conservadores
 * 5. Inicializa arreglos auxiliares a cero
 *
 * @param params Puntero a estructura de parámetros a inicializar
 */
void initialize_default_parameters(SimulationParams* params);  // DONE

/**
 * @brief Calcula parámetros derivados basados en los parámetros básicos
 *
 * Paso a paso:
 * 1. Calcula espaciado espacial dx = L/(n_volumes-1)
 * 2. Calcula difusividad térmica alpha = k/rho_c
 * 3. Determina número de pasos de tiempo n_time_steps = total_time/dt
 * 4. Verifica consistencia dimensional
 *
 * @param params Puntero a estructura de parámetros (calcula dx, alpha,
 * n_time_steps)
 */
void calculate_derived_parameters(SimulationParams* params);  // DONE

/**
 * @brief Asigna memoria para el campo de temperaturas
 *
 * Paso a paso:
 * 1. Verifica que n_volumes esté dentro de límites razonables
 * 2. Asigna memoria para n_volumes elementos double
 * 3. Inicializa todos los elementos a cero
 * 4. Retorna puntero al arreglo asignado
 *
 * @param n_volumes Número de volúmenes de control
 * @return Puntero al arreglo de temperaturas asignado
 */
double* allocate_temperature_field(int n_volumes);  // DONE

/**
 * @brief
 *
 *
 * @param n_profiles
 * @param n_volumes
 * @return double**
 */
double** allocate_temperature_profiles(int n_profiles, int n_volumes);  // DONE

/**
 * @brief Libera memoria del campo de temperaturas
 *
 * Paso a paso:
 * 1. Verifica que el puntero no sea NULL
 * 2. Libera la memoria asignada
 * 3. Establece el puntero a NULL (buena práctica)
 *
 * @param T Puntero al arreglo de temperaturas a liberar
 */
void free_temperature_field(double* T);  // DONE

/**
 * @brief
 *
 * @param TT
 * @param profiles
 */
void free_temperature_profiles(double** TT, int profiles);

/**
 * @brief Valida que los parámetros de simulación sean físicamente
 * consistentes
 *
 * Paso a paso:
 * 1. Verifica parámetros físicos positivos (k, rho_c, L)
 * 2. Valida número de volúmenes dentro de límites
 * 3. Chequea condiciones de estabilidad numérica para esquema explícito
 * 4. Verifica consistencia dimensional
 * 5. Retorna código de error específico si falla
 *
 * @param params Puntero a estructura de parámetros a validar
 * @return 1 si parámetros son válidos, 0 en caso contrario
 */
int validate_parameters(const SimulationParams* params);  // DONE

// ============================================================================
// CÁLCULOS DE ESTABILIDAD Y COEFICIENTES (EXPLÍCITO)
// ============================================================================

/**
 * @brief Calcula el límite de estabilidad para el esquema explícito
 *
 * Paso a paso:
 * 1. Calcula número de Fourier Fo = alpha*dt/dx²
 * 2. Determina límite de estabilidad para esquema explícito (Fo <= 0.5)
 * 3. Considera condiciones de frontera específicas
 * 4. Aplica factor de seguridad conservador
 *
 * @param params Puntero a parámetros de simulación
 * @return Límite máximo de paso de tiempo para estabilidad
 */
double calculate_stability_limit(const SimulationParams* params);  // DONE

/**
 * @brief Verifica si el paso de tiempo cumple condición de estabilidad
 * explícita
 *
 * Paso a paso:
 * 1. Calcula número de Fourier actual Fo = alpha*dt/dx²
 * 2. Compara con límite de estabilidad teórico (Fo <= 0.5)
 * 3. Considera esquema explícito forward Euler
 * 4. Retorna resultado de la verificación
 *
 * @param params Puntero a parámetros de simulación
 * @return 1 si es estable, 0 si no cumple condición
 */
int check_stability_condition(
    const SimulationParams*
        params);  // DONE depends on calculate_fourier_number

/**
 * @brief Calcula la difusividad térmica basada en k y rho_c
 *
 * Paso a paso:
 * 1. Verifica que rho_c no sea cero
 * 2. Calcula alpha = k / rho_c
 * 3. Almacena resultado en estructura params
 * 4. Verifica resultado físicamente razonable
 *
 * @param params Puntero a estructura donde se almacenará alpha
 */
void calculate_thermal_diffusivity(SimulationParams* params);  // DONE

/**
 * @brief Calcula el número de Fourier para el esquema explícito
 *
 * Paso a paso:
 * 1. Calcula Fo = alpha * dt / (dx * dx)
 * 2. Verifica que esté dentro de límites de estabilidad
 * 3. Retorna el número de Fourier calculado
 *
 * @param params Puntero a parámetros de simulación
 * @return Número de Fourier calculado
 */
double calculate_fourier_number(const SimulationParams* params);  // DONE

// ============================================================================
// SIMULACIÓN SECUENCIAL (EXPLÍCITA)
// ============================================================================

/**
 * @brief Resuelve la ecuación de calor de forma secuencial con método explícito
 *
 * Paso a paso:
 * 1. Inicializa campo de temperaturas con condición inicial
 * 2. Aplica condiciones de frontera
 * 3. Para cada paso de tiempo:
 *    a. Calcula nuevas temperaturas usando esquema explícito
 *    b. Aplica condiciones de frontera
 *    c. Actualiza campo de temperaturas
 * 4. Calcula métricas de convergencia
 *
 * @param T Arreglo de temperaturas (entrada/salida)
 * @param params Parámetros de simulación
 */
void solve_heat_equation_sequential(double* T,
                                    const SimulationParams* params);  // DONE

/**
 * @brief Simulación transitoria secuencial explícita guardando perfiles
 * específicos
 *
 * Paso a paso:
 * 1. Identifica tiempos donde guardar perfiles
 * 2. Para cada paso de tiempo:
 *    a. Avanza solución un paso con método explícito
 *    b. Si es tiempo de guardar perfil, copia T a T_profiles
 *    c. Actualiza contadores de progreso
 * 3. Guarda perfiles en estructura organizada
 *
 * @param T Arreglo de temperaturas inicial
 * @param params Parámetros de simulación
 */
void solve_transient_sequential(double* T,
                                const SimulationParams* params);  // TODO

/**
 * @brief Integración temporal secuencial explícita del campo de temperaturas
 *
 * Paso a paso:
 * 1. Asigna memoria para campo temporal
 * 2. Para cada paso de tiempo:
 *    a. Calcula nuevas temperaturas con esquema explícito
 *    b. Aplica condiciones de frontera
 *    c. Intercambia punteros para siguiente iteración
 * 3. Libera memoria temporal
 *
 * @param T Arreglo de temperaturas (se actualiza)
 * @param params Parámetros de simulación
 */
void time_integration_sequential(double* T,
                                 const SimulationParams* params);  // TODO

/**
 * @brief Calcula nuevo paso de tiempo usando esquema explícito
 *
 * Paso a paso:
 * 1. Para cada volumen interno:
 *    a. Calcula termino independiente de la nueva temperatura (b)
 *    b. Aplica esquema backward Euler: T_new = b/T_old
 *    c. Almacena en arreglo temporal
 * 2. Maneja volúmenes frontera separadamente
 *
 * @param T_new Temperaturas del nuevo paso de tiempo (salida)
 * @param T_old Temperaturas del paso anterior
 * @param params Parámetros de simulación
 */
void calculate_explicit_step_sequential(
    double* T_new, const double* T_old,
    const SimulationParams* params);  // DONE

/**
 * @brief Aplica condiciones de frontera en implementación secuencial explícita
 *
 * Paso a paso:
 * 1. Aplica condición Dirichlet en frontera izquierda (T = T_cooled)
 * 2. Aplica condición Neumann en frontera derecha (dT/dx = 0, aislada)
 * 3. Actualiza temperaturas de frontera directamente
 *
 * @param T Campo de temperaturas (modificado)
 * @param params Parámetros de simulación
 */
void apply_boundary_conditions_sequential(
    double* T,
    const SimulationParams* params);  // DONE

// ============================================================================
// SIMULACIÓN PARALELA (OPENMP - EXPLÍCITA)
// ============================================================================

/**
 * @brief Resuelve la ecuación de calor usando paralelización OpenMP con método
 * explícito
 *
 * Paso a paso:
 * 1. Configura entorno OpenMP
 * 2. Particiona dominio entre hilos
 * 3. Paraleliza cálculo de paso explícito
 * 4. Sincroniza hilos en actualización de temperatura
 * 5. Aplica condiciones de frontera
 *
 * @param T Arreglo de temperaturas (entrada/salida)
 * @param params Parámetros de simulación
 */
void solve_heat_equation_parallel(double* T,
                                  const SimulationParams* params);  // TODO

/**
 * @brief Simulación transitoria paralela explícita guardando perfiles
 * específicos
 *
 *
 * @param T Arreglo de temperaturas inicial
 * @param params Parámetros de simulación
 */
void solve_transient_parallel(double* T,
                              const SimulationParams* params);  // TODO

/**
 * @brief Integración temporal paralela explícita del campo de temperaturas
 *
 * Paso a paso:
 * 1. Paraleliza cálculo de paso explícito con OpenMP
 * 2. Usa directivas for para paralelizar volúmenes
 * 3. Aplica condiciones frontera después de sincronización
 * 4. Minimiza contención en acceso a memoria
 *
 * @param T Arreglo de temperaturas (se actualiza)
 * @param params Parámetros de simulación
 */
void time_integration_parallel(double* T,
                               const SimulationParams* params);  // TODO

/**
 * @brief Calcula nuevo paso de tiempo usando esquema explícito en paralelo
 *
 * Paso a paso:
 * 1. Divide volúmenes internos entre hilos
 * 2. Cada hilo calcula sus volúmenes asignados
 * 3. Usa schedule estático para mejor localidad de datos
 * 4. Combina resultados automáticamente
 *
 * @param T_new Temperaturas del nuevo paso de tiempo (salida)
 * @param T_old Temperaturas del paso anterior
 * @param params Parámetros de simulación
 */
void calculate_explicit_step_parallel(double* T_new, const double* T_old,
                                      const SimulationParams* params);  // TODO

/**
 * @brief Aplica condiciones de frontera en implementación paralela explícita
 *
 * Paso a paso:
 * 1. Aplica condiciones frontera izquierda en hilo maestro
 * 2. Aplica condiciones frontera derecha en último hilo
 * 3. Sincroniza hilos después de aplicar condiciones
 * 4. Usa directivas single para condiciones globales
 *
 * @param T Campo de temperaturas (modificado)
 * @param params Parámetros de simulación
 */
void apply_boundary_conditions_parallel(
    double* T,
    const SimulationParams* params);  // TODO

// ============================================================================
// ANÁLISIS DE PERFORMANCE Y BENCHMARKING
// ============================================================================

/**
 * @brief Compara performance entre implementación secuencial y paralela
 *
 * Paso a paso:
 * 1. Ejecuta simulación secuencial y mide tiempo
 * 2. Ejecuta simulación paralela y mide tiempo
 * 3. Calcula speedup y eficiencia
 * 4. Compara resultados numéricos
 * 5. Genera reporte de performance
 *
 * @param params Parámetros de simulación
 * @return Estructura con métricas de performance comparativas
 */
PerformanceMetrics compare_sequential_vs_parallel(
    const SimulationParams* params);  // TODO

/**
 * @brief Barrido de parámetros para análisis de sensibilidad de performance
 *
 * Paso a paso:
 * 1. Varía parámetros clave (MAX_TIME_PROFILES)
 * 2. Para cada combinación, mide performance
 * 3. Identifica cuellos de botella
 * 4. Genera superficies de respuesta
 *
 * @param base_params Parámetros base para la simulación
 */
void performance_sweep_parameters(const SimulationParams* base_params);  // TODO

/**
 * @brief Calcula ratio de speedup entre versiones
 *
 * Paso a paso:
 * 1. Verifica tiempos válidos (positivos, no cero)
 * 2. Calcula ratio seq_time/par_time
 * 3. Aplica límites razonables (máximo teórico)
 * 4. Retorna speedup calculado
 *
 * @param seq_time Tiempo de ejecución secuencial
 * @param par_time Tiempo de ejecución paralelo
 * @return Ratio de speedup (seq_time/par_time)
 */
double calculate_speedup_ratio(double seq_time, double par_time);  // DONE

/**
 * @brief Calcula eficiencia paralela
 *
 * Paso a paso:
 * 1. Verifica inputs válidos
 * 2. Calcula eficiencia = speedup/n_threads
 * 3. Aplica normalización (0-100%)
 * 4. Retorna eficiencia calculada
 *
 * @param speedup Ratio de speedup obtenido
 * @param n_threads Número de hilos utilizados
 * @return Eficiencia paralela (speedup/n_threads)
 */
double calculate_parallel_efficiency(double speedup, int n_threads);  // DONE

/**
 * @brief Mide tiempo de ejecución de una función de solver
 *
 * Paso a paso:
 * 1. Tiempo inicial con omp_get_wtime()
 * 2. Ejecuta función solver proporcionada
 * 3. Tiempo final con omp_get_wtime()
 * 4. Calcula diferencia y retorna
 *
 * @param solver Función solver a medir
 * @param T Arreglo de temperaturas
 * @param params Parámetros de simulación
 * @return Tiempo de ejecución en segundos
 */
double measure_execution_time(void (*solver)(double*, const SimulationParams*),
                              double* T,
                              const SimulationParams* params);  // DONE

// ============================================================================
// VALIDACIÓN Y VERIFICACIÓN (EXPLÍCITA)
// ============================================================================

/**
 * @brief Verifica equivalencia entre soluciones secuencial y paralela
 *
 * Paso a paso:
 * 1. Calcula diferencia punto a punto
 * 2. Encuentra máxima diferencia
 * 3. Compara con tolerancia especificada
 * 4. Verifica patrones de error sistemáticos
 *
 * @param T_seq Solución secuencial
 * @param T_par Solución paralela
 * @param n_volumes Número de volúmenes
 * @param tolerance Tolerancia para comparación
 * @return 1 si son equivalentes, 0 en caso contrario
 */
int verify_solution_equivalence(const double* T_seq, const double* T_par,
                                int n_volumes, double tolerance);  // DONE

/**
 * @brief Calcula error numérico comparando con solución analítica (si
 * disponible)
 *
 * Paso a paso:
 * 1. Calcula solución analítica para caso simple
 * 2. Calcula normas de error (L1, L2, Linf)
 * 3. Considera comportamiento asintótico
 * 4. Retorna métrica de error principal
 *
 * @param T_numeric Solución numérica
 * @param params Parámetros de simulación
 * @param current_time Tiempo actual de simulación
 * @return Magnitud del error numérico
 */
double calculate_numerical_error(const double* T_numeric,
                                 const SimulationParams* params,
                                 double current_time);  // !!!

/**
 * @brief Valida historia de convergencia de la simulación explícita
 *
 * Paso a paso:
 * 1. Analiza estabilidad del esquema explícito
 * 2. Verifica que temperaturas permanecen dentro de límites físicos
 * 3. Identifica inestabilidades numéricas
 * 4. Genera diagnóstico de convergencia
 *
 * @param params Parámetros de simulación
 */
void validate_convergence_history(const SimulationParams* params);  // TODO

/**
 * @brief Encuentra temperaturas máxima y mínima en el dominio
 *
 * Paso a paso:
 * 1. Inicializa max y min con primer valor
 * 2. Recorre todos los volúmenes
 * 3. Actualiza max y min
 * 4. Retorna en estructura de métricas
 *
 * @param T Campo de temperaturas
 * @param n_volumes Número de volúmenes
 * @param max_temp Temperatura máxima (salida)
 * @param min_temp Temperatura mínima (salida)
 */
void find_temperature_extremes(const double* T, int n_volumes, double* max_temp,
                               double* min_temp);  // TODO

// ============================================================================
// GESTIÓN DE DATOS Y ARCHIVOS
// ============================================================================

/**
 * @brief Guarda perfil de temperatura en archivo CSV
 *
 * Paso a paso:
 * 1. Abre archivo con manejo de errores
 * 2. Escribe headers descriptivos
 * 3. Para cada volumen, escribe posición y temperatura
 * 4. Incluye metadatos de simulación
 * 5. Cierra archivo correctamente
 *
 * @param T Arreglo de temperaturas
 * @param params Parámetros de simulación
 * @param current_time Tiempo actual de simulación
 * @param filename Nombre del archivo CSV
 */
void save_temperature_profile_csv(const double* T,
                                  const SimulationParams* params,
                                  double current_time,
                                  const char* filename);  // TODO

/**
 * @brief Guarda métricas de performance en archivo CSV
 *
 * Paso a paso:
 * 1. Crea/abre archivo de métricas
 * 2. Escribe headers de columnas
 * 3. Formatea números con precisión adecuada
 * 4. Incluye información de configuración
 * 5. Cierra archivo de forma segura
 *
 * @param metrics Métricas a guardar
 * @param filename Nombre del archivo CSV
 */
void save_performance_metrics_csv(const PerformanceMetrics* metrics,
                                  const char* filename);  // TODO

/**
 * @brief Guarda múltiples perfiles transitorios en CSV
 *
 * Paso a paso:
 * 1. Organiza datos en formato matriz
 * 2. Escribe tiempos como primera columna
 * 3. Para cada tiempo, escribe perfil completo
 * 4. Usa formato compatible con herramientas de visualización
 *
 * @param T_profiles Arreglo 2D de perfiles de temperatura
 * @param params Parámetros de simulación
 * @param filename Nombre del archivo CSV
 */
void save_transient_profiles_csv(const double* T_profiles,
                                 const SimulationParams* params,
                                 const char* filename);  // TODO

/**
 * @brief Guarda datos de escalabilidad en CSV
 *
 * Paso a paso:
 * 1. Para cada configuración de hilos
 * 2. Escribe número de hilos, tiempo, speedup, eficiencia
 * 3. Incluye metadatos del sistema
 * 4. Formato para fácil análisis posterior
 *
 * @param metrics_array Arreglo de métricas para diferentes configuraciones
 * @param num_configs Número de configuraciones
 * @param filename Nombre del archivo CSV
 */
void save_scalability_data_csv(const PerformanceMetrics* metrics_array,
                               int num_configs, const char* filename);  // TODO

/**
 * @brief Abre archivo de forma segura con verificación de errores
 *
 * Paso a paso:
 * 1. Intenta abrir archivo en modo especificado
 * 2. Verifica que apertura fue exitosa
 * 3. Maneja errores específicos (permisos, espacio, etc.)
 * 4. Retorna NULL en caso de error con mensaje
 *
 * @param filename Nombre del archivo
 * @param mode Modo de apertura
 * @return Puntero FILE* o NULL en caso de error
 */
FILE* safe_file_open(const char* filename, const char* mode);  // TODO

/**
 * @brief Cierra archivo de forma segura
 *
 * Paso a paso:
 * 1. Verifica que puntero no sea NULL
 * 2. Cierra archivo
 * 3. Verifica que cierre fue exitoso
 * 4. Maneja errores de cierre
 *
 * @param file Puntero al archivo a cerrar
 */
void safe_file_close(FILE* file);  // TODO

/**
 * @brief Escribe headers en archivo CSV
 *
 * Paso a paso:
 * 1. Escribe primer header sin coma
 * 2. Para headers subsiguientes, añade coma
 * 3. Termina línea con newline
 * 4. Flushea buffer para persistencia
 *
 * @param file Puntero al archivo
 * @param headers Arreglo de strings con headers
 * @param n_headers Número de headers
 */
void write_csv_headers(FILE* file, const char* headers[],
                       int n_headers);  // TODO

// ============================================================================
// CONFIGURACIÓN Y CONTROL OPENMP
// ============================================================================

/**
 * @brief Configura entorno OpenMP con número específico de hilos
 *
 * Paso a paso:
 * 1. Establece número de hilos con omp_set_num_threads()
 * 2. Configura variables de entorno si es necesario
 * 3. Verifica que configuración fue aplicada
 * 4. Reporta configuración resultante
 *
 * @param num_threads Número de hilos a utilizar
 */
void configure_omp_environment(int num_threads);

/**
 * @brief Imprime información de configuración OpenMP
 *
 * Paso a paso:
 * 1. Obtiene número máximo de hilos disponibles
 * 2. Obtiene número actual de hilos
 * 3. Imprime información del procesador
 * 4. Muestra configuración de scheduling
 */
void print_omp_configuration_info(void);

/**
 * @brief Obtiene número de hilos paralelos disponibles
 *
 * Paso a paso:
 * 1. Consulta número de procesadores lógicos
 * 2. Considera límites del sistema
 * 3. Ajusta por disponibilidad actual
 * 4. Retorna número usable de hilos
 *
 * @return Número de hilos disponibles
 */
int get_available_parallel_threads(void);

/**
 * @brief Configura scheduling dinámico de OpenMP
 *
 * Paso a paso:
 * 1. Establece política de scheduling
 * 2. Configura chunk size óptimo
 * 3. Ajusta según características del problema
 * 4. Aplica a regiones paralelas
 *
 * @param chunk_size Tamaño de chunk para scheduling
 */
void set_omp_dynamic_scheduling(int chunk_size);

// ============================================================================
// UTILIDADES DE VISUALIZACIÓN Y DEBUG
// ============================================================================

/**
 * @brief Imprime campo de temperaturas con etiqueta
 *
 * Paso a paso:
 * 1. Imprime etiqueta descriptiva
 * 2. Para cada volumen, imprime posición y temperatura
 * 3. Formatea números para legibilidad
 * 4. Incluye separadores visuales
 *
 * @param T Arreglo de temperaturas
 * @param n_volumes Número de volúmenes
 * @param label Etiqueta descriptiva
 */
void print_temperature_field(const double* T, int n_volumes, const char* label);

/**
 * @brief Imprime progreso de la simulación
 *
 * Paso a paso:
 * 1. Calcula porcentaje completado
 * 2. Formatea barra de progreso ASCII
 * 3. Muestra temperatura máxima y mínima actual
 * 4. Actualiza en misma línea (carriage return)
 *
 * @param iteration Iteración actual
 * @param total Iteraciones totales
 * @param max_temp Temperatura máxima actual
 * @param min_temp Temperatura mínima actual
 */
void print_simulation_progress(int iteration, int total, double max_temp,
                               double min_temp);

/**
 * @brief Visualiza particionamiento de dominio para debug
 *
 * Paso a paso:
 * 1. Calcula rango de cada hilo
 * 2. Dibuja representación visual del dominio
 * 3. Muestra volúmenes asignados a cada hilo
 * 4. Identifica desbalance de carga
 *
 * @param n_volumes Número de volúmenes
 * @param n_threads Número de hilos
 */
void visualize_domain_partitioning(int n_volumes, int n_threads);

/**
 * @brief Imprime resumen de métricas de performance
 *
 * Paso a paso:
 * 1. Formatea tiempos con unidades apropiadas
 * 2. Calcula porcentajes de mejora
 * 3. Destaca métricas clave
 * 4. Proporciona interpretación cualitativa
 *
 * @param metrics Métricas a imprimir
 */
void print_performance_summary(const PerformanceMetrics* metrics);

// ============================================================================
// FUNCIONES DE PRUEBA Y VERIFICACIÓN
// ============================================================================

/**
 * @brief Ejecuta tests de corrección de la implementación explícita
 * 1. test boundary conditions
 * 2. calculate numerical error
 */
void run_correctness_test(void);

/**
 * @brief Verifica correcta aplicación de condiciones de frontera
 *
 * Paso a paso:
 * 1. Configura condiciones de frontera específicas
 * 2. Ejecuta un paso de simulación
 * 3. Verifica que condiciones se mantienen
 * 4. Chequea consistencia con solución interna
 */
void test_boundary_conditions(void);

/**
 * @brief Verifica corrección de implementación paralela explícita
 *
 * Paso a paso:
 * 1. Compara con solución secuencial de referencia
 * 2. Verifica que no hay race conditions
 * 3. Chequea consistencia en interfaces entre hilos
 * 4. Valida para diferentes números de hilos
 *
 * @param params Parámetros de simulación
 */
void verify_parallel_correctness(const SimulationParams* params);

#endif  // FUNCTIONS_H
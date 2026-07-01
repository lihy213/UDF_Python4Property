#include "udf.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*
 * Fast uniform-grid version of CO2_propertiesV3_Python4Property.c.
 *
 * Use this file only with strictly uniform CSV tables:
 *   P(MPa)/T(K), 220.0, 220.5, ...
 *   0.05, value, value, ...
 *
 * Compared with the nonuniform-grid UDF:
 *   - Table reading is still flexible and reads CSV dimensions automatically.
 *   - Runtime interpolation uses direct row/column indexing instead of binary search.
 *   - Phase diagnostics are controlled by ENABLE_PHASE_DIAGNOSTICS.
 *   - Cp enthalpy uses h = cp * (T - Tref), following the common Fluent UDF example.
 */

#define CSV_MAX_LINE_CHARS 1048576
#define INITIAL_ROW_CAPACITY 512
#define CP_REFERENCE_PRESSURE_PA 14000000.0
#define CO2_MISSING_VALUE 0.0

/* Set to 1 during debugging, 0 for long production runs. */
#define ENABLE_PHASE_DIAGNOSTICS 0
#define ENABLE_RANGE_WARNINGS 0
#define WARNING_LIMIT 20

/* Relative tolerance used to judge whether the input table is uniform. */
#define UNIFORM_REL_TOL 1.0e-8

typedef struct {
    int nt;
    int np;
    double *t_grid;
    double *p_grid;
    double *values;
    double t_min;
    double t_max;
    double p_min;
    double p_max;
    double t_step;
    double p_step;
    int loaded;
    int uniform_ready;
} Table2D;

typedef struct {
    int nt;
    int np;
    int *phase;
    int loaded;
} PhaseTable2D;

enum {
    PHASE_UNKNOWN = 0,
    PHASE_GAS = 1,
    PHASE_LIQUID = 2,
    PHASE_SUPERCRITICAL = 3,
    PHASE_TWO_PHASE = 4,
    PHASE_INVALID = 5
};

static Table2D co2_density = {0};
static Table2D co2_viscosity = {0};
static Table2D co2_conductivity = {0};
static Table2D co2_cp_2d = {0};
static PhaseTable2D co2_phase = {0};

static double *cp_ref = NULL;
static int cp_ref_nt = 0;
static double cp_t_min = 0.0;
static double cp_t_max = 0.0;
static double cp_t_step = 0.0;

static double operating_pressure = 0.0;
static int range_warning_count = 0;
static int phase_warning_count = 0;

static void free_table(Table2D *table)
{
    if (table->t_grid) free(table->t_grid);
    if (table->p_grid) free(table->p_grid);
    if (table->values) free(table->values);
    table->t_grid = NULL;
    table->p_grid = NULL;
    table->values = NULL;
    table->nt = 0;
    table->np = 0;
    table->loaded = 0;
    table->uniform_ready = 0;
}

static void free_phase_table(PhaseTable2D *table)
{
    if (table->phase) free(table->phase);
    table->phase = NULL;
    table->nt = 0;
    table->np = 0;
    table->loaded = 0;
}

static int count_csv_tokens(char *line)
{
    int count = 0;
    char *token = strtok(line, ",\n\r");
    while (token != NULL) {
        count++;
        token = strtok(NULL, ",\n\r");
    }
    return count;
}

static int parse_phase_family(const char *text)
{
    if (text == NULL) return PHASE_UNKNOWN;
    if (strstr(text, "invalid") != NULL) return PHASE_INVALID;
    if (strstr(text, "two") != NULL || strstr(text, "phase_twophase") != NULL) return PHASE_TWO_PHASE;
    if (strstr(text, "supercritical") != NULL) {
        if (strstr(text, "gas") != NULL) return PHASE_GAS;
        if (strstr(text, "liquid") != NULL) return PHASE_LIQUID;
        return PHASE_SUPERCRITICAL;
    }
    if (strstr(text, "gas") != NULL) return PHASE_GAS;
    if (strstr(text, "liquid") != NULL) return PHASE_LIQUID;
    if (strstr(text, "critical") != NULL) return PHASE_SUPERCRITICAL;
    return PHASE_UNKNOWN;
}

static int grid_is_uniform(const double *grid, int n, double *step)
{
    int i;
    double expected_step;
    double diff;
    double tol;

    if (grid == NULL || n < 2) return 0;
    expected_step = grid[1] - grid[0];
    if (expected_step <= 0.0) return 0;

    tol = fabs(expected_step) * UNIFORM_REL_TOL;
    if (tol < 1.0e-12) tol = 1.0e-12;

    for (i = 1; i < n - 1; i++) {
        diff = fabs((grid[i + 1] - grid[i]) - expected_step);
        if (diff > tol) return 0;
    }

    *step = expected_step;
    return 1;
}

static int read_table_2d(const char *filename, Table2D *table)
{
    FILE *fp;
    char *line;
    char *header_copy;
    char *token;
    int nt, row_capacity, row_count, j;
    double *new_values;

    fp = fopen(filename, "r");
    if (fp == NULL) {
        Message("\n[CO2-UDF Error] Cannot open %s\n", filename);
        return 0;
    }

    line = (char*)malloc(CSV_MAX_LINE_CHARS);
    if (line == NULL) {
        Message("\n[CO2-UDF Error] Failed to allocate CSV line buffer.\n");
        fclose(fp);
        return 0;
    }

    if (fgets(line, CSV_MAX_LINE_CHARS, fp) == NULL) {
        Message("\n[CO2-UDF Error] Empty CSV file: %s\n", filename);
        free(line);
        fclose(fp);
        return 0;
    }

    header_copy = (char*)malloc(strlen(line) + 1);
    if (header_copy == NULL) {
        free(line);
        fclose(fp);
        return 0;
    }
    strcpy(header_copy, line);
    nt = count_csv_tokens(header_copy) - 1;
    free(header_copy);

    if (nt < 2) {
        Message("\n[CO2-UDF Error] Invalid temperature header in %s\n", filename);
        free(line);
        fclose(fp);
        return 0;
    }

    free_table(table);
    table->nt = nt;
    table->t_grid = (double*)malloc(sizeof(double) * nt);

    token = strtok(line, ",\n\r");
    for (j = 0; j < nt; j++) {
        token = strtok(NULL, ",\n\r");
        if (token == NULL) {
            Message("\n[CO2-UDF Error] Temperature header length mismatch in %s\n", filename);
            free(line);
            fclose(fp);
            free_table(table);
            return 0;
        }
        table->t_grid[j] = atof(token);
    }

    row_capacity = INITIAL_ROW_CAPACITY;
    row_count = 0;
    table->p_grid = (double*)malloc(sizeof(double) * row_capacity);
    table->values = (double*)malloc(sizeof(double) * row_capacity * nt);

    if (table->t_grid == NULL || table->p_grid == NULL || table->values == NULL) {
        Message("\n[CO2-UDF Error] Failed to allocate table arrays for %s\n", filename);
        free(line);
        fclose(fp);
        free_table(table);
        return 0;
    }

    while (fgets(line, CSV_MAX_LINE_CHARS, fp) != NULL) {
        token = strtok(line, ",\n\r");
        if (token == NULL) continue;

        if (row_count >= row_capacity) {
            row_capacity *= 2;
            table->p_grid = (double*)realloc(table->p_grid, sizeof(double) * row_capacity);
            new_values = (double*)realloc(table->values, sizeof(double) * row_capacity * nt);
            if (table->p_grid == NULL || new_values == NULL) {
                Message("\n[CO2-UDF Error] Failed to grow table arrays for %s\n", filename);
                free(line);
                fclose(fp);
                free_table(table);
                return 0;
            }
            table->values = new_values;
        }

        table->p_grid[row_count] = atof(token) * 1.0e6;
        for (j = 0; j < nt; j++) {
            token = strtok(NULL, ",\n\r");
            table->values[row_count * nt + j] = (token != NULL) ? atof(token) : CO2_MISSING_VALUE;
        }
        row_count++;
    }

    table->np = row_count;
    table->loaded = (table->np >= 2 && table->nt >= 2);
    free(line);
    fclose(fp);

    if (!table->loaded) {
        Message("\n[CO2-UDF Error] Insufficient rows in %s\n", filename);
        free_table(table);
        return 0;
    }

    table->t_min = table->t_grid[0];
    table->t_max = table->t_grid[table->nt - 1];
    table->p_min = table->p_grid[0];
    table->p_max = table->p_grid[table->np - 1];
    table->uniform_ready = grid_is_uniform(table->t_grid, table->nt, &table->t_step)
                        && grid_is_uniform(table->p_grid, table->np, &table->p_step);

    Message("\n[CO2-UDF] Loaded %s: np=%d, nt=%d, P=[%.6g, %.6g] Pa, T=[%.6g, %.6g] K\n",
            filename, table->np, table->nt, table->p_min, table->p_max, table->t_min, table->t_max);

    if (!table->uniform_ready) {
        Message("\n[CO2-UDF Error] %s is not a uniform grid table. Use the nonuniform UDF instead.\n", filename);
        return 0;
    }

    Message("[CO2-UDF] Uniform steps: dP=%.9g Pa, dT=%.9g K\n", table->p_step, table->t_step);
    return 1;
}

#if ENABLE_PHASE_DIAGNOSTICS
static int read_phase_table(const char *filename, PhaseTable2D *phase_table, const Table2D *reference)
{
    FILE *fp;
    char *line;
    char *token;
    int i, j;

    if (reference == NULL || !reference->uniform_ready) return 0;

    fp = fopen(filename, "r");
    if (fp == NULL) {
        Message("\n[CO2-UDF Warning] Cannot open %s. Phase diagnostics disabled.\n", filename);
        return 0;
    }

    line = (char*)malloc(CSV_MAX_LINE_CHARS);
    if (line == NULL) {
        fclose(fp);
        return 0;
    }

    free_phase_table(phase_table);
    phase_table->nt = reference->nt;
    phase_table->np = reference->np;
    phase_table->phase = (int*)malloc(sizeof(int) * reference->nt * reference->np);

    if (phase_table->phase == NULL) {
        free(line);
        fclose(fp);
        return 0;
    }

    fgets(line, CSV_MAX_LINE_CHARS, fp);
    for (i = 0; i < reference->np; i++) {
        if (fgets(line, CSV_MAX_LINE_CHARS, fp) == NULL) break;
        token = strtok(line, ",\n\r");
        for (j = 0; j < reference->nt; j++) {
            token = strtok(NULL, ",\n\r");
            phase_table->phase[i * reference->nt + j] = parse_phase_family(token);
        }
    }

    phase_table->loaded = 1;
    free(line);
    fclose(fp);
    Message("\n[CO2-UDF] Loaded %s for phase diagnostics.\n", filename);
    return 1;
}

static int cell_crosses_phase_boundary(const PhaseTable2D *phase_table, int ip, int it)
{
    int nt, p00, p10, p01, p11;
    if (phase_table == NULL || !phase_table->loaded) return 0;

    nt = phase_table->nt;
    p00 = phase_table->phase[ip * nt + it];
    p10 = phase_table->phase[ip * nt + it + 1];
    p01 = phase_table->phase[(ip + 1) * nt + it];
    p11 = phase_table->phase[(ip + 1) * nt + it + 1];

    return !(p00 == p10 && p00 == p01 && p00 == p11);
}
#endif

static double interpolate_table_uniform(double t, double p, const Table2D *table, const char *property_name)
{
    int it, ip;
    double calc_t, calc_p;
    double t0, p0;
    double u, v;
    double f00, f10, f01, f11;
    int nt;

    if (table == NULL || !table->uniform_ready) {
        return 0.0;
    }

    calc_t = t;
    calc_p = p;

    if (calc_t <= table->t_min) {
        calc_t = table->t_min;
#if ENABLE_RANGE_WARNINGS
        if (range_warning_count < WARNING_LIMIT) {
            Message("\n[CO2-UDF Warning] %s T below table range: %.6g K\n", property_name, t);
            range_warning_count++;
        }
#endif
    } else if (calc_t >= table->t_max) {
        calc_t = table->t_max;
#if ENABLE_RANGE_WARNINGS
        if (range_warning_count < WARNING_LIMIT) {
            Message("\n[CO2-UDF Warning] %s T above table range: %.6g K\n", property_name, t);
            range_warning_count++;
        }
#endif
    }

    if (calc_p <= table->p_min) {
        calc_p = table->p_min;
#if ENABLE_RANGE_WARNINGS
        if (range_warning_count < WARNING_LIMIT) {
            Message("\n[CO2-UDF Warning] %s P below table range: %.6g Pa\n", property_name, p);
            range_warning_count++;
        }
#endif
    } else if (calc_p >= table->p_max) {
        calc_p = table->p_max;
#if ENABLE_RANGE_WARNINGS
        if (range_warning_count < WARNING_LIMIT) {
            Message("\n[CO2-UDF Warning] %s P above table range: %.6g Pa\n", property_name, p);
            range_warning_count++;
        }
#endif
    }

    it = (int)floor((calc_t - table->t_min) / table->t_step);
    ip = (int)floor((calc_p - table->p_min) / table->p_step);

    if (it < 0) it = 0;
    if (ip < 0) ip = 0;
    if (it >= table->nt - 1) it = table->nt - 2;
    if (ip >= table->np - 1) ip = table->np - 2;

#if ENABLE_PHASE_DIAGNOSTICS
    if (cell_crosses_phase_boundary(&co2_phase, ip, it) && phase_warning_count < WARNING_LIMIT) {
        Message("\n[CO2-UDF Warning] %s interpolation cell crosses phase boundary near T=%.6g K, P=%.6g Pa.\n",
                property_name, t, p);
        phase_warning_count++;
    }
#endif

    t0 = table->t_min + it * table->t_step;
    p0 = table->p_min + ip * table->p_step;
    u = (calc_t - t0) / table->t_step;
    v = (calc_p - p0) / table->p_step;

    if (u < 0.0) u = 0.0;
    if (u > 1.0) u = 1.0;
    if (v < 0.0) v = 0.0;
    if (v > 1.0) v = 1.0;

    nt = table->nt;
    f00 = table->values[ip * nt + it];
    f10 = table->values[ip * nt + it + 1];
    f01 = table->values[(ip + 1) * nt + it];
    f11 = table->values[(ip + 1) * nt + it + 1];

    return (1.0 - u) * (1.0 - v) * f00
         + u * (1.0 - v) * f10
         + (1.0 - u) * v * f01
         + u * v * f11;
}

static double interpolate_cp_uniform(double t)
{
    int it;
    double calc_t;
    double t0;
    double u;

    if (cp_ref == NULL || cp_ref_nt < 2 || cp_t_step <= 0.0) return 0.0;

    calc_t = t;
    if (calc_t <= cp_t_min) calc_t = cp_t_min;
    else if (calc_t >= cp_t_max) calc_t = cp_t_max;

    it = (int)floor((calc_t - cp_t_min) / cp_t_step);
    if (it < 0) it = 0;
    if (it >= cp_ref_nt - 1) it = cp_ref_nt - 2;

    t0 = cp_t_min + it * cp_t_step;
    u = (calc_t - t0) / cp_t_step;
    if (u < 0.0) u = 0.0;
    if (u > 1.0) u = 1.0;

    return cp_ref[it] * (1.0 - u) + cp_ref[it + 1] * u;
}

static int nearest_pressure_row(const Table2D *table, double p_ref)
{
    int row;
    double x;

    if (table == NULL || !table->uniform_ready) return 0;
    x = (p_ref - table->p_min) / table->p_step;
    row = (int)floor(x + 0.5);
    if (row < 0) row = 0;
    if (row >= table->np) row = table->np - 1;
    return row;
}

static void build_cp_reference_curve(void)
{
    int i, row;

    if (!co2_cp_2d.uniform_ready) return;

    if (cp_ref) free(cp_ref);
    cp_ref_nt = co2_cp_2d.nt;
    cp_t_min = co2_cp_2d.t_min;
    cp_t_max = co2_cp_2d.t_max;
    cp_t_step = co2_cp_2d.t_step;
    cp_ref = (double*)malloc(sizeof(double) * cp_ref_nt);

    if (cp_ref == NULL) {
        cp_ref_nt = 0;
        return;
    }

    row = nearest_pressure_row(&co2_cp_2d, CP_REFERENCE_PRESSURE_PA);
    for (i = 0; i < cp_ref_nt; i++) {
        cp_ref[i] = co2_cp_2d.values[row * co2_cp_2d.nt + i];
    }

    Message("\n[CO2-UDF] Cp(T) reference curve uses nearest pressure row: %.6g Pa.\n",
            co2_cp_2d.p_grid[row]);
}

DEFINE_ON_DEMAND(CO2_import_property_tables_uniform)
{
    operating_pressure = RP_Get_Real("operating-pressure");
    range_warning_count = 0;
    phase_warning_count = 0;

    read_table_2d("co2_density.csv", &co2_density);
    read_table_2d("co2_viscosity.csv", &co2_viscosity);
    read_table_2d("co2_conductivity.csv", &co2_conductivity);
    read_table_2d("co2_cp.csv", &co2_cp_2d);

#if ENABLE_PHASE_DIAGNOSTICS
    read_phase_table("co2_phase.csv", &co2_phase, &co2_density);
#else
    Message("\n[CO2-UDF] Phase diagnostics disabled by ENABLE_PHASE_DIAGNOSTICS=0.\n");
#endif

    build_cp_reference_curve();
    Message("\n[CO2-UDF] Uniform import complete. Operating pressure = %.6g Pa.\n", operating_pressure);
}

DEFINE_PROPERTY(CO2_uniform_density, c, thread)
{
    double t = C_T(c, thread);
    double p_abs = C_P(c, thread) + operating_pressure;
    return interpolate_table_uniform(t, p_abs, &co2_density, "density");
}

DEFINE_PROPERTY(CO2_uniform_viscosity, c, thread)
{
    double t = C_T(c, thread);
    double p_abs = C_P(c, thread) + operating_pressure;
    return interpolate_table_uniform(t, p_abs, &co2_viscosity, "viscosity");
}

DEFINE_PROPERTY(CO2_uniform_thermal_conductivity, c, thread)
{
    double t = C_T(c, thread);
    double p_abs = C_P(c, thread) + operating_pressure;
    return interpolate_table_uniform(t, p_abs, &co2_conductivity, "thermal conductivity");
}

DEFINE_SPECIFIC_HEAT(CO2_uniform_cp, T, Tref, h, yi)
{
    double cp = interpolate_cp_uniform(T);
    *h = cp * (T - Tref);
    return cp;
}

DEFINE_ON_DEMAND(CO2_test_property_tables_uniform)
{
    double test_t = 623.15;
    double test_p = 14.0e6;
    double rho, mu, k, cp;

    if (!co2_density.uniform_ready || cp_ref == NULL) {
        Message("\n[CO2-UDF Warning] Please run CO2_import_property_tables_uniform first with uniform tables.\n");
        return;
    }

    rho = interpolate_table_uniform(test_t, test_p, &co2_density, "density");
    mu = interpolate_table_uniform(test_t, test_p, &co2_viscosity, "viscosity");
    k = interpolate_table_uniform(test_t, test_p, &co2_conductivity, "thermal conductivity");
    cp = interpolate_cp_uniform(test_t);

    Message("\n============ CO2 uniform table test ============\n");
    Message("Test point: T = %.3f K, P = %.6f MPa\n", test_t, test_p / 1.0e6);
    Message("Density              = %.9g kg/m3\n", rho);
    Message("Viscosity            = %.9g Pa s\n", mu);
    Message("Thermal conductivity = %.9g W/m/K\n", k);
    Message("Cp(reference P)      = %.9g J/kg/K\n", cp);
    Message("================================================\n\n");
}

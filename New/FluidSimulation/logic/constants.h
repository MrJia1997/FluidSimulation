#define MATH_PI 3.14159265358979323846

#define DELTA_T 0.016
#define SOLVER_ITERATIONS 4

#define ADJ_DISTANCE 0.1
#define KERNEL_H ADJ_DISTANCE * 1.5
#define EPS 1e-7

#define REST_DENSITY 960
#define CFM_EPSILON 300.0
// #define DRAG_COEFF 1

#define TENSILE_DELTA_Q (0.2 * KERNEL_H)
#define TENSILE_K 0.001
#define TENSILE_N 4

#define VISCOSITY_C 0.01
#define VORTICITY_EPSILON 0.0006
#define VORTICITY_DELTA KERNEL_H / 10.0

#define X_MIN -1.0
#define X_MAX 1.0
#define Y_MIN -1.5
#define Y_MAX 1.5
#define Z_MIN 0.0
#define Z_MAX 1.0

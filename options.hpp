#ifndef DDA_OPTIONS_H
#define DDA_OPTIONS_H

#include <limits>
#include <string>
#include <vector>

#include "complex.hpp"      /* real, cplx */

/* Program version, printed by --version. */
#define DDA_VERSION "7.1"

/* Values that mark an option as "not given". After parse_options() returns
   PARSE_OK every field of sim_options holds a concrete value. */

constexpr double OPT_UNSET_DOUBLE = std::numeric_limits<double>::quiet_NaN();
constexpr int OPT_UNSET_INT = std::numeric_limits<int>::min();

/* Configuration of one run. The fields correspond one-to-one to the
   variables that main() used to read with atof()/atoi() and scanf(). */

struct sim_options
{                               /* phase angles, degrees                     */
    double theta_min = OPT_UNSET_DOUBLE;
    double theta_max = OPT_UNSET_DOUBLE;
    double theta_step = OPT_UNSET_DOUBLE;

                                /* files                                     */
    std::string state_file;     /* binary state, resumed if it exists        */
    std::string output_file;    /* text result                               */
    bool restart = false;       /* delete the state file before starting     */

                                /* samples and orientation                   */
    int samples = OPT_UNSET_INT;
    int orientation = OPT_UNSET_INT;    /* 0 fixed, 1 Monte-Carlo            */
    double psi = OPT_UNSET_DOUBLE;      /* degrees, fixed orientation        */
    double iota = OPT_UNSET_DOUBLE;
    double kappa = OPT_UNSET_DOUBLE;
    int planes = OPT_UNSET_INT;         /* scattering planes                 */

                                /* target                                    */
    int target = OPT_UNSET_INT;         /* 1..14                             */
    int grid[3] = {0, 0, 0};            /* i_max, j_max, l_max               */
    double size_parameter = OPT_UNSET_DOUBLE;
    double aspect_b = OPT_UNSET_DOUBLE; /* b/a                               */
    double aspect_c = OPT_UNSET_DOUBLE; /* c/a                               */
    double radius_ratio = OPT_UNSET_DOUBLE;     /* R1/R2, target 12          */
    int hierarchical_levels = OPT_UNSET_INT;    /* target 13                 */
    int bites = OPT_UNSET_INT;                  /* target 10                 */
    double depth = OPT_UNSET_DOUBLE;            /* target 10                 */
    double packing = OPT_UNSET_DOUBLE;          /* rho, targets 6,7,8,14     */
    double monomer_radius = OPT_UNSET_DOUBLE;   /* targets 6,7,8             */
    double monomer_length_ratio = OPT_UNSET_DOUBLE;
    std::string shape_file;                     /* target 11                 */
    std::vector<cplx> materials;                /* refractive indices        */
    std::vector<int> seeds;                     /* per-material seeds        */

                                /* numerics                                  */
    int polarizability = OPT_UNSET_INT; /* 1..5                              */
    double tolerance = OPT_UNSET_DOUBLE;        /* Error_max                 */
    int max_steps = OPT_UNSET_INT;              /* Max_Step                  */
    int fft_mode = OPT_UNSET_INT;               /* 0 Cooley-Tukey, 1 R-B     */
    int solver = OPT_UNSET_INT;                 /* 0..3                      */
    double tau = OPT_UNSET_DOUBLE;              /* weight factor, solver 3   */
    int threads = OPT_UNSET_INT;                /* 0 = auto                  */

                                /* reproducibility                           */
    bool has_seed = false;
    unsigned long long seed = 0;
};

enum parse_status
{
    PARSE_OK = 0,       /* options resolved; run the simulation   */
    PARSE_EXIT_OK = 1,  /* --help or --version printed; exit 0    */
    PARSE_ERROR = 2     /* invalid usage; the message is printed  */
};

/* Fills *opt from the command line. When the program is started without
   arguments, or with --interactive, the values that are still missing are
   asked for on the standard input. */

int parse_options(int argc, char **argv, sim_options *opt);

#endif

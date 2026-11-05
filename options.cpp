/****************************************************************************/
/************************** Command-line arguments **************************/
/****************************************************************************/

/* The program is configured with named options (argtable3). Everything that
   is not given on the command line is asked for interactively when the
   program is started without arguments or with --interactive.

   argtable3 takes exactly one value per option occurrence, so the list-like
   options (--material, --grid, --aspect, --seeds) accept either comma
   separated values or several values after a single option:
       --material 1.5 0.01      --material 1.5,0.01
       --grid 64 64 64          --grid 64,64,64
   The second form is folded into the first before arg_parse() is called. */

#include "options.hpp"

#include <cerrno>
#include <climits>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <strings.h>            /* strcasecmp */

#include "third_party/argtable3/argtable3.h"

/****************************************************************************/
/****************************** Name tables *********************************/
/****************************************************************************/

typedef
    struct
    {   const char * name;
        int code;
    }
    name_code;

static const name_code target_names[] =
{   {"sphere",                 1},
    {"ellipsoid",              2},
    {"parallelepiped",         3},
    {"random-sphere",          4},
    {"random-ellipsoid",       5},
    {"anthrax",                6},
    {"hedgehog",               7},
    {"dead-anthrax",           8},
    {"random-parallelepiped",  9},
    {"nibble-ellipsoid",      10},
    {"file",                  11},
    {"bisphere",              12},
    {"hierarchical",          13},
    {"sphere-spec",           14},
    {NULL, 0}
};

static const name_code orientation_names[] =
{   {"fixed",         0},
    {"mc",            1},
    {"monte-carlo",   1},
    {NULL, 0}
};

static const name_code polarizability_names[] =
{   {"cm",       1},
    {"cmrr",     2},
    {"cmnr",     3},
    {"ildr",     4},
    {"ldr",      5},
    {NULL, 0}
};

static const name_code fft_names[] =
{   {"ct",            0},
    {"cooley-tukey",  0},
    {"rb",            1},
    {"rader-brenner", 1},
    {NULL, 0}
};

static const name_code solver_names[] =
{   {"cg",       0},
    {"aa",       1},
    {"bicgstab", 2},
    {"simple",   3},
    {NULL, 0}
};

/****************************************************************************/
/************************** Value conversions *******************************/
/****************************************************************************/

/* An integer, or one of the names in the table. */

static int parse_code(const name_code * names, const char * s, bool * ok)
{   char * end;
    long v;
    int i;

    *ok=false;
    if(s==NULL) return OPT_UNSET_INT;

    errno=0;
    v=strtol(s,&end,10);
    if(end!=s&&*end==0&&errno==0)
    {   *ok=true;
        return (int)v;
    }

    for(i=0;names[i].name!=NULL;i++)
        if(strcasecmp(names[i].name,s)==0)
        {   *ok=true;
            return names[i].code;
        }

    return OPT_UNSET_INT;
}

/* Comma separated lists of numbers, used by --material, --grid, --aspect and
   --seeds. The values are trimmed and must be complete numbers. */

static bool split_values(const char * s, std::vector<double> * out)
{   const char * p=s;

    out->clear();
    if(s==NULL||*s==0) return false;

    for(;;)
    {   const char * comma=strchr(p,',');
        size_t len=(comma==NULL)?strlen(p):(size_t)(comma-p);
        std::string item(p,len);
        size_t b=item.find_first_not_of(" \t");
        size_t e=item.find_last_not_of(" \t");
        char * end;
        double v;

        if(b==std::string::npos) return false;
        item=item.substr(b,e-b+1);

        errno=0;
        v=strtod(item.c_str(),&end);
        if(end==item.c_str()||*end!=0||errno!=0) return false;
        out->push_back(v);

        if(comma==NULL) break;
        p=comma+1;
    }

    return !out->empty();
}

static bool split_ints(const char * s, std::vector<int> * out)
{   const char * p=s;

    out->clear();
    if(s==NULL||*s==0) return false;

    for(;;)
    {   const char * comma=strchr(p,',');
        size_t len=(comma==NULL)?strlen(p):(size_t)(comma-p);
        std::string item(p,len);
        size_t b=item.find_first_not_of(" \t");
        size_t e=item.find_last_not_of(" \t");
        char * end;
        long v;

        if(b==std::string::npos) return false;
        item=item.substr(b,e-b+1);

        errno=0;
        v=strtol(item.c_str(),&end,10);
        if(end==item.c_str()||*end!=0||errno!=0) return false;
        if(v<INT_MIN||v>INT_MAX) return false;
        out->push_back((int)v);

        if(comma==NULL) break;
        p=comma+1;
    }

    return !out->empty();
}

/****************************************************************************/
/******************* Folding of the list-option values **********************/
/****************************************************************************/

/* argtable3 consumes exactly one value per occurrence of an option. To keep
   "--grid 64 64 64" working, the numeric values that follow such an option
   are glued into one comma separated token. */

typedef
    struct
    {   const char * token;      /* what the user types                     */
        const char * emit;       /* long form used in the folded argument   */
        int min_values;
        int max_values;
    }
    list_option;

static const list_option list_options[] =
{   {"-m",        "--material",  2,  2},
    {"--material", "--material", 2,  2},
    {"-g",        "--grid",      1,  3},
    {"--grid",    "--grid",      1,  3},
    {"--aspect",  "--aspect",    2,  2},
    {"--seeds",   "--seeds",     1, 64},
    {NULL, NULL, 0, 0}
};

static const list_option * find_list_option(const char * token)
{   int i;

    for(i=0;list_options[i].token!=NULL;i++)
        if(strcmp(list_options[i].token,token)==0)
            return &list_options[i];

    return NULL;
}

static bool looks_numeric(const char * s)
{   char * end;

    if(s==NULL||*s==0) return false;
    errno=0;
    (void)strtod(s,&end);

    return end!=s&&*end==0&&errno==0;
}

/* Builds the argv that is handed to arg_parse(). The strings live in *storage
   and the pointers in *out_argv; both must outlive the arg_parse() call. */

static void fold_list_values(int argc, char ** argv,
                             std::vector<std::string> * storage,
                             std::vector<char *> * out_argv)
{   int i,j;

    storage->clear();
    out_argv->clear();
    storage->reserve((size_t)argc+1);
    out_argv->reserve((size_t)argc+1);

    storage->push_back(std::string(argv[0]));

    for(i=1;i<argc;)
    {   const list_option * lo=find_list_option(argv[i]);

        if(lo==NULL)
        {   storage->push_back(std::string(argv[i]));
            i++;
            continue;
        }

        {   std::string values;
            int taken=0;

            for(j=i+1;j<argc&&taken<lo->max_values;j++)
            {   if(!looks_numeric(argv[j])) break;
                if(taken>0) values+=',';
                values+=argv[j];
                taken++;
            }

            if(taken>=lo->min_values)
            {   storage->push_back(std::string(lo->emit)+"="+values);
                i+=taken+1;
            }
            else
            {   storage->push_back(std::string(argv[i]));
                i++;
            }
        }
    }

    for(i=0;i<(int)storage->size();i++)
        out_argv->push_back(&(*storage)[(size_t)i][0]);
}

/****************************************************************************/
/*************************** Interactive input ******************************/
/****************************************************************************/

static bool read_line(char * buf, size_t size)
{   size_t len;

    if(fgets(buf,(int)size,stdin)==NULL)
    {   buf[0]=0;
        return false;
    }

    len=strlen(buf);
    while(len>0&&(buf[len-1]=='\n'||buf[len-1]=='\r'))
        buf[--len]=0;

    return true;
}

/* def == OPT_UNSET_DOUBLE means "no default, a value is required". */

static double ask_double(const char * prompt, double def, double given)
{   char buf[256];
    bool has_default=!std::isnan(def);

    if(!std::isnan(given)) return given;

    for(;;)
    {   if(has_default) printf("\n %s [%g] ",prompt,def);
        else            printf("\n %s ",prompt);
        fflush(stdout);

        if(!read_line(buf,sizeof(buf))||buf[0]==0)
        {   if(has_default) return def;
            printf(" (a value is required)");
            continue;
        }

        {   char * end;
            double v;

            errno=0;
            v=strtod(buf,&end);
            if(end!=buf&&*end==0&&errno==0) return v;
        }

        printf(" (expected a number)");
    }
}

/* def == OPT_UNSET_INT means "no default, a value is required". */

static int ask_int(const char * prompt, int def, int given)
{   char buf[256];
    bool has_default=(def!=OPT_UNSET_INT);

    if(given!=OPT_UNSET_INT) return given;

    for(;;)
    {   if(has_default) printf("\n %s [%d] ",prompt,def);
        else            printf("\n %s ",prompt);
        fflush(stdout);

        if(!read_line(buf,sizeof(buf))||buf[0]==0)
        {   if(has_default) return def;
            printf(" (a value is required)");
            continue;
        }

        {   char * end;
            long v;

            errno=0;
            v=strtol(buf,&end,10);
            if(end!=buf&&*end==0&&errno==0&&v>=INT_MIN&&v<=INT_MAX)
                return (int)v;
        }

        printf(" (expected an integer)");
    }
}

static std::string ask_string(const char * prompt, const char * def,
                              const std::string & given)
{   char buf[1024];

    if(!given.empty()) return given;

    printf("\n %s [%s] ",prompt,def);
    fflush(stdout);

    if(!read_line(buf,sizeof(buf))||buf[0]==0) return std::string(def);

    return std::string(buf);
}

/* def == OPT_UNSET_INT means "no default, a value is required". */

static int ask_choice(const char * prompt, const name_code * names, int def,
                      int given)
{   char buf[256];

    if(given!=OPT_UNSET_INT) return given;

    for(;;)
    {   printf("\n %s",prompt);
        if(def!=OPT_UNSET_INT) printf(" [%d]",def);
        printf(" ");
        fflush(stdout);

        if(!read_line(buf,sizeof(buf))||buf[0]==0)
        {   if(def!=OPT_UNSET_INT) return def;
            printf(" (a value is required)");
            continue;
        }

        {   bool ok=false;
            int v=parse_code(names,buf,&ok);

            if(ok) return v;
        }

        printf(" (expected a number or one of the names)");
    }
}

static const char * target_menu=
    "\n Please, choose the morphology type of target"
    "\n Input  1 for sphere"
    "\n Input  2 for ellipsoid"
    "\n Input  3 for parallelepiped"
    "\n Input  4 for random particle in spherical matrix"
    "\n Input  5 for random particle in ellipsoidal matrix"
    "\n Input  6 for random antrax-like particle"
    "\n Input  7 for random hedgehog-like particle"
    "\n Input  8 for random dead antrax-like particle"
    "\n Input  9 for random particle in parallelepiped matrix"
    "\n Input 10 for random nibble particle in ellipsoidal matrix"
    "\n Input 11 for particle from file"
    "\n Input 12 for two spheres in contact (along X-axis)"
    "\n Input 13 for hierarchical random particle in ellipsoidal matrix"
    "\n Input 14 for sphere with arbitrary radius";

static const char * polarizability_menu=
    "\n Please, choose the relation between polarizability and refractive index."
    "\n Input 1 for Clausius-Mossotti"
    "\n Input 2 for Clausius-Mossotti plus radiative reaction (by Draine)"
    "\n Input 3 for Clausius-Mossotti plus nonradiative reaction (by Goedecke&O'Brien)"
    "\n Input 4 for Isotropized Lattice Dispersion Relation (by Draine)"
    "\n Input 5 for Lattice Dispersion Relation (by Draine)"
    "\n Note: Do not choose 4 if the wave-vector of incident field is parallel"
    "\n to one of the coordinate axes, as it reduces the accuracy of computations!";

static const char * fft_menu=
    "\n Please, choose the algorithm of the FFT"
    "\n Input 0 for the Cooley-Tukey method"
    "\n Input 1 for the Rader-Brenner method";

static const char * solver_menu=
    "\n Please, choose the solution method of system of linear algebraic equations"
    "\n Input 0 for Conjugate Gradient Method (by Voevodin)"
    "\n Input 1 for AA* Minimal Iterations Method (by Voevodin)"
    "\n Input 2 for BiConjugate Gradient Stabilized Method"
    "\n Input 3 for Simple Iteration Method";

/****************************************************************************/
/****************************** Validation **********************************/
/****************************************************************************/

static bool target_needs_aspect(int t)
{   return t==2||t==3||t==5||t==9||t==10||t==13;
}

static bool target_needs_seeds(int t)
{   return t==4||t==5||t==9||t==10||t==13;
}

static bool target_is_worm(int t)
{   return t==6||t==7||t==8;
}

/* Targets built from a single constituent material. The random particles
   (4,5,9,10), the file target (11) and the two-material targets (12,13) are
   handled separately. */

static bool target_single_material(int t)
{   return t==1||t==2||t==3||t==6||t==7||t==8||t==14;
}

static bool is_power_of_two(int n)
{   return n>0&&(n&(n-1))==0;
}

static std::string join_errors(const std::vector<std::string> & errors)
{   std::string res;
    size_t i;

    for(i=0;i<errors.size();i++)
    {   if(i>0) res+="\n";
        res+="  ";
        res+=errors[i];
    }

    return res;
}

static std::string validate(const sim_options & o)
{   std::vector<std::string> errors;
    char buf[256];
    int i;

    if(!(o.theta_step>0.))
        errors.push_back("--theta-step must be positive");
    if(o.theta_max<o.theta_min)
        errors.push_back("--theta-max must not be less than --theta-min");
    if(o.samples<1)
        errors.push_back("--samples must be at least 1");
    if(o.planes<1)
        errors.push_back("--planes must be at least 1");
    if(o.orientation!=0&&o.orientation!=1)
        errors.push_back("--orientation must be fixed or mc");
    if(o.target<1||o.target>14)
        errors.push_back("--target must be 1..14 or a target name");

    {   int bad_sign=0,bad_pow2=0;

        for(i=0;i<3;i++)
        {   if(o.grid[i]<=0)                        bad_sign=1;
            else if(!is_power_of_two(o.grid[i]))    bad_pow2=1;
        }

        if(bad_sign)
            errors.push_back("--grid needs 1 or 3 positive dimensions");
        else if(bad_pow2)
            errors.push_back("--grid dimensions must be powers of two");
    }

    if(!(o.size_parameter>0.))
        errors.push_back("--size-parameter must be positive");
    if(o.polarizability<1||o.polarizability>5)
        errors.push_back("--polarizability must be 1..5");
    if(!(o.tolerance>0.))
        errors.push_back("--tolerance must be positive");
    if(o.max_steps<1)
        errors.push_back("--max-steps must be at least 1");
    if(o.threads<0)
        errors.push_back("--threads must be 0 (auto) or positive");

    if(o.output_file.empty())
        errors.push_back("--output must not be empty");
    if(o.state_file.empty())
        errors.push_back("--state must not be empty");
    if(o.output_file==o.state_file)
        errors.push_back("--output and --state must name different files");
    if(o.output_file.size()>1023)
        errors.push_back("--output is too long (1024 bytes at most)");
    if(o.state_file.size()>1023)
        errors.push_back("--state is too long (1024 bytes at most)");

    if(o.solver==3&&!(o.tau>0.))
        errors.push_back("solver 3 needs --tau > 0");

    /* Target specific parameters. */

    if(o.target==11)
    {   if(o.shape_file.empty())
            errors.push_back("target 11 needs --shape-file");
        if(o.shape_file.size()>1023)
            errors.push_back("--shape-file is too long (1024 bytes at most)");
        if(o.materials.empty())
            errors.push_back("target 11 needs at least one --material");
    }
    else if(o.target==12||o.target==13)
    {   if(o.materials.size()!=2)
        {   snprintf(buf,sizeof(buf),
                     "target %d needs exactly two --material options",
                     o.target);
            errors.push_back(buf);
        }
    }
    else if(target_single_material(o.target))
    {   if(o.materials.size()!=1)
            errors.push_back("this target needs exactly one --material");
    }
    else if(o.target>=4&&o.target<=10)
    {   if(o.materials.empty())
            errors.push_back("this target needs at least one --material");
    }

    if(target_needs_seeds(o.target)
       &&o.seeds.size()!=o.materials.size())
    {   snprintf(buf,sizeof(buf),
                 "--seeds needs one value per material (%d expected, %d given)",
                 (int)o.materials.size(),(int)o.seeds.size());
        errors.push_back(buf);
    }

    if(target_needs_aspect(o.target))
    {   if(std::isnan(o.aspect_b)||!(o.aspect_b>0.))
            errors.push_back("this target needs --aspect <b/a,c/a>");
        if(std::isnan(o.aspect_c)||!(o.aspect_c>0.))
            errors.push_back("this target needs --aspect <b/a,c/a>");
    }

    if(o.target==12&&!(o.radius_ratio>0.))
        errors.push_back("target 12 needs --radius-ratio > 0");

    if(o.target==13&&o.hierarchical_levels<1)
        errors.push_back("target 13 needs --levels >= 1");

    if(o.target==10)
    {   if(o.bites<0)
            errors.push_back("target 10 needs --bites >= 0");
        if(std::isnan(o.depth)||!(o.depth>=0.))
            errors.push_back("target 10 needs --depth >= 0");
    }

    if(target_is_worm(o.target)||o.target==14)
    {   if(std::isnan(o.packing)||!(o.packing>=0.))
            errors.push_back("this target needs --packing >= 0");
    }

    if(target_is_worm(o.target))
    {   if(std::isnan(o.monomer_radius)||!(o.monomer_radius>0.))
            errors.push_back("this target needs --monomer-radius > 0");
        if(std::isnan(o.monomer_length_ratio)||!(o.monomer_length_ratio>=0.))
            errors.push_back("this target needs --monomer-length-ratio >= 0");
    }

    /* Options that this run ignores are reported but do not stop it. */

    if(o.solver!=3&&o.tau>0.)
        fprintf(stderr,"warning: --tau is only used by solver 3, ignored\n");

    if(o.orientation==1
       &&(o.psi!=0.||o.iota!=0.||o.kappa!=0.))
        fprintf(stderr,"warning: --psi/--iota/--kappa are only used by "
                       "--orientation fixed, ignored\n");

    return join_errors(errors);
}

/****************************************************************************/
/****************************** Defaults ************************************/
/****************************************************************************/

static void apply_defaults(sim_options * o)
{   if(std::isnan(o->theta_min))    o->theta_min=0.;
    if(std::isnan(o->theta_max))    o->theta_max=180.;
    if(std::isnan(o->theta_step))   o->theta_step=10.;

    if(o->state_file.empty())       o->state_file="out.dst";
    if(o->output_file.empty())      o->output_file="out.dat";

    if(o->samples==OPT_UNSET_INT)   o->samples=1;
    if(o->orientation==OPT_UNSET_INT)   o->orientation=0;
    if(std::isnan(o->psi))          o->psi=0.;
    if(std::isnan(o->iota))         o->iota=0.;
    if(std::isnan(o->kappa))        o->kappa=0.;
    if(o->planes==OPT_UNSET_INT)    o->planes=1;

    if(std::isnan(o->aspect_b))     o->aspect_b=1.;
    if(std::isnan(o->aspect_c))     o->aspect_c=1.;
    if(std::isnan(o->radius_ratio)) o->radius_ratio=1.;
    if(o->hierarchical_levels==OPT_UNSET_INT)   o->hierarchical_levels=2;
    if(o->bites==OPT_UNSET_INT)     o->bites=0;
    if(std::isnan(o->depth))        o->depth=0.;

    if(o->polarizability==OPT_UNSET_INT)    o->polarizability=1;
    if(std::isnan(o->tolerance))    o->tolerance=1e-3;
    if(o->max_steps==OPT_UNSET_INT) o->max_steps=20;
    if(o->fft_mode==OPT_UNSET_INT)  o->fft_mode=0;
    if(o->solver==OPT_UNSET_INT)    o->solver=0;
    if(std::isnan(o->tau))          o->tau=0.;
    if(o->threads==OPT_UNSET_INT)   o->threads=0;
}

/****************************************************************************/
/*************************** Interactive prompts ****************************/
/****************************************************************************/

static void prompt_target_parameters(sim_options * o)
{   int i;

    if(o->materials.empty())
    {   int n;

        if(o->target==12||o->target==13)
            n=2;
        else
            n=ask_int("Input number of the constituent materials",1,OPT_UNSET_INT);

        if(n<1) n=1;

        for(i=0;i<n;i++)
        {   cplx m;
            char prompt[128];

            snprintf(prompt,sizeof(prompt),
                     "Input the real part of refractive index of %d material",
                     i+1);
            m.r=(real)ask_double(prompt,OPT_UNSET_DOUBLE,OPT_UNSET_DOUBLE);

            snprintf(prompt,sizeof(prompt),
                     "Input the imaginary part of refractive index of %d material",
                     i+1);
            m.i=(real)ask_double(prompt,OPT_UNSET_DOUBLE,OPT_UNSET_DOUBLE);

            o->materials.push_back(m);
        }
    }

    if(target_needs_seeds(o->target))
    {   while(o->seeds.size()<o->materials.size())
        {   char prompt[128];

            snprintf(prompt,sizeof(prompt),"Input number of %d material seeds",
                     (int)o->seeds.size()+1);
            o->seeds.push_back(ask_int(prompt,OPT_UNSET_INT,OPT_UNSET_INT));
        }
    }

    if(target_needs_aspect(o->target))
    {   o->aspect_b=ask_double("Input the aspect ratio b/a",1.,o->aspect_b);
        o->aspect_c=ask_double("Input the aspect ratio c/a",1.,o->aspect_c);
    }

    if(o->target==12)
    {   o->radius_ratio=ask_double("Input ratio of the first sphere radius "
                                   "over the second sphere radius",
                                   1.,o->radius_ratio);
    }

    if(o->target==13)
    {   o->hierarchical_levels=ask_int("Input number of hierarchical levels",
                                       2,o->hierarchical_levels);
    }

    if(o->target==10)
    {   o->bites=ask_int("Input number of the outer ceed cells for empty space",
                         0,o->bites);
        o->depth=ask_double("Input depth of the surface layer with these ceed "
                            "cells (in radii of the initial matrix)",
                            0.,o->depth);
    }

    if(target_is_worm(o->target))
    {   o->packing=ask_double("Input the volume density of the target particle",
                              OPT_UNSET_DOUBLE,o->packing);
        o->monomer_radius=ask_double("Input radius of the monomer (in cells)",
                                     OPT_UNSET_DOUBLE,o->monomer_radius);
        o->monomer_length_ratio=ask_double("Input the ratio length over radius "
                                           "for the monomer",
                                           OPT_UNSET_DOUBLE,
                                           o->monomer_length_ratio);
    }

    if(o->target==14)
        o->packing=ask_double("Input the density rho (the part of the initial "
                              "spherical volume filled by dipoles)",
                              OPT_UNSET_DOUBLE,o->packing);

    if(o->target==11)
        o->shape_file=ask_string("Input name of the file containing data on "
                                 "the target","",o->shape_file);
}

static void prompt_missing(sim_options * o)
{   o->theta_min=ask_double("Input the bottom limit of range of phase angles",
                             0.,o->theta_min);
    o->theta_max=ask_double("Input the upper limit of range of phase angles",
                            180.,o->theta_max);
    o->theta_step=ask_double("Input the increment of phase angles",
                             10.,o->theta_step);

    o->state_file=ask_string("Input name for the binary result file (dst)",
                             "out.dst",o->state_file);
    o->output_file=ask_string("Input name for the text result file (dat)",
                              "out.dat",o->output_file);

    o->samples=ask_int("Input number of samples/orientations",1,o->samples);

    o->orientation=ask_choice("Please, choose the orientation mode"
                              "\n 0 - no averaging"
                              "\n 1 - Monte-Carlo averaging (random orientations)",
                              orientation_names,0,o->orientation);

    if(o->orientation==0)
    {   o->psi=ask_double("Input the angle psi for the first rotation around "
                          "axis Z (in degree)",0.,o->psi);
        o->iota=ask_double("Input the angle iota for the second rotation around "
                           "axis Y (in degree)",0.,o->iota);
        o->kappa=ask_double("Input the angle kappa for the third rotation around "
                            "axis X (in degree)",0.,o->kappa);
    }

    o->planes=ask_int("Input number of the scattering planes",1,o->planes);

    o->target=ask_choice(target_menu,target_names,1,o->target);

    prompt_target_parameters(o);

    if(o->grid[0]<=0)
    {   o->grid[0]=ask_int("Input the maximum number of cells along the axis X "
                           "(should be 2^integer)",OPT_UNSET_INT,OPT_UNSET_INT);
        o->grid[1]=ask_int("Input the maximum number of cells along the axis Y "
                           "(should be 2^integer)",o->grid[0],OPT_UNSET_INT);
        o->grid[2]=ask_int("Input the maximum number of cells along the axis Z "
                           "(should be 2^integer)",o->grid[1],OPT_UNSET_INT);
    }

    o->size_parameter=ask_double("Input the size parameter x of the target "
                                 "(in respect to the side X)",
                                 OPT_UNSET_DOUBLE,o->size_parameter);

    o->polarizability=ask_choice(polarizability_menu,polarizability_names,1,
                                 o->polarizability);
    o->tolerance=ask_double("Input the parameter Error_max (it specifies the "
                            "termination of the iterative process)",
                            1e-3,o->tolerance);
    o->max_steps=ask_int("Input the maximum number of steps in the iterative "
                         "process",20,o->max_steps);
    o->fft_mode=ask_choice(fft_menu,fft_names,0,o->fft_mode);
    o->solver=ask_choice(solver_menu,solver_names,0,o->solver);

    if(o->solver==3)
        o->tau=ask_double("Input the weight factor for Simple Iteration Method",
                          OPT_UNSET_DOUBLE,o->tau);
}

/****************************************************************************/
/****************************** Parsing *************************************/
/****************************************************************************/

int parse_options(int argc, char ** argv, sim_options * opt)
{   arg_lit_t * o_help;
    arg_lit_t * o_version;
    arg_lit_t * o_interactive;
    arg_file_t * o_output;
    arg_file_t * o_state;
    arg_lit_t * o_restart;
    arg_int_t * o_samples;
    arg_dbl_t * o_tmin;
    arg_dbl_t * o_tmax;
    arg_dbl_t * o_tstep;
    arg_str_t * o_orientation;
    arg_dbl_t * o_psi;
    arg_dbl_t * o_iota;
    arg_dbl_t * o_kappa;
    arg_int_t * o_planes;
    arg_str_t * o_target;
    arg_str_t * o_grid;
    arg_dbl_t * o_size;
    arg_str_t * o_aspect;
    arg_dbl_t * o_rratio;
    arg_int_t * o_levels;
    arg_int_t * o_bites;
    arg_dbl_t * o_depth;
    arg_dbl_t * o_packing;
    arg_dbl_t * o_mradius;
    arg_dbl_t * o_mlr;
    arg_file_t * o_shape;
    arg_str_t * o_material;
    arg_str_t * o_seeds;
    arg_str_t * o_polarizability;
    arg_dbl_t * o_tolerance;
    arg_int_t * o_steps;
    arg_str_t * o_fft;
    arg_str_t * o_solver;
    arg_dbl_t * o_tau;
    arg_int_t * o_threads;
    arg_str_t * o_seed;
    arg_end_t * o_end;

    void * argtable[] =
    {   o_help=arg_lit0("h","help","print this help and exit"),
        o_version=arg_lit0("V","version","print version information and exit"),
        o_interactive=arg_lit0("i","interactive","ask for every value that is not given"),
        o_output=arg_file0("o","output","<file>","text result file (default out.dat)"),
        o_state=arg_file0("D","state","<file>","binary state file, resumed if it exists (default out.dst)"),
        o_restart=arg_lit0(NULL,"restart","delete the state file and start from scratch"),
        o_samples=arg_int0("n","samples","<N>","number of samples/orientations (default 1)"),
        o_tmin=arg_dbl0("a","theta-min","<deg>","lower limit of the phase angles (default 0)"),
        o_tmax=arg_dbl0("b","theta-max","<deg>","upper limit of the phase angles (default 180)"),
        o_tstep=arg_dbl0("s","theta-step","<deg>","phase-angle increment (default 10)"),
        o_orientation=arg_str0("O","orientation","<mode>","fixed|mc (default fixed)"),
        o_psi=arg_dbl0(NULL,"psi","<deg>","first rotation, fixed orientation (default 0)"),
        o_iota=arg_dbl0(NULL,"iota","<deg>","second rotation, fixed orientation (default 0)"),
        o_kappa=arg_dbl0(NULL,"kappa","<deg>","third rotation, fixed orientation (default 0)"),
        o_planes=arg_int0("p","planes","<N>","number of scattering planes (default 1)"),
        o_target=arg_str0("t","target","<name|1..14>","target type, see the target list"),
        o_grid=arg_str0("g","grid","<i[,j,l]>","lattice dimensions, powers of two; one value means a cube"),
        o_size=arg_dbl0("x","size-parameter","<x>","size parameter of the target"),
        o_aspect=arg_str0(NULL,"aspect","<b/a,c/a>","aspect ratios"),
        o_rratio=arg_dbl0(NULL,"radius-ratio","<R1/R2>","radius ratio, target 12"),
        o_levels=arg_int0(NULL,"levels","<N>","number of hierarchical levels, target 13"),
        o_bites=arg_int0(NULL,"bites","<N>","outer ceed cells for empty space, target 10"),
        o_depth=arg_dbl0(NULL,"depth","<r>","depth of the surface layer, target 10"),
        o_packing=arg_dbl0(NULL,"packing","<rho>","volume density, targets 6,7,8,14"),
        o_mradius=arg_dbl0(NULL,"monomer-radius","<r>","monomer radius in cells, targets 6,7,8"),
        o_mlr=arg_dbl0(NULL,"monomer-length-ratio","<l/r>","monomer length over radius, targets 6,7,8"),
        o_shape=arg_file0("f","shape-file","<file>","file with the target shape, target 11"),
        o_material=arg_strn("m","material","<mre,mi>",0,64,"refractive index of one material; repeat for a mixture"),
        o_seeds=arg_strn(NULL,"seeds","<seed[,seed...]>",0,64,"random seeds, one per material"),
        o_polarizability=arg_str0("P","polarizability","<model>","1..5 or cm,cmrr,cmnr,ildr,ldr (default 1)"),
        o_tolerance=arg_dbl0("e","tolerance","<err>","convergence threshold (default 1e-3)"),
        o_steps=arg_int0("S","max-steps","<N>","maximum number of iterations (default 20)"),
        o_fft=arg_str0("F","fft","<method>","0|1 or ct,rb (default 0)"),
        o_solver=arg_str0("M","solver","<method>","0..3 or cg,aa,bicgstab,simple (default 0)"),
        o_tau=arg_dbl0(NULL,"tau","<factor>","weight factor, required by solver 3"),
        o_threads=arg_int0("j","threads","<N>","worker threads, 0 = auto (default 0)"),
        o_seed=arg_str0(NULL,"seed","<u64>","random seed for reproducible runs"),
        o_end=arg_end(20)
    };

    const char * prog=strrchr(argv[0],'/');
    std::vector<std::string> fstorage;
    std::vector<char *> fargv;
    std::vector<std::string> errors;
    std::string err;
    int nerrors;
    int i;

    prog=(prog==NULL)?argv[0]:prog+1;

    *opt=sim_options();

    fold_list_values(argc,argv,&fstorage,&fargv);
    nerrors=arg_parse((int)fargv.size(),&fargv[0],argtable);

    if(o_help->count>0||o_version->count>0)
    {   if(o_help->count>0)
        {   void * common[]={o_help,o_version,o_interactive,o_end};
            void * run[]={o_output,o_state,o_restart,o_samples,o_tmin,o_tmax,
                          o_tstep,o_orientation,o_psi,o_iota,o_kappa,o_planes,
                          o_end};
            void * target[]={o_target,o_grid,o_size,o_aspect,o_rratio,o_levels,
                             o_bites,o_depth,o_packing,o_mradius,o_mlr,o_shape,
                             o_material,o_seeds,o_end};
            void * numerics[]={o_polarizability,o_tolerance,o_steps,o_fft,
                               o_solver,o_tau,o_threads,o_seed,o_end};

            printf("\nUsage: %s [options]\n",prog);
            printf("\n DDA-FFT computes light scattering by an arbitrary target with the"
                   "\n discrete dipole approximation. Every option is optional: the values"
                   "\n that are missing are asked for interactively when the program is"
                   "\n started without arguments (or with --interactive)."
                   "\n\n List options accept '--option a b' and '--option a,b'."
                   "\n The worker count defaults to $DDA_THREADS, or to one per core when"
                   "\n the variable is not set.\n");

            printf("\nCommon options:\n");
            arg_print_glossary(stdout,common,"  %-32s %s\n");
            printf("\nRun and output:\n");
            arg_print_glossary(stdout,run,"  %-32s %s\n");
            printf("\nTarget:\n");
            arg_print_glossary(stdout,target,"  %-32s %s\n");
            printf("\nNumerics:\n");
            arg_print_glossary(stdout,numerics,"  %-32s %s\n");

            printf("\nTargets: 1 sphere, 2 ellipsoid, 3 parallelepiped,"
                   "\n 4 random-sphere, 5 random-ellipsoid, 6 anthrax, 7 hedgehog,"
                   "\n 8 dead-anthrax, 9 random-parallelepiped, 10 nibble-ellipsoid,"
                   "\n 11 file, 12 bisphere, 13 hierarchical, 14 sphere-spec."
                   "\n Names from the --target column may be used instead of numbers."
                   "\n\nMaterials: --material <mre> <mi>, repeated once per material;"
                   "\n --seeds gives one random seed per material for the random targets."
                   "\n\nExample:"
                   "\n  %s --state out.dst --output out.dat --target sphere"
                   "\n      --material 1.5 0.01 --grid 64 --size-parameter 20\n",
                   prog);
        }
        else
            printf("\n%s %s (argtable3)\n",prog,DDA_VERSION);

        arg_freetable(argtable,sizeof(argtable)/sizeof(argtable[0]));
        return PARSE_EXIT_OK;
    }

    if(nerrors>0)
    {   arg_print_errors(stderr,o_end,prog);
        fprintf(stderr,"Try '%s --help' for the list of options.\n",prog);
        fprintf(stderr,"Note: list options take one value per option; repeat "
                       "the option or separate the values with commas.\n");

        arg_freetable(argtable,sizeof(argtable)/sizeof(argtable[0]));
        return PARSE_ERROR;
    }

    /* Map the parsed values into the option structure. */

    if(o_output->count>0)   opt->output_file=o_output->filename[0];
    if(o_state->count>0)    opt->state_file=o_state->filename[0];
    opt->restart=o_restart->count>0;

    if(o_samples->count>0)      opt->samples=o_samples->ival[0];
    if(o_tmin->count>0)         opt->theta_min=o_tmin->dval[0];
    if(o_tmax->count>0)         opt->theta_max=o_tmax->dval[0];
    if(o_tstep->count>0)        opt->theta_step=o_tstep->dval[0];
    if(o_psi->count>0)          opt->psi=o_psi->dval[0];
    if(o_iota->count>0)         opt->iota=o_iota->dval[0];
    if(o_kappa->count>0)        opt->kappa=o_kappa->dval[0];
    if(o_planes->count>0)       opt->planes=o_planes->ival[0];

    if(o_orientation->count>0)
    {   bool ok=false;
        int v=parse_code(orientation_names,o_orientation->sval[0],&ok);

        if(ok&&(v==0||v==1))
            opt->orientation=v;
        else
            errors.push_back(std::string("--orientation: expected fixed or mc, got \"")
                             +o_orientation->sval[0]+"\"");
    }

    if(o_target->count>0)
    {   bool ok=false;
        int v=parse_code(target_names,o_target->sval[0],&ok);

        if(ok&&v>=1&&v<=14)
            opt->target=v;
        else
            errors.push_back(std::string("--target: unknown target \"")
                             +o_target->sval[0]+"\"");
    }

    if(o_grid->count>0)
    {   std::vector<int> v;

        if(o_grid->count>1)
            errors.push_back("--grid must be given once");
        else if(!split_ints(o_grid->sval[0],&v))
            errors.push_back(std::string("--grid: expected integers, got \"")
                             +o_grid->sval[0]+"\"");
        else if(v.size()==1)
        {   opt->grid[0]=v[0];  opt->grid[1]=v[0];  opt->grid[2]=v[0];
        }
        else if(v.size()==3)
        {   opt->grid[0]=v[0];  opt->grid[1]=v[1];  opt->grid[2]=v[2];
        }
        else
            errors.push_back("--grid takes either 1 or 3 values");
    }

    if(o_size->count>0)         opt->size_parameter=o_size->dval[0];

    if(o_aspect->count>0)
    {   std::vector<double> v;

        if(o_aspect->count>1)
            errors.push_back("--aspect must be given once");
        else if(!split_values(o_aspect->sval[0],&v)||v.size()!=2)
            errors.push_back(std::string("--aspect: expected b/a,c/a, got \"")
                             +o_aspect->sval[0]+"\"");
        else
        {   opt->aspect_b=v[0];
            opt->aspect_c=v[1];
        }
    }

    if(o_rratio->count>0)       opt->radius_ratio=o_rratio->dval[0];
    if(o_levels->count>0)       opt->hierarchical_levels=o_levels->ival[0];
    if(o_bites->count>0)        opt->bites=o_bites->ival[0];
    if(o_depth->count>0)        opt->depth=o_depth->dval[0];
    if(o_packing->count>0)      opt->packing=o_packing->dval[0];
    if(o_mradius->count>0)      opt->monomer_radius=o_mradius->dval[0];
    if(o_mlr->count>0)          opt->monomer_length_ratio=o_mlr->dval[0];
    if(o_shape->count>0)        opt->shape_file=o_shape->filename[0];

    for(i=0;i<o_material->count;i++)
    {   std::vector<double> v;

        if(!split_values(o_material->sval[i],&v)||v.size()!=2)
        {   errors.push_back(std::string("--material: expected mre,mi, got \"")
                             +o_material->sval[i]+"\"");
            continue;
        }

        {   cplx m;
            m.r=(real)v[0];     m.i=(real)v[1];
            opt->materials.push_back(m);
        }
    }

    for(i=0;i<o_seeds->count;i++)
    {   std::vector<int> v;

        if(!split_ints(o_seeds->sval[i],&v))
        {   errors.push_back(std::string("--seeds: expected integers, got \"")
                             +o_seeds->sval[i]+"\"");
            continue;
        }

        opt->seeds.insert(opt->seeds.end(),v.begin(),v.end());
    }

    if(o_polarizability->count>0)
    {   bool ok=false;
        int v=parse_code(polarizability_names,o_polarizability->sval[0],&ok);

        if(ok&&v>=1&&v<=5)
            opt->polarizability=v;
        else
            errors.push_back(std::string("--polarizability: expected 1..5 or a name, got \"")
                             +o_polarizability->sval[0]+"\"");
    }

    if(o_tolerance->count>0)    opt->tolerance=o_tolerance->dval[0];
    if(o_steps->count>0)        opt->max_steps=o_steps->ival[0];

    if(o_fft->count>0)
    {   bool ok=false;
        int v=parse_code(fft_names,o_fft->sval[0],&ok);

        if(ok&&(v==0||v==1))
            opt->fft_mode=v;
        else
            errors.push_back(std::string("--fft: expected 0|1 or ct|rb, got \"")
                             +o_fft->sval[0]+"\"");
    }

    if(o_solver->count>0)
    {   bool ok=false;
        int v=parse_code(solver_names,o_solver->sval[0],&ok);

        if(ok&&v>=0&&v<=3)
            opt->solver=v;
        else
            errors.push_back(std::string("--solver: expected 0..3 or a name, got \"")
                             +o_solver->sval[0]+"\"");
    }

    if(o_tau->count>0)          opt->tau=o_tau->dval[0];
    if(o_threads->count>0)      opt->threads=o_threads->ival[0];

    if(o_seed->count>0)
    {   const char * s=o_seed->sval[0];
        char * end;
        unsigned long long v;

        errno=0;
        v=strtoull(s,&end,10);

        if(end==s||*end!=0||errno!=0)
            errors.push_back(std::string("--seed: expected an unsigned integer, got \"")
                             +s+"\"");
        else
        {   opt->has_seed=true;
            opt->seed=v;
        }
    }

    /* The argtable storage is not needed any more. Capture everything that is
       still taken from it before freeing the table. */

    bool interactive=(o_interactive->count>0||argc==1);

    arg_freetable(argtable,sizeof(argtable)/sizeof(argtable[0]));

    if(!errors.empty())
    {   fprintf(stderr,"error:\n%s\n",join_errors(errors).c_str());
        return PARSE_ERROR;
    }

    if(interactive)
        prompt_missing(opt);

    apply_defaults(opt);

    err=validate(*opt);
    if(!err.empty())
    {   fprintf(stderr,"error:\n%s\n",err.c_str());
        return PARSE_ERROR;
    }

    return PARSE_OK;
}

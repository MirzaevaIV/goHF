#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cint.h"
#include "libcint_wrap.h"

int wrap_cint2e_sph(double *buf, int *shls, int *atm, int natm, int *bas, int nbas, double *env){
    CINTOpt *opt = NULL;
    cint2e_sph_optimizer(&opt, atm, natm, bas, nbas, env);
    cint2e_sph(buf, shls, atm, natm, bas, nbas, env, opt);
    CINTdel_optimizer(&opt);
    return 0;
};

int wrap_cint1e_ovlp_sph(double *buf, int *shls, int *atm, int natm, int *bas, int nbas, double *env){
    cint1e_ovlp_sph(buf, shls, atm, natm, bas, nbas, env);
    return 0;
}




/* Data transform - Gaussian function transform */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/datastructures/array.h"
#include "ViennaRNA/math/functions.h"

#ifdef __GNUC__
# define INLINE inline
#else
# define INLINE
#endif


typedef struct {
  vrna_array(vrna_math_fun_dbl_f)           f;
  vrna_array(vrna_math_fun_dbl_opt_t)       f_opt;
  vrna_array(vrna_math_fun_dbl_opt_free_f)  f_opt_free;
} fun_chain_opt_t;


/*
 #################################
 # GLOBAL VARIABLES              #
 #################################
 */

/*
 #################################
 # PRIVATE VARIABLES             #
 #################################
 */

/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */
PRIVATE double
fun_chain_opt(double                  v,
              vrna_math_fun_dbl_opt_t options);


PRIVATE void
chain_option_free(vrna_math_fun_dbl_opt_t options);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC vrna_math_fun_dbl_f
vrna_math_fun_dbl_chain(vrna_math_fun_dbl_f           first_f,
                        vrna_math_fun_dbl_opt_t       first_f_opt,
                        vrna_math_fun_dbl_opt_free_f  first_f_opt_free,
                        vrna_math_fun_dbl_opt_t       *fun_options_p,
                        vrna_math_fun_dbl_opt_free_f  *fun_options_free)
{
  vrna_math_fun_dbl_f cb = NULL;

  if ((fun_options_p) &&
      (fun_options_free)) {
    fun_chain_opt_t *o = (fun_chain_opt_t *)vrna_alloc(sizeof(fun_chain_opt_t));

    vrna_array_init(o->f);
    vrna_array_init(o->f_opt);
    vrna_array_init(o->f_opt_free);

    if (first_f) {
      vrna_array_append(o->f, first_f);
      vrna_array_append(o->f_opt, first_f_opt);
      vrna_array_append(o->f_opt_free, first_f_opt_free);
    }

    cb                = fun_chain_opt;
    *fun_options_p    = (vrna_math_fun_dbl_opt_t)o;
    *fun_options_free = chain_option_free;
  }

  return cb;
}


PUBLIC size_t
vrna_math_fun_dbl_chain_append(vrna_math_fun_dbl_opt_t      chain,
                               vrna_math_fun_dbl_f          f,
                               vrna_math_fun_dbl_opt_t      f_opt,
                               vrna_math_fun_dbl_opt_free_f f_opt_free)
{
  if ((chain) &&
      (f)) {
    fun_chain_opt_t *o = (fun_chain_opt_t *)chain;

    vrna_array_append(o->f, f);
    vrna_array_append(o->f_opt, f_opt);
    vrna_array_append(o->f_opt_free, f_opt_free);

    return vrna_array_size(o->f);
  }

  return 0;
}


PUBLIC size_t
vrna_math_fun_dbl_chain_size(vrna_math_fun_dbl_opt_t chain)
{
  if (chain)
    return vrna_array_size(((fun_chain_opt_t *)chain)->f);

  return 0;
}


PUBLIC vrna_math_fun_dbl_f
vrna_math_fun_dbl_chain_at(vrna_math_fun_dbl_opt_t  chain,
                           size_t                   pos,
                           vrna_math_fun_dbl_opt_t  *fun_options_p)
{
  vrna_math_fun_dbl_f cb = NULL;

  if ((chain) &&
      (fun_options_p)) {
    fun_chain_opt_t *o = (fun_chain_opt_t *)chain;

    if (pos < vrna_array_size(o->f)) {
      cb              = o->f[pos];
      *fun_options_p  = o->f_opt[pos];
    } else {
      vrna_log_error("Requested position %ld exceeds total size (%ld) of function chain",
                     pos,
                     vrna_array_size(o->f));
    }
  }

  return cb;
}


/*
 #####################################
 # BEGIN OF STATIC HELPER FUNCTIONS  #
 #####################################
 */
PRIVATE void
chain_option_free(vrna_math_fun_dbl_opt_t options)
{
  if (options) {
    fun_chain_opt_t *o = (fun_chain_opt_t *)options;

    /* release memory for chained function objects */
    for (size_t i = 0; i < vrna_array_size(o->f); ++i)
      if (o->f_opt_free[i])
        o->f_opt_free[i](o->f_opt[i]);

    vrna_array_free(o->f);
    vrna_array_free(o->f_opt);
    vrna_array_free(o->f_opt_free);
  }

  free(options);
}


PRIVATE double
fun_chain_opt(double                  value,
              vrna_math_fun_dbl_opt_t options)
{
  double          result  = value;
  fun_chain_opt_t *o      = (fun_chain_opt_t *)options;

  for (size_t i = 0; i < vrna_array_size(o->f); ++i)
    result = o->f[i](result, o->f_opt[i]);

  return result;
}

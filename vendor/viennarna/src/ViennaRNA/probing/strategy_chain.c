/* SHAPE reactivity data handling */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/data/transform.h"
#include "ViennaRNA/probing/basic.h"

typedef struct {
  vrna_array(vrna_probing_strategy_f)     strategy;
  vrna_array(vrna_probing_strategy_opt_t) strategy_opt;
  vrna_array(vrna_auxdata_free_f)         strategy_opt_free;
} strategy_chain_t;


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


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC double *
vrna_probing_strategy_chain(vrna_fold_compound_t        *fc,
                            const double                *data,
                            size_t                      data_size,
                            unsigned int                target,
                            vrna_probing_strategy_opt_t options)
{
  unsigned char     converted;
  double            *pseudo_energies, *e;
  strategy_chain_t  *opt;

  pseudo_energies = NULL;

  if ((data) &&
      (data_size > 0) &&
      (options)) {
    opt             = (strategy_chain_t *)options;
    converted       = 0;
    pseudo_energies = (double *)vrna_alloc(sizeof(double) * data_size);
    pseudo_energies = memcpy(pseudo_energies, data, sizeof(double) * data_size);

    /* transform data into actual pseudo-energies */
    for (size_t s = 0; s < vrna_array_size(opt->strategy); ++s) {
      if (opt->strategy[s]) {
        e = opt->strategy[s](fc,
                             pseudo_energies,
                             data_size,
                             target,
                             opt->strategy_opt[s]);

        /* add data from this strategy */
        for (size_t i = 0; i < data_size; i++)
          pseudo_energies[i] += e[i];

        converted = 1;
      }
    }

    if (converted == 0) {
      free(pseudo_energies);
      pseudo_energies = NULL;
    }
  }

  return pseudo_energies;
}


PUBLIC vrna_probing_strategy_opt_t
vrna_probing_strategy_chain_data(vrna_probing_strategy_f      strategy_cb,
                                 vrna_probing_strategy_opt_t  strategy_cb_options,
                                 vrna_auxdata_free_f          strategy_cb_options_free)
{
  strategy_chain_t *d = (strategy_chain_t *)vrna_alloc(sizeof(strategy_chain_t));

  vrna_array_init(d->strategy);
  vrna_array_init(d->strategy_opt);
  vrna_array_init(d->strategy_opt_free);

  if (strategy_cb) {
    vrna_array_append(d->strategy, strategy_cb);
    vrna_array_append(d->strategy_opt, strategy_cb_options);
    vrna_array_append(d->strategy_opt_free, strategy_cb_options_free);
  }

  return (vrna_probing_strategy_opt_t)d;
}


PUBLIC size_t
vrna_probing_strategy_chain_data_append(vrna_probing_strategy_opt_t chain,
                                        vrna_probing_strategy_f     strategy_cb,
                                        vrna_probing_strategy_opt_t strategy_cb_options,
                                        vrna_auxdata_free_f         strategy_cb_options_free)
{
  if ((chain) &&
      (strategy_cb)) {
    strategy_chain_t *d = (strategy_chain_t *)chain;

    vrna_array_append(d->strategy, strategy_cb);
    vrna_array_append(d->strategy_opt, strategy_cb_options);
    vrna_array_append(d->strategy_opt_free, strategy_cb_options_free);

    return vrna_array_size(d->strategy);
  }

  return 0;
}


PUBLIC size_t
vrna_probing_strategy_chain_data_size(vrna_probing_strategy_opt_t chain)
{
  if (chain)
    return vrna_array_size(((strategy_chain_t *)chain)->strategy);

  return 0;
}


PUBLIC vrna_probing_strategy_f
vrna_probing_strategy_chain_data_at(vrna_probing_strategy_opt_t chain,
                                    size_t                      pos,
                                    vrna_probing_strategy_opt_t *strategy_cb_options)
{
  vrna_probing_strategy_f cb = NULL;

  if ((chain) &&
      (strategy_cb_options)) {
    strategy_chain_t *d = (strategy_chain_t *)chain;

    if (pos < vrna_array_size(d->strategy)) {
      cb                    = d->strategy[pos];
      *strategy_cb_options  = d->strategy_opt[pos];
    } else {
      vrna_log_error("Requested position %ld exceeds total size (%ld) of strategy chain",
                     pos,
                     vrna_array_size(d->strategy));
    }
  }

  return cb;
}


PUBLIC void
vrna_probing_strategy_chain_data_free(vrna_probing_strategy_opt_t options)
{
  strategy_chain_t *d = (strategy_chain_t *)options;

  /* release memory for chained strategy objects */
  for (size_t i = 0; i < vrna_array_size(d->strategy); ++i)
    if (d->strategy_opt_free[i])
      d->strategy_opt_free[i](d->strategy_opt[i]);

  vrna_array_free(d->strategy);
  vrna_array_free(d->strategy_opt);
  vrna_array_free(d->strategy_opt_free);

  free(options);
}


/*
 #####################################
 # BEGIN OF STATIC HELPER FUNCTIONS  #
 #####################################
 */

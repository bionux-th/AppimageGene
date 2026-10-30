#ifndef VIENNA_RNA_PACKAGE_PROBING_STRATEGY_NLOGP_H
#define VIENNA_RNA_PACKAGE_PROBING_STRATEGY_NLOGP_H


#include <ViennaRNA/probing/basic.h>
#include <ViennaRNA/math/functions.h>

/**
 *  @file     ViennaRNA/probing/strategy_nlogp.h
 *  @ingroup  probing_data_strategy
 *  @brief    This module provides the API for a generic strategy to convert structure probing data using its information content
 */

/**
 *  @addtogroup probing_data_strategy_nlogp
 *  @{
 */


/**
 *  @brief  Default options for the *nlogp* probing data conversion strategy
 *
 *  @see  vrna_probing_strategy_nlogp_options(), vrna_probing_data_nlogp(), vrna_probing_data_nlogp_trans(),
 *        vrna_probing_data_nlogp_comparative(), vrna_probing_data_nlogp_trans_comparative()
 */
#define VRNA_PROBING_STRATEGY_NLOGP_OPTIONS_DEFAULT              0

#define VRNA_PROBING_STRATEGY_NLOGP_NONSTD_BASE                 (1 << 0)
#define VRNA_PROBING_STRATEGY_NLOGP_FACTOR_CELCIUS              (1 << 1)
#define VRNA_PROBING_STRATEGY_NLOGP_FACTOR_KELVIN               (1 << 2)


double *
vrna_probing_strategy_nlogp(vrna_fold_compound_t        *fc,
                            const double                *data,
                            size_t                      data_size,
                            unsigned int                target,
                            vrna_probing_strategy_opt_t options);


vrna_probing_strategy_opt_t
vrna_probing_strategy_nlogp_options(double                        base,
                                    double                        factor,
                                    unsigned int                  options,
                                    vrna_math_fun_dbl_f           cb_preprocess,
                                    vrna_math_fun_dbl_opt_t       cb_preprocess_opt,
                                    vrna_math_fun_dbl_opt_free_f  cb_preprocess_opt_free);


void
vrna_probing_strategy_nlogp_options_free(vrna_probing_strategy_opt_t options);


vrna_probing_data_t
vrna_probing_data_nlogp(const double  *data,
                        unsigned int  n,
                        double        base,
                        double        factor,
                        unsigned int  options);


vrna_probing_data_t
vrna_probing_data_nlogp_trans(const double                  *data,
                              unsigned int                  n,
                              double                        base,
                              double                        factor,
                              unsigned int                  options,
                              vrna_math_fun_dbl_f           trans,
                              vrna_math_fun_dbl_opt_t       trans_options,
                              vrna_math_fun_dbl_opt_free_f  trans_options_free);


vrna_probing_data_t
vrna_probing_data_nlogp_comparative(const double        **datas,
                                    const unsigned int  *n,
                                    unsigned int        n_seq,
                                    double              base,
                                    double              factor,
                                    unsigned int        options,
                                    unsigned int        multi_params);


vrna_probing_data_t
vrna_probing_data_nlogp_trans_comparative(const double                  **datas,
                                          const unsigned int            *n,
                                          unsigned int                  n_seq,
                                          double                        base,
                                          double                        factor,
                                          unsigned int                  options,
                                          unsigned int                  multi_params,
                                          vrna_math_fun_dbl_f           *trans,
                                          vrna_math_fun_dbl_opt_t       *trans_options,
                                          vrna_math_fun_dbl_opt_free_f  *trans_options_free);


/**
 *  @}
 */

#endif

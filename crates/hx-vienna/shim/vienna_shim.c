/*
 * HelixBench 与内嵌 ViennaRNA 之间的最小 C 接口。
 *
 * 只暴露三个动作：折叠（MFE + 配分函数 + 配对概率）、naview 布局、释放。
 * 上游的 `vrna_md_t` 等结构体不跨 FFI：全部在这个文件里构造和释放。
 */

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "ViennaRNA/fold_compound.h"
#include "ViennaRNA/mfe.h"
#include "ViennaRNA/params/basic.h"
#include "ViennaRNA/part_func.h"
#include "ViennaRNA/probabilities/structures.h"
#include "ViennaRNA/plotting/layouts.h"
#include "ViennaRNA/subopt/wuchty.h"
#include "ViennaRNA/eval.h"

static void
md_init(vrna_md_t *md, double temperature, int no_gu_close, int no_lp)
{
  vrna_md_set_default(md);
  md->temperature   = temperature;
  md->noGUclosure   = no_gu_close ? 1 : 0;
  md->noLP         = no_lp ? 1 : 0;
}

/*
 * 折叠。`db` 由调用方提供（长度 >= 序列长度 + 1）。
 * 概率通过回调返回（只报 p > 5e-4 的配对，避免传一个 n² 矩阵）。
 *
 * 返回 0 成功，负数为错误码。
 */
int
hx_vienna_fold(const char         *seq,
               double             temperature,
               int                no_gu_close,
               int                no_lp,
               char               *db,
               int                db_len,
               double             *mfe_energy,
               double             *ensemble,
               double             *diversity,
               void               (*cb)(int i, int j, double p, void *ud),
               void               *ud)
{
  int len = (int)strlen(seq);

  if (len <= 0 || db_len < len + 1)
    return -1;

  vrna_md_t md;
  md_init(&md, temperature, no_gu_close, no_lp);

  vrna_fold_compound_t *fc = vrna_fold_compound(seq,
                                                &md,
                                                VRNA_OPTION_MFE | VRNA_OPTION_PF);
  if (!fc)
    return -2;

  db[0] = '\0';
  if (mfe_energy)
    *mfe_energy = (double)vrna_mfe(fc, db);

  if (ensemble)
    *ensemble = (double)vrna_pf(fc, NULL);

  if (diversity)
    *diversity = vrna_mean_bp_distance(fc);

  if (cb) {
    int           *idx  = fc->iindx;
    FLT_OR_DBL    *P    = fc->exp_matrices->probs;
    for (int i = 1; i < len; i++) {
      for (int j = i + 1; j <= len; j++) {
        double p = (double)P[idx[i] - j];
        if (p > 5e-4)
          cb(i - 1, j - 1, p, ud);
      }
    }
  }

  vrna_fold_compound_free(fc);
  return 0;
}

/*
 * 次优结构（vrna_subopt_cb）。`delta` 单位是 10 cal/mol（200 = 2.0 kcal/mol）。
 * 回调逐条返回 dot-bracket 与能量。
 */
typedef void (*hx_subopt_cb)(double energy, const char *structure, void *ud);

typedef struct {
  hx_subopt_cb cb;
  void        *ud;
} hx_subopt_ctx;

static void
subopt_trampoline(const char *structure, float energy, void *data)
{
  hx_subopt_ctx *c = (hx_subopt_ctx *)data;
  if (c->cb)
    c->cb((double)energy, structure, c->ud);
}

double
hx_vienna_eval(const char *seq,
               const char *db,
               double     temperature,
               int        no_gu_close,
               int        no_lp)
{
  if (!seq || !db)
    return 0.0;

  vrna_md_t md;
  md_init(&md, temperature, no_gu_close, no_lp);

  vrna_fold_compound_t *fc = vrna_fold_compound(seq, &md, VRNA_OPTION_DEFAULT);
  if (!fc)
    return 0.0;

  double e = vrna_eval_structure(fc, db);
  vrna_fold_compound_free(fc);
  return e;
}

int
hx_vienna_subopt(const char *seq,
                 double     temperature,
                 int        no_gu_close,
                 int        no_lp,
                 int        delta,
                 hx_subopt_cb cb,
                 void        *ud)
{
  if (!seq || !cb)
    return -1;

  vrna_md_t md;
  md_init(&md, temperature, no_gu_close, no_lp);
  /* Wuchty 需要唯一的 multibranch 回溯矩阵 */
  md.uniq_ML = 1;

  vrna_fold_compound_t *fc = vrna_fold_compound(seq, &md, VRNA_OPTION_DEFAULT);
  if (!fc)
    return -2;

  hx_subopt_ctx ctx = { cb, ud };
  vrna_subopt_cb(fc, delta, subopt_trampoline, &ctx);
  vrna_fold_compound_free(fc);
  return 0;
}

/*
 * RNAplot 的 naview 布局。x/y 由调用方分配（长度 >= 结构长度）。
 * 返回碱基个数，负数为错误码。
 */
int
hx_vienna_naview(const char *structure,
                 double     *x,
                 double     *y,
                 int        max_len)
{
  if (!structure)
    return -1;

  int len = (int)strlen(structure);
  if (len > max_len)
    return -2;

  vrna_plot_layout_t *pl = vrna_plot_layout(structure, VRNA_PLOT_TYPE_NAVIEW);
  if (!pl)
    return -3;

  for (int i = 0; i < len; i++) {
    x[i] = (double)pl->x[i];
    y[i] = (double)pl->y[i];
  }

  vrna_plot_layout_free(pl);
  return len;
}

# HelixBench 对 ViennaRNA 源码做的改动

上游：**ViennaRNA 2.7.2**（PyPI sdist `viennarna-2.7.2.tar.gz`，含完整 C 源码）。
许可证见同目录 `COPYING`（研究/教育/商业使用免费，但**再分发不得收费**，且必须署名；
若要把它并入商业产品，需联系维也纳大学的作者）。本仓库只收录源码并在构建时静态编译。

## 保留了什么

- `src/ViennaRNA/**` 的全部 `.c` / `.h` / `.inc`（`zscore/` 例外，见下）。
- `src/cephes/`（`kn.c` `const.c` `mtherr.c` `expn.c` 及其头文件）、`src/json/`。
- `src/ViennaRNA/static/` 只保留绘图模板（`templates_svg.h`、`templates_postscript.h`、
  `svg/`、`postscript/`）与 SHAPE 先验（`probing_data_priors.h` + `misc/prior_probing_*.hex`）。

## 删掉了什么

| 删除项 | 原因 |
|---|---|
| `src/ViennaRNA/zscore/` | 依赖 dlib（C++），HelixBench 不用 z-score |
| `src/ViennaRNA/utils/svm_utils.c` | 依赖 libsvm，未使用 |
| `src/ViennaRNA/utils/higher_order_functions_{sse41,avx512}.c` | SIMD 特化，未启用 |
| `src/ViennaRNA/constraints/sc_cb_mod_wrappers.c` | 只包装修饰碱基参数集（见下），无人引用 |
| `src/ViennaRNA/static/misc/*.hex`（除 4 个 `prior_probing_*`） | 内置参数集转储，约 12 MB；默认 Turner 2004 参数来自 `params/default.c` 的编译表 |
| `src/ViennaRNA/static/energy_parameter_sets.h` | 同上（只被 `params/io.c` 与已删除的 `sc_cb_mod_wrappers.c` 引用） |
| 文档、测试数据、Python/SWIG 接口、`src/dlib-*`、`src/libsvm-*` | 构建不需要 |

## 唯一的源码改动

`src/ViennaRNA/params/io.c`：7 个「从内置字符串加载参数集」的函数
（`vrna_params_load_RNA_Turner2004` / `Turner1999` / `Andronescu2007` / `Langdon2018` /
`misc_special_hairpins` / `DNA_Mathews2004` / `DNA_Mathews1999`）里，
把读取被删掉的 `parameter_set_*` 数组的调用换成 `return -1;`，并加了 `/* HelixBench: ... */` 注释。

**这不影响折叠**：ViennaRNA 的默认参数来自 `params/default.c` 里编译好的 Turner 2004 表，
以上函数只服务于「运行时显式切换参数集」的旧 API，HelixBench 不调用它们。
（`vrna_params_load_*` 若在将来需要，把对应的 `.hex` 文件放回 `static/misc/` 即可。）

## 可移植性

- **不使用**上游的 `config.h`（autoconf 在构建机生成的 `HAVE_*` 会绑定平台）。
  `crates/hx-vienna/build.rs` 里放了一份只含版本号的 `shim/config.h`；不定义
  `HAVE_CONFIG_H` 时源码走自带的可移植回退分支（`strdup`、`erand48` 等都有 fallback）。
- 上游源码本身对 Windows 有 `_WIN32` 分支（`ViennaRNA/intern/unistd_win.h`），
  macOS/Linux 用 POSIX 分支。
- `-DVRNA_WITH_NAVIEW_LAYOUT=1` 打开 naview（上游 `config.h` 里默认是注释掉的，
  不打开时 `vrna_plot_layout(..., VRNA_PLOT_TYPE_NAVIEW)` 会返回无效布局）。

## 升级步骤

1. 用新版本的 sdist 覆盖 `src/`（只保留上表列出的部分）。
2. 重新应用 `params/io.c` 的补丁（搜索 `parameter_set_`）。
3. 更新 `crates/hx-vienna/build.rs` 与 `crates/hx-vienna/src/lib.rs` 里的版本号。
4. `cargo test -p hx-vienna`（tRNA 的 MFE / 点括号 / 系综自由能 / 多样性是回归基准）。

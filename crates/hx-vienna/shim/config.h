/*
 * HelixBench 为内嵌 ViennaRNA 提供的可移植 config.h。
 *
 * 上游的 config.h 是 autoconf 在某个具体平台上生成的，含大量 HAVE_* 探测结果
 * （HAVE_ASPRINTF、HAVE_ERAND48、HAVE_UNISTD_H ……），直接拿来会在
 * Windows / macOS 上编错。这里只保留版本号，把所有 HAVE_* 都留空，
 * 让 ViennaRNA 源码走它自带的可移植回退实现。
 *
 * 需要打开的功能（naview 布局、版本字符串）由 build.rs 的 -D 传入。
 */
#ifndef HELIXBENCH_VIENNA_CONFIG_H
#define HELIXBENCH_VIENNA_CONFIG_H

#define VRNA_VERSION "2.7.2"
#define VRNA_VERSION_MAJOR 2
#define VRNA_VERSION_MINOR 7
#define VRNA_VERSION_PATCH 2

/* 数学库在 Linux/macOS 需要 -lm（build.rs 已经处理），Windows 不需要。 */
#define HAVE_MATH_H 1
#define HAVE_FLOAT_H 1
#define HAVE_LIMITS_H 1
#define HAVE_STDINT_H 1
#define HAVE_INTTYPES_H 1
#define HAVE_STDDEF_H 1
#define HAVE_STDLIB_H 1
#define HAVE_STRING_H 1
#define HAVE_CTYPE_H 1

#endif

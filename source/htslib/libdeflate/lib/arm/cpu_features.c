/*
 * arm/cpu_features.c - feature detection for ARM CPUs
 *
 * Copyright 2018 Eric Biggers
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

/*
 * ARM CPUs don't have a standard way for unprivileged programs to detect CPU
 * features.  But an OS-specific way can be used when available.
 */

#ifdef __APPLE__
#  undef _ANSI_SOURCE
#  undef _DARWIN_C_SOURCE
#  define _DARWIN_C_SOURCE /* for sysctlbyname() */
#endif

#include "../cpu_features_common.h" /* must be included first */
#include "cpu_features.h"

#ifdef ARM_CPU_FEATURES_KNOWN
/* Runtime ARM CPU feature detection is supported. */

#ifdef __linux__
/*
 * On Linux, arm32 and arm64 CPU features can be detected by reading the
 * AT_HWCAP and AT_HWCAP2 values from /proc/self/auxv.
 *
 * Ideally we'd use the C library function getauxval(), but it's not guaranteed
 * to be available: it was only added to glibc in 2.16, and in Android it was
 * added to API level 18 for arm32 and level 21 for arm64.
 */

#include <errno.h>
#include <fcntl.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#define AT_HWCAP	16
#define AT_HWCAP2	26

static void scan_auxv(unsigned long *hwcap, unsigned long *hwcap2)
{
	int fd;
	unsigned long auxbuf[32];
	int filled = 0;
	int i;

	fd = open("/proc/self/auxv", O_RDONLY | O_CLOEXEC);
	if (fd < 0)
		return;

	for (;;) {
		do {
			int ret = read(fd, &((char *)auxbuf)[filled],
				       sizeof(auxbuf) - filled);
			if (ret <= 0) {
				if (ret < 0 && errno == EINTR)
					continue;
				goto out;
			}
			filled += ret;
		} while (filled < 2 * sizeof(long));

		i = 0;
		do {
			unsigned long type = auxbuf[i];
			unsigned long value = auxbuf[i + 1];

			if (type == AT_HWCAP)
				*hwcap = value;
			else if (type == AT_HWCAP2)
				*hwcap2 = value;
			i += 2;
			filled -= 2 * sizeof(long);
		} while (filled >= 2 * sizeof(long));

		memmove(auxbuf, &auxbuf[i], filled);
	}
out:
	close(fd);
}

static u32 query_arm_cpu_features(void)
{
	u32 features = 0;
	unsigned long hwcap = 0;
	unsigned long hwcap2 = 0;

	scan_auxv(&hwcap, &hwcap2);

#ifdef ARCH_ARM32
	STATIC_ASSERT(sizeof(long) == 4);
	if (hwcap & (1 << 12))	/* HWCAP_NEON */
		features |= ARM_CPU_FEATURE_NEON;
#else
	STATIC_ASSERT(sizeof(long) == 8);
	if (hwcap & (1 << 1))	/* HWCAP_ASIMD */
		features |= ARM_CPU_FEATURE_NEON;
	if (hwcap & (1 << 4))	/* HWCAP_PMULL */
		features |= ARM_CPU_FEATURE_PMULL;
	if (hwcap & (1 << 7))	/* HWCAP_CRC32 */
		features |= ARM_CPU_FEATURE_CRC32;
	if (hwcap & (1 << 17))	/* HWCAP_SHA3 */
		features |= ARM_CPU_FEATURE_SHA3;
	if (hwcap & (1 << 20))	/* HWCAP_ASIMDDP */
		features |= ARM_CPU_FEATURE_DOTPROD;
#endif
	return features;
}

#ifdef ARCH_ARM64
/*
 * Return whether cpu0's MIDR_EL1 identifies one of the Arm Neoverse
 * V-class server cores (V1 / V2 / V3 / V3AE).  MIDR_EL1 is exposed
 * unprivileged via sysfs (added in Linux 4.7).  Reading cpu0 only is fine
 * in practice: no Neoverse V-class server SKU has shipped as part of a
 * big.LITTLE cluster.  Any failure (file missing, read error, parse
 * failure, unrecognized CPU) returns false.
 */
static bool arm64_cpu_is_neoverse_v_class(void)
{
	int fd;
	char buf[32];
	ssize_t n;
	unsigned long midr;
	u32 part;

	fd = open("/sys/devices/system/cpu/cpu0/regs/identification/midr_el1",
		  O_RDONLY | O_CLOEXEC);
	if (fd < 0)
		return false;
	do {
		n = read(fd, buf, sizeof(buf) - 1);
	} while (n < 0 && errno == EINTR);
	close(fd);
	if (n <= 0)
		return false;
	buf[n] = '\0';
	midr = strtoul(buf, NULL, 0); /* sysfs prints "0x%016llx\n" */

	/* MIDR_EL1: [31:24]=Implementer, [15:4]=PartNum. */
	if (((midr >> 24) & 0xff) != 0x41) /* Implementer must be Arm Ltd. */
		return false;
	part = (midr >> 4) & 0xfff;
	switch (part) {
	case 0xd40: /* Neoverse V1   (e.g. AWS Graviton 3) */
	case 0xd4f: /* Neoverse V2   (e.g. AWS Graviton 4) */
	case 0xd83: /* Neoverse V3AE */
	case 0xd84: /* Neoverse V3 */
		return true;
	}
	return false;
}
#endif /* ARCH_ARM64 */

#elif defined(__APPLE__)
/* On Apple platforms, arm64 CPU features can be detected via sysctlbyname(). */

#include <sys/types.h>
#include <sys/sysctl.h>
#include <TargetConditionals.h>

static const struct {
	const char *name;
	u32 feature;
} feature_sysctls[] = {
	{ "hw.optional.neon",		  ARM_CPU_FEATURE_NEON },
	{ "hw.optional.AdvSIMD",	  ARM_CPU_FEATURE_NEON },
	{ "hw.optional.arm.FEAT_PMULL",	  ARM_CPU_FEATURE_PMULL },
	{ "hw.optional.armv8_crc32",	  ARM_CPU_FEATURE_CRC32 },
	{ "hw.optional.armv8_2_sha3",	  ARM_CPU_FEATURE_SHA3 },
	{ "hw.optional.arm.FEAT_SHA3",	  ARM_CPU_FEATURE_SHA3 },
	{ "hw.optional.arm.FEAT_DotProd", ARM_CPU_FEATURE_DOTPROD },
};

static u32 query_arm_cpu_features(void)
{
	u32 features = 0;
	size_t i;

	for (i = 0; i < ARRAY_LEN(feature_sysctls); i++) {
		const char *name = feature_sysctls[i].name;
		u32 val = 0;
		size_t valsize = sizeof(val);

		if (sysctlbyname(name, &val, &valsize, NULL, 0) == 0 &&
		    valsize == sizeof(val) && val == 1)
			features |= feature_sysctls[i].feature;
	}
	return features;
}
#elif defined(_WIN32)

#include <windows.h>

#ifndef PF_ARM_V82_DP_INSTRUCTIONS_AVAILABLE /* added in Windows SDK 20348 */
#  define PF_ARM_V82_DP_INSTRUCTIONS_AVAILABLE 43
#endif

static u32 query_arm_cpu_features(void)
{
	u32 features = ARM_CPU_FEATURE_NEON;

	if (IsProcessorFeaturePresent(PF_ARM_V8_CRYPTO_INSTRUCTIONS_AVAILABLE))
		features |= ARM_CPU_FEATURE_PMULL;
	if (IsProcessorFeaturePresent(PF_ARM_V8_CRC32_INSTRUCTIONS_AVAILABLE))
		features |= ARM_CPU_FEATURE_CRC32;
	if (IsProcessorFeaturePresent(PF_ARM_V82_DP_INSTRUCTIONS_AVAILABLE))
		features |= ARM_CPU_FEATURE_DOTPROD;

	/* FIXME: detect SHA3 support too. */

	return features;
}
#else
#error "unhandled case"
#endif

static const struct cpu_feature arm_cpu_feature_table[] = {
	{ARM_CPU_FEATURE_NEON,		"neon"},
	{ARM_CPU_FEATURE_PMULL,		"pmull"},
	{ARM_CPU_FEATURE_PREFER_PMULL,  "prefer_pmull"},
	{ARM_CPU_FEATURE_CRC32,		"crc32"},
	{ARM_CPU_FEATURE_SHA3,		"sha3"},
	{ARM_CPU_FEATURE_DOTPROD,	"dotprod"},
};

/*
 * Whether to set ARM_CPU_FEATURE_PREFER_PMULL on this CPU.  This is the
 * right choice on CPUs whose pmull pipes have more aggregate throughput
 * than the crc32 unit -- in measured cases by a wide margin: the Apple M
 * series sustains ~68 GB/s on pmull vs ~25 GB/s on crc32 (M1), and the Arm
 * Neoverse V class sustains ~40 GB/s vs ~22 GB/s (Graviton 4 / V2).
 *
 * We detect Apple at compile time, and Neoverse V-class cores at runtime
 * via MIDR_EL1 on Linux.  Elsewhere we leave this unset, and the dispatcher
 * picks the crc32-instruction path which is the right default for most
 * other modern ARM CPUs.
 */
static bool arm_cpu_prefers_pmull(void)
{
#if defined(__APPLE__) && TARGET_OS_OSX
	return true;
#elif defined(__linux__) && defined(ARCH_ARM64)
	if (arm64_cpu_is_neoverse_v_class())
		return true;
#endif
#ifdef TEST_SUPPORT__DO_NOT_USE
	return true;
#endif
	return false;
}

volatile u32 libdeflate_arm_cpu_features = 0;

void libdeflate_init_arm_cpu_features(void)
{
	u32 features = query_arm_cpu_features();

	if (arm_cpu_prefers_pmull())
		features |= ARM_CPU_FEATURE_PREFER_PMULL;

	disable_cpu_features_for_testing(&features, arm_cpu_feature_table,
					 ARRAY_LEN(arm_cpu_feature_table));

	libdeflate_arm_cpu_features = features | ARM_CPU_FEATURES_KNOWN;
}

#endif /* ARM_CPU_FEATURES_KNOWN */

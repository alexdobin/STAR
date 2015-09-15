#pragma once

#ifndef PRNG48_H_
#define PRNG48_H_

#ifdef __cplusplus
extern "C" {
#endif

	// [0.0, 1.0)
	double drand48(void);

	// 0 and 2^31
	long int lrand48(void);

	// initialization
	void srand48(long int seedval);

#ifdef __cplusplus
}
#endif

#endif

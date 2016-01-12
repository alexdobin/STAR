#include "PRNG48.h"
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>

boost::mt19937 gen;
boost::uniform_real<> uni_dist(0,1);
boost::variate_generator<boost::mt19937&, boost::uniform_real<> > drand_gen(gen, uni_dist);

// [0.0, 1.0)
double drand48(void) {
	return drand_gen();
}

// 0 and 2^31
long int lrand48(void) {
	return gen();
}

// initialization
void srand48(long int seedval) {
	gen.seed((uint32_t)seedval);
}

#ifndef BAM_ENDIAN_H
#define BAM_ENDIAN_H

#include <stdint.h>

static inline int bam_is_big_endian()
{
	long one= 1;
	return !(*((char *)(&one)));
}
static inline uint16_t bam_swap_endian_2(uint16_t v)
{
	return (uint16_t)(((v & 0x00FF00FFU) << 8) | ((v & 0xFF00FF00U) >> 8));
}
static inline void *bam_swap_endian_2p(void *x)
{
	*(uint16_t*)x = bam_swap_endian_2(*(uint16_t*)x);
	return x;
}
static inline uint32_t bam_swap_endian_4(uint32_t v)
{
	v = ((v & 0x0000FFFFU) << 16) | (v >> 16);
	return ((v & 0x00FF00FFU) << 8) | ((v & 0xFF00FF00U) >> 8);
}
static inline void *bam_swap_endian_4p(void *x)
{
	*(uint32_t*)x = bam_swap_endian_4(*(uint32_t*)x);
	return x;
}
static inline uint64_t bam_swap_endian_8(uint64_t v)
{
	v = ((v & 0x00000000FFFFFFFFLLU) << 32) | (v >> 32);
	v = ((v & 0x0000FFFF0000FFFFLLU) << 16) | ((v & 0xFFFF0000FFFF0000LLU) >> 16);
	return ((v & 0x00FF00FF00FF00FFLLU) << 8) | ((v & 0xFF00FF00FF00FF00LLU) >> 8);
}
static inline void *bam_swap_endian_8p(void *x)
{
	*(uint64_t*)x = bam_swap_endian_8(*(uint64_t*)x);
	return x;
}

#endif

/*
    Time will show if this module will be merged into others 
    or perhaps removed completely.
*/
#ifndef VCF_UTILS_H
#define VCF_UTILS_H

#include "vcf.h"


/**
 *  bcf_trim_alleles() - remove ALT alleles unused in genotype fields
 *  @header:  for access to BCF_DT_ID dictionary
 *  @line:    VCF line obtain from vcf_parse1
 *
 *  Returns the number of removed alleles.
 *
 *  todo: BCF output
 */
int bcf_trim_alleles(const bcf_hdr_t *header, bcf1_t *line);


/**
 *  bcf_remove_alleles() - remove ALT alleles according to bitmask @mask
 *  @header:  for access to BCF_DT_ID dictionary
 *  @line:    VCF line obtained from vcf_parse1
 *  @mask:    alleles to remove
 *
 *  todo: BCF output
 */
void bcf_remove_alleles(const bcf_hdr_t *header, bcf1_t *line, int mask);


/**
 *  bcf_calc_ac() - calculate the number of REF and ALT alleles
 *  @header:  for access to BCF_DT_ID dictionary
 *  @line:    VCF line obtained from vcf_parse1
 *  @ac:      array of length line->n_allele
 *  @which:   determine if INFO/AN,AC and indv fields be used
 *
 *  Returns 1 if the call succeeded, or 0 if the value could not 
 *  be determined.
 *
 *  The value of @which determines if existing INFO/AC,AN can be 
 *  used (BCF_UN_INFO) and and if indv fields can be splitted 
 *  (BCF_UN_FMT). 
 */
int bcf_calc_ac(const bcf_hdr_t *header, bcf1_t *line, int *ac, int which);


/**
 * bcf_gt_type() - determines type of the genotype
 * @fmt_ptr:  the GT format field as set for example by set_fmt_ptr
 * @isample:  sample index (starting from 0)
 * @ial:      index of the 1st non-reference allele (starting from 1)
 * @jal:      index of the 2nd non-reference allele (starting from 1)
 *
 * Returns the type of the genotype (one of GT_HOM_RR, GT_HET_RA,
 * GT_HOM_AA, GT_HET_AA, GT_HAPL_R, GT_HAPL_A or GT_UNKN). If $ial 
 * is not NULL and the genotype has one or more non-reference 
 * alleles, $ial will be set. In case of GT_HET_AA, $ial is the 
 * position of the allele which appeared first in ALT. If $jal is 
 * not null and the genotype is GT_HET_AA, $jal will be set and is 
 * the position of the second allele in ALT.
 */
#define GT_HOM_RR 0 // note: the actual value of GT_* matters, used in dosage r2 calculation
#define GT_HOM_AA 1
#define GT_HET_RA 2
#define GT_HET_AA 3
#define GT_HAPL_R 4
#define GT_HAPL_A 5
#define GT_UNKN   6
int bcf_gt_type(bcf_fmt_t *fmt_ptr, int isample, int *ial, int *jal);

static inline int bcf_acgt2int(char c)
{
    if ( (int)c>96 ) c -= 32;
    if ( c=='A' ) return 0;
    if ( c=='C' ) return 1;
    if ( c=='G' ) return 2;
    if ( c=='T' ) return 3;
    return -1;
}
#define bcf_int2acgt(i) "ACGT"[i]

/**
  * bcf_ij2G() - common task: allele indexes to Number=G index (diploid)
  * @i,j:  allele indexes, 0-based, i<=j
  * 
  * Returns index to the Number=G diploid array
  */
#define bcf_ij2G(i, j) ((j)*((j)+1)/2+(i))

#endif

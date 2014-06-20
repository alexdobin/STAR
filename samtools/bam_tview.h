#ifndef BAM_TVIEW_H
#define BAM_TVIEW_H

#include <ctype.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <stdarg.h>
#include "bam.h"
#include "faidx.h"
#include "bam2bcf.h"
#include "sam_header.h"
#include "khash.h"

KHASH_MAP_INIT_STR(kh_rg, const char *)

typedef struct AbstractTview {
	int mrow, mcol;
	
	bam_index_t *idx;
	bam_lplbuf_t *lplbuf;
	bam_header_t *header;
	bamFile fp;
	int curr_tid, left_pos;
	faidx_t *fai;
	bcf_callaux_t *bca;

	int ccol, last_pos, row_shift, base_for, color_for, is_dot, l_ref, ins, no_skip, show_name;
	char *ref;
    khash_t(kh_rg) *rg_hash;
    /* callbacks */
    void (*my_destroy)(struct AbstractTview* );
    void (*my_mvprintw)(struct AbstractTview* ,int,int,const char*,...);
    void (*my_mvaddch)(struct AbstractTview*,int,int,int);
    void (*my_attron)(struct AbstractTview*,int);
    void (*my_attroff)(struct AbstractTview*,int);
    void (*my_clear)(struct AbstractTview*);
    int (*my_colorpair)(struct AbstractTview*,int);
    int (*my_drawaln)(struct AbstractTview*,int,int);
    int (*my_loop)(struct AbstractTview*);
    int (*my_underline)(struct AbstractTview*);
} tview_t;


char bam_aux_getCEi(bam1_t *b, int i);
char bam_aux_getCSi(bam1_t *b, int i);
char bam_aux_getCQi(bam1_t *b, int i);

#define TV_MIN_ALNROW 2
#define TV_MAX_GOTO  40
#define TV_LOW_MAPQ  10

#define TV_COLOR_MAPQ   0
#define TV_COLOR_BASEQ  1
#define TV_COLOR_NUCL   2
#define TV_COLOR_COL    3
#define TV_COLOR_COLQ   4

#define TV_BASE_NUCL 0
#define TV_BASE_COLOR_SPACE 1

int tv_pl_func(uint32_t tid, uint32_t pos, int n, const bam_pileup1_t *pl, void *data);
int base_tv_init(tview_t*,const char *fn, const char *fn_fa, const char *samples);
void base_tv_destroy(tview_t*);
int base_draw_aln(tview_t *tv, int tid, int pos);

typedef struct Tixel
	{
	int ch;
	int attributes;
	}tixel_t;

#endif


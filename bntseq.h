/* The MIT License

   Copyright (c) 2008 Genome Research Ltd (GRL).

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

/* Contact: Heng Li <lh3@sanger.ac.uk> */

#ifndef BWT_BNTSEQ_H
#define BWT_BNTSEQ_H

#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <zlib.h>

#ifndef BWA_UBYTE
#define BWA_UBYTE
typedef uint8_t ubyte_t;
#endif

typedef struct {
	int64_t offset; // 序列在pac上的偏移。例如，两条序列分别长100bp 200bp，那么它们的偏移分别是0bp 100bp
	int32_t len; // 序列长度，即bp数
	int32_t n_ambs; // 该序列的holes的数目
	uint32_t gi;
	int32_t is_alt;
	char *name, *anno; // 序列的名字，和名字的注释。
} bntann1_t; // 序列的注释信息

typedef struct {
	int64_t offset; // hole在pac上的偏移
	int32_t len; // hole的长度，即N的个数
	char amb; // 填充hole的字符，也就是N
} bntamb1_t; // hole的信息

typedef struct {
	int64_t l_pac; // pac序列的bp数（值为fasta文件中所有序列的bp数相加，如果pac序列也包含了反向互补序列，那么还要乘以二）
	int32_t n_seqs; // fasta文件中序列的条数
	uint32_t seed;	// 定值，为11
	bntann1_t *anns; // n_seqs elements，保存序列的注释信息
	int32_t n_holes;
	bntamb1_t *ambs; // n_holes elements，保存序列的holes信息，hole指的是序列中为N的区域。例如ACTNGGCANNNTCGACGT中就有两个hole。这里保存的是所有序列的holes，不同的序列的holes的边界保存在anns中。
	FILE *fp_pac;
} bntseq_t;

extern unsigned char nst_nt4_table[256];

#ifdef __cplusplus
extern "C" {
#endif

	void bns_dump(const bntseq_t *bns, const char *prefix);
	bntseq_t *bns_restore(const char *prefix);
	bntseq_t *bns_restore_core(const char *ann_filename, const char* amb_filename, const char* pac_filename);
	void bns_destroy(bntseq_t *bns);
	int64_t bns_fasta2bntseq(gzFile fp_fa, const char *prefix, int for_only);
	int bns_pos2rid(const bntseq_t *bns, int64_t pos_f);
	int bns_cnt_ambi(const bntseq_t *bns, int64_t pos_f, int len, int *ref_id);
	uint8_t *bns_get_seq(int64_t l_pac, const uint8_t *pac, int64_t beg, int64_t end, int64_t *len);
	uint8_t *bns_fetch_seq(const bntseq_t *bns, const uint8_t *pac, int64_t *beg, int64_t mid, int64_t *end, int *rid);
	int bns_intv2rid(const bntseq_t *bns, int64_t rb, int64_t re);

#ifdef __cplusplus
}
#endif

static inline int64_t bns_depos(const bntseq_t *bns, int64_t pos, int *is_rev)
{
	return (*is_rev = (pos >= bns->l_pac))? (bns->l_pac<<1) - 1 - pos : pos;
}

#endif

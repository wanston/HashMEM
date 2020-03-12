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

#ifndef BWA_BWT_H
#define BWA_BWT_H

#include <stdint.h>
#include <stddef.h>

// requirement: (OCC_INTERVAL%16 == 0); please DO NOT change this line because some part of the code assume OCC_INTERVAL=0x80
#define OCC_INTV_SHIFT 7
#define OCC_INTERVAL   (1LL<<OCC_INTV_SHIFT)
#define OCC_INTV_MASK  (OCC_INTERVAL - 1)

#ifndef BWA_UBYTE
#define BWA_UBYTE
typedef unsigned char ubyte_t;
#endif

typedef uint64_t bwtint_t;

typedef struct {
	bwtint_t primary; // S^{-1}(0), or the primary index of BWT，$在BWT中应该插入的位置
	bwtint_t L2[5]; // C(), cumulative count，C[0]到C[4]分别表示小于ACGTN的字符数
	bwtint_t seq_len; // sequence length, fasta中所有序列的长度相加得到
	bwtint_t bwt_size; // size of bwt,表示bwt指针指向的buffer含有多少uint32_t元素, 值是 seq_len/16向上取整
	uint32_t *bwt; // BWT，不含$字符，其中每个碱基用2bit表示
	// occurance array, separated to two parts
	uint32_t cnt_table[256];
	// suffix array
	int sa_intv; // 默认是32
	bwtint_t n_sa; // 后缀数组sa的大小，值为 (seq_len+1+sa_intv-1) / sa_intv
	bwtint_t *sa; // 后缀数组
} bwt_t;

typedef struct {
	bwtint_t x[3];
} bwtintv_x_t;

typedef struct {
	bwtint_t x[3], info; // x[0]表示low边界，x[1]表示反向互补序列的low边界，x[2]表示interval的大小，info高32bit和低32bit分别表示此interval对应的起止位置
} bwtintv_t; // 用于表示bwt的interval

typedef struct { size_t n, m; bwtintv_t *a; } bwtintv_v; // vector结构，n表示size，m表示capacity，a则是指向数组

typedef struct {
	bwtintv_t e1;
	bwtintv_t e2;
} bwtintv_t2;

typedef struct { size_t n, m; bwtintv_t2 *a; } bwtintv_v2;

/* For general OCC_INTERVAL, the following is correct:
#define bwt_bwt(b, k) ((b)->bwt[(k)/OCC_INTERVAL * (OCC_INTERVAL/(sizeof(uint32_t)*8/2) + sizeof(bwtint_t)/4*4) + sizeof(bwtint_t)/4*4 + (k)%OCC_INTERVAL/16])
#define bwt_occ_intv(b, k) ((b)->bwt + (k)/OCC_INTERVAL * (OCC_INTERVAL/(sizeof(uint32_t)*8/2) + sizeof(bwtint_t)/4*4)
*/

// The following two lines are ONLY correct when OCC_INTERVAL==0x80
#define bwt_bwt(b, k) ((b)->bwt[((k)>>7<<4) + sizeof(bwtint_t) + (((k)&0x7f)>>4)])
#define bwt_occ_intv(b, k) ((b)->bwt + ((k)>>7<<4))

/* retrieve a character from the $-removed BWT string. Note that
 * bwt_t::bwt is not exactly the BWT string and therefore this macro is
 * called bwt_B0 instead of bwt_B */
#define bwt_B0(b, k) (bwt_bwt(b, k)>>((~(k)&0xf)<<1)&3)

#define bwt_set_intv(bwt, c, ik) ((ik).x[0] = (bwt)->L2[(int)(c)]+1, (ik).x[2] = (bwt)->L2[(int)(c)+1]-(bwt)->L2[(int)(c)], (ik).x[1] = (bwt)->L2[3-(c)]+1, (ik).info = 0)

#ifdef __cplusplus
extern "C" {
#endif

	void bwt_dump_bwt(const char *fn, const bwt_t *bwt);
	void bwt_dump_sa(const char *fn, const bwt_t *bwt);

	bwt_t *bwt_restore_bwt(const char *fn);
	void bwt_restore_sa(const char *fn, bwt_t *bwt);

	void bwt_destroy(bwt_t *bwt);

	void bwt_bwtgen(const char *fn_pac, const char *fn_bwt); // from BWT-SW
	void bwt_bwtgen2(const char *fn_pac, const char *fn_bwt, int block_size); // from BWT-SW
	void bwt_cal_sa(bwt_t *bwt, int intv);

	void bwt_bwtupdate_core(bwt_t *bwt);

	bwtint_t bwt_occ(const bwt_t *bwt, bwtint_t k, ubyte_t c);
	void bwt_occ4(const bwt_t *bwt, bwtint_t k, bwtint_t cnt[4]);
	bwtint_t bwt_sa(const bwt_t *bwt, bwtint_t k);

	// more efficient version of bwt_occ/bwt_occ4 for retrieving two close Occ values
	void bwt_gen_cnt_table(bwt_t *bwt);
	void bwt_2occ(const bwt_t *bwt, bwtint_t k, bwtint_t l, ubyte_t c, bwtint_t *ok, bwtint_t *ol);
	void bwt_2occ4(const bwt_t *bwt, bwtint_t k, bwtint_t l, bwtint_t cntk[4], bwtint_t cntl[4]);

	int bwt_match_exact(const bwt_t *bwt, int len, const ubyte_t *str, bwtint_t *sa_begin, bwtint_t *sa_end);
	int bwt_match_exact_alt(const bwt_t *bwt, int len, const ubyte_t *str, bwtint_t *k0, bwtint_t *l0);

	/**
	 * Extend bi-SA-interval _ik_
	 */
	void bwt_extend(const bwt_t *bwt, const bwtintv_t *ik, bwtintv_t ok[4], int is_back);

	/**
	 * Given a query _q_, collect potential SMEMs covering position _x_ and store them in _mem_.
	 * Return the end of the longest exact match starting from _x_.
	 */
	int bwt_smem1(const bwt_t *bwt, int len, const uint8_t *q, int x, int min_intv, bwtintv_v *mem, bwtintv_v *tmpvec[2]);
	int bwt_smem1a(const bwt_t *bwt, int len, const uint8_t *q, int x, int min_intv, uint64_t max_intv, bwtintv_v *mem, bwtintv_v *tmpvec[2]);
	int bwt_smem2(const bwt_t *bwt, int len, const uint8_t *q, int min_intv, bwtintv_v *mem, bwtintv_v *tmpvec[2], bwtintv_t kmer_intv);

	int bwt_seed_strategy1(const bwt_t *bwt, int len, const uint8_t *q, int x, int min_len, int max_intv, bwtintv_t *mem);

#ifdef __cplusplus
}
#endif

#endif

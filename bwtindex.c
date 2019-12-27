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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#include <zlib.h>
#include "bntseq.h"
#include "bwa.h"
#include "bwt.h"
#include "utils.h"
#include "rle.h"
#include "rope.h"
#include "mmp2_minimap.h"


#ifdef _DIVBWT
#include "divsufsort.h"
#endif

#ifdef USE_MALLOC_WRAPPERS
#  include "malloc_wrap.h"
#endif


int is_bwt(ubyte_t *T, int n);

int64_t bwa_seq_len(const char *fn_pac)
{
	FILE *fp;
	int64_t pac_len;
	ubyte_t c;
	fp = xopen(fn_pac, "rb");
	err_fseek(fp, -1, SEEK_END);
	pac_len = err_ftell(fp);
	err_fread_noeof(&c, 1, 1, fp);
	err_fclose(fp);
	return (pac_len - 1) * 4 + (int)c;
}

/**
 * 计算pac序列的bwt变换。
 *
 * @param fn_pac        pac文件的路径
 * @param use_is        决定了计算bwt变换的算法
 * @return
 */
bwt_t *bwt_pac2bwt(const char *fn_pac, int use_is)
{
	bwt_t *bwt;
	ubyte_t *buf, *buf2; // buf用来存储原始序列，buf2用来存储pac序列。
	int64_t i, pac_size;
	FILE *fp;

	// initialization
	bwt = (bwt_t*)calloc(1, sizeof(bwt_t));
	bwt->seq_len = bwa_seq_len(fn_pac);
	bwt->bwt_size = (bwt->seq_len + 15) >> 4;
	fp = xopen(fn_pac, "rb");

	// prepare sequence
	pac_size = (bwt->seq_len>>2) + ((bwt->seq_len&3) == 0? 0 : 1); // 存放pac序列的buffer的大小
	buf2 = (ubyte_t*)calloc(pac_size, 1);
	err_fread_noeof(buf2, 1, pac_size, fp);
	err_fclose(fp);
	memset(bwt->L2, 0, 5 * 4);
	buf = (ubyte_t*)calloc(bwt->seq_len + 1, 1);
	for (i = 0; i < bwt->seq_len; ++i) {// 从pac序列中解析出原始序列存在buf中，求解C表（需要注意pac序列中不会保存N，N在pac中是用任意的某一碱基代替的）。
		buf[i] = buf2[i>>2] >> ((3 - (i&3)) << 1) & 3;
		++bwt->L2[1+buf[i]];
	}
	for (i = 2; i <= 4; ++i) bwt->L2[i] += bwt->L2[i-1];
	free(buf2);

	// Burrows-Wheeler Transform
	if (use_is) {
		bwt->primary = is_bwt(buf, bwt->seq_len);
	} else {
		rope_t *r;
		int64_t x;
		rpitr_t itr;
		const uint8_t *blk;

		r = rope_init(ROPE_DEF_MAX_NODES, ROPE_DEF_BLOCK_LEN);
		for (i = bwt->seq_len - 1, x = 0; i >= 0; --i) {
			int c = buf[i] + 1;
			x = rope_insert_run(r, x, c, 1, 0) + 1;
			while (--c >= 0) x += r->c[c];
		}
		bwt->primary = x;
		rope_itr_first(r, &itr);
		x = 0;
		while ((blk = rope_itr_next_block(&itr)) != 0) {
			const uint8_t *q = blk + 2, *end = blk + 2 + *rle_nptr(blk);
			while (q < end) {
				int c = 0;
				int64_t l;
				rle_dec1(q, c, l);
				for (i = 0; i < l; ++i)
					buf[x++] = c - 1;
			}
		}
		rope_destroy(r);
	}
	bwt->bwt = (uint32_t*)calloc(bwt->bwt_size, 4);
	for (i = 0; i < bwt->seq_len; ++i)
		bwt->bwt[i>>4] |= buf[i] << ((15 - (i&15)) << 1);
	free(buf);
	return bwt;
}

int bwa_pac2bwt(int argc, char *argv[]) // the "pac2bwt" command; IMPORTANT: bwt generated at this step CANNOT be used with BWA. bwtupdate is required!
{
	bwt_t *bwt;
	int c, use_is = 1;
	while ((c = getopt(argc, argv, "d")) >= 0) {
		switch (c) {
		case 'd': use_is = 0; break;
		default: return 1;
		}
	}
	if (optind + 2 > argc) {
		fprintf(stderr, "Usage: bwa pac2bwt [-d] <in.pac> <out.bwt>\n");
		return 1;
	}
	bwt = bwt_pac2bwt(argv[optind], use_is);
	bwt_dump_bwt(argv[optind+1], bwt);
	bwt_destroy(bwt);
	return 0;
}

#define bwt_B00(b, k) ((b)->bwt[(k)>>4]>>((~(k)&0xf)<<1)&3)

/**
 * 更新bwt，原来的bwt序列是单纯的序列的bwt变换，每2bit表示一个碱基。
 * 该函数把occurence信息添加到bwt序列中，128bp记录一次occurence信息。
 * 具体的格式如下：
 *      32 Byte: O表
 *      32 Byte: 128bp的bwt序列
 *      32 Byte: O表
 *      32 Byte: 128bp的bwt序列
 *      ...
 *      32 Byte: O表
 *      <32Byte: <128bp的bwt序列
 *      32 Byte: O表
 *
 * @param bwt       bwt变换
 */
void bwt_bwtupdate_core(bwt_t *bwt)
{
	bwtint_t i, k, c[4], n_occ;
	uint32_t *buf;

	n_occ = (bwt->seq_len + OCC_INTERVAL - 1) / OCC_INTERVAL + 1;
	bwt->bwt_size += n_occ * sizeof(bwtint_t); // the new size
	buf = (uint32_t*)calloc(bwt->bwt_size, 4); // will be the new bwt
	c[0] = c[1] = c[2] = c[3] = 0;
	for (i = k = 0; i < bwt->seq_len; ++i) {
		if (i % OCC_INTERVAL == 0) {
			memcpy(buf + k, c, sizeof(bwtint_t) * 4);
			k += sizeof(bwtint_t); // in fact: sizeof(bwtint_t)=4*(sizeof(bwtint_t)/4)
		}
		if (i % 16 == 0) buf[k++] = bwt->bwt[i/16]; // 16 == sizeof(uint32_t)/2
		++c[bwt_B00(bwt, i)];
	}
	// the last element, 每128bp记录一次c，记录的c是前128bp的c，因此最后<=128bp的c在这里记录。
	memcpy(buf + k, c, sizeof(bwtint_t) * 4);
	xassert(k + sizeof(bwtint_t) == bwt->bwt_size, "inconsistent bwt_size");
	// update bwt
	free(bwt->bwt); bwt->bwt = buf;
}

int bwa_bwtupdate(int argc, char *argv[]) // the "bwtupdate" command
{
	bwt_t *bwt;
	if (argc != 2) {
		fprintf(stderr, "Usage: bwa bwtupdate <the.bwt>\n");
		return 1;
	}
	bwt = bwt_restore_bwt(argv[1]);
	bwt_bwtupdate_core(bwt);
	bwt_dump_bwt(argv[1], bwt);
	bwt_destroy(bwt);
	return 0;
}

int bwa_bwt2sa(int argc, char *argv[]) // the "bwt2sa" command
{
	bwt_t *bwt;
	int c, sa_intv = 32;
	while ((c = getopt(argc, argv, "i:")) >= 0) {
		switch (c) {
		case 'i': sa_intv = atoi(optarg); break;
		default: return 1;
		}
	}
	if (optind + 2 > argc) {
		fprintf(stderr, "Usage: bwa bwt2sa [-i %d] <in.bwt> <out.sa>\n", sa_intv);
		return 1;
	}
	bwt = bwt_restore_bwt(argv[optind]);
	bwt_cal_sa(bwt, sa_intv);
	bwt_dump_sa(argv[optind+1], bwt);
	bwt_destroy(bwt);
	return 0;
}

int bwa_index(int argc, char *argv[]) // the "index" command
{
	int c, algo_type = BWTALGO_AUTO, is_64 = 0, block_size = 10000000;
	char *prefix = 0, *str;
	while ((c = getopt(argc, argv, "6a:p:b:")) >= 0) {
		switch (c) {
		case 'a': // if -a is not set, algo_type will be determined later
			if (strcmp(optarg, "rb2") == 0) algo_type = BWTALGO_RB2;
			else if (strcmp(optarg, "bwtsw") == 0) algo_type = BWTALGO_BWTSW;
			else if (strcmp(optarg, "is") == 0) algo_type = BWTALGO_IS;
			else err_fatal(__func__, "unknown algorithm: '%s'.", optarg);
			break;
		case 'p': prefix = strdup(optarg); break;
		case '6': is_64 = 1; break;
		case 'b':
			block_size = strtol(optarg, &str, 10);
			if (*str == 'G' || *str == 'g') block_size *= 1024 * 1024 * 1024;
			else if (*str == 'M' || *str == 'm') block_size *= 1024 * 1024;
			else if (*str == 'K' || *str == 'k') block_size *= 1024;
			break;
		default: return 1;
		}
	}

	if (optind + 1 > argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   bwa index [options] <in.fasta>\n\n");
		fprintf(stderr, "Options: -a STR    BWT construction algorithm: bwtsw, is or rb2 [auto]\n");
		fprintf(stderr, "         -p STR    prefix of the index [same as fasta name]\n");
		fprintf(stderr, "         -b INT    block size for the bwtsw algorithm (effective with -a bwtsw) [%d]\n", block_size);
		fprintf(stderr, "         -6        index files named as <in.fasta>.64.* instead of <in.fasta>.* \n");
		fprintf(stderr, "\n");
		fprintf(stderr,	"Warning: `-a bwtsw' does not work for short genomes, while `-a is' and\n");
		fprintf(stderr, "         `-a div' do not work not for long genomes.\n\n");
		return 1;
	}
	if (prefix == 0) { // prefix 存储的是fasta文件的路径
		prefix = malloc(strlen(argv[optind]) + 4);
		strcpy(prefix, argv[optind]);
		if (is_64) strcat(prefix, ".64");
	}
	bwa_idx_build(argv[optind], prefix, algo_type, block_size); // 生成索引
	free(prefix);
	return 0;
}

/**
 * 生成索引，过程主要分为5步。
 * step1: 读取fasta文件，生成.ann .amb 和 带反向互补序列的.pac文件
 * step2: 根据step1生成的.pac文件计算序列的bwt变换，生成.bwt文件
 * step3: 向step2中得到的bwt变换中添加occ信息，更新.bwt文件
 * step4: 读取fasta文件，重新生成不带反向互补序列的.pac .ann .amb文件，覆盖掉step1中生成的文件
 * step5: 从.bwt文件计算后缀数组sa，生成.sa文件
 *
 * @param fa			fasta文件的文件路径
 * @param prefix		待生成的索引文件的前缀
 * @param algo_type		algorithm type，指定生成BWT的算法；默认为0，表示不指定算法，而是根据fasta文件中序列的长度决定采用哪种算法，
 * @param block_size	默认为10,000,000
 * @return 0
 */
int bwa_idx_build(const char *fa, const char *prefix, int algo_type, int block_size)
{
	extern void bwa_pac_rev_core(const char *fn, const char *fn_rev);

	char *str, *str2, *str3;
	clock_t t;
	int64_t l_pac;

	str  = (char*)calloc(strlen(prefix) + 10, 1);
	str2 = (char*)calloc(strlen(prefix) + 10, 1);
	str3 = (char*)calloc(strlen(prefix) + 10, 1);

	// step 1
	{ // nucleotide indexing
		gzFile fp = xzopen(fa, "r");// 尝试调用zlib库打开FASTA的压缩文件，同时也可以打开非压缩文件，fp是gzFile类型的
		t = clock();
		if (bwa_verbose >= 3) fprintf(stderr, "[bwa_index] Pack FASTA... ");
		l_pac = bns_fasta2bntseq(fp, prefix, 0);// 生成.pac, .ann和.amb文件，返回pac文件中bp的长度。
		if (bwa_verbose >= 3) fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
		err_gzclose(fp);
	}
	if (algo_type == 0) algo_type = l_pac > 50000000? 2 : 3; // set the algorithm for generating BWT，algo_type值为0表示没有指定算法，那么就根据bp的长度来确定算法

	// step 2
	{
		strcpy(str, prefix); strcat(str, ".pac");
		strcpy(str2, prefix); strcat(str2, ".bwt");
		t = clock();
		if (bwa_verbose >= 3) fprintf(stderr, "[bwa_index] Construct BWT for the packed sequence...\n");
		if (algo_type == 2)
		    bwt_bwtgen2(str, str2, block_size);
		else if (algo_type == 1 || algo_type == 3) {
			bwt_t *bwt;
			bwt = bwt_pac2bwt(str, algo_type == 3);
			bwt_dump_bwt(str2, bwt);
			bwt_destroy(bwt);
		}
		if (bwa_verbose >= 3) fprintf(stderr, "[bwa_index] %.2f seconds elapse.\n", (float)(clock() - t) / CLOCKS_PER_SEC);
	}

	// step 3
	{
		bwt_t *bwt;
		strcpy(str, prefix); strcat(str, ".bwt");
		t = clock();
		if (bwa_verbose >= 3) fprintf(stderr, "[bwa_index] Update BWT... ");
		bwt = bwt_restore_bwt(str);
		bwt_bwtupdate_core(bwt);
		bwt_dump_bwt(str, bwt);
		bwt_destroy(bwt);
		if (bwa_verbose >= 3) fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
	}

	// step 4
	{
		gzFile fp = xzopen(fa, "r");
		t = clock();
		if (bwa_verbose >= 3) fprintf(stderr, "[bwa_index] Pack forward-only FASTA... ");
		l_pac = bns_fasta2bntseq(fp, prefix, 1);
		if (bwa_verbose >= 3) fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
		err_gzclose(fp);
	}

	// step 5
    bwt_t *bwt;
	{

		strcpy(str, prefix); strcat(str, ".bwt");
		strcpy(str3, prefix); strcat(str3, ".sa");
		t = clock();
		if (bwa_verbose >= 3) fprintf(stderr, "[bwa_index] Construct SA from BWT and Occ... ");
		bwt = bwt_restore_bwt(str);
		bwt_cal_sa(bwt, 32);
		bwt_dump_sa(str3, bwt);
		if (bwa_verbose >= 3) fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
	}

	// step 6
    {
        /************新建哈希索引*************/
        char *fnw;
        mm_idx_reader_t *idx_rdr;
        mm_idxopt_t ipt;
        mm_idx_t *mi;
        int n_threads = 1;

        ipt.batch_size = 4000000000;
        ipt.bucket_bits = 14;
        ipt.flag = 0;
        ipt.k = 21;
        ipt.mini_batch_size = 50000000;
        ipt.w = 11;

        fnw = (char*)malloc(strlen(fa) + 10);
        strcpy(fnw, fa);
        strcat(fnw, ".mmi2");

        idx_rdr = mm_idx_reader_open(fa, &ipt, fnw);
        if (idx_rdr == 0) {
            fprintf(stderr, "[ERROR] failed to open file '%s'\n", fa);
            return 1;
        }
        mi = mm_idx_reader_read(idx_rdr, n_threads, bwt);
        mm_idx_reader_close(idx_rdr);
        /***********************************/
    }
    bwt_destroy(bwt);
	free(str3); free(str2); free(str);
	return 0;
}

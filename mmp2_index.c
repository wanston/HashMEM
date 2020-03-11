#include <stdlib.h>
#include <assert.h>
#if defined(WIN32) || defined(_WIN32)
#include <io.h> // for open(2)
#else
#include <unistd.h>
#endif
#include <fcntl.h>
#include <stdio.h>
#define __STDC_LIMIT_MACROS
#include "mmp2_kthread.h"
#include "mmp2_bseq.h"
#include "mmp2_minimap.h"
#include "mmp2_mmpriv.h"
#include "mmp2_kvec.h"
#include "mmp2_khash.h"
#include "bwt.h"

#define idx_hash(a) ((a))
#define idx_eq(a, b) ((a) == (b))
KHASH_INIT(idx, uint64_t, bwtintv_x_t, 1, idx_hash, idx_eq)
typedef khash_t(idx) idxhash_t;

KHASH_MAP_INIT_STR(str, uint32_t)

#define kroundup64(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, (x)|=(x)>>32, ++(x))


/**
 * 关键的索引结构
 * a    (minimizer, position) array，fasta序列生成的mm128_t结构体存放于此，存放于此的只是所有(minmizer, position)的一部分，
 * 		因为有多个bucket，存放与此的minimizer只是所有minimizer的一部分。
 *      当p、h建立好之后，(minimizer，position)就会存在p、h中，此后a的内容就会变空。
 * n    position array的大小
 * p    position array，记录的是出现次数大于一次的minimizer的position；同一minimizer的positions被升序连续存放，并且minimizer
 *      之间也是升序的。比如，minimizer 1 的位置有77 89， minimizer 2的位置有35 70，那么position数组的内容就是77、89、35、70。
 * h    哈希表，int64_t到int64_t的映射
 *      key:  高63bit：mm128_t::x>>8>>mi->b<<1，低1bit：表示minimizer是否仅出现一次。
 *      val:  如果minimizer仅出现一次，val记录的是minimizer的position，即mm128_t::y；
 *            如果minimizer出现多次，position则有多个值，存储在p中，此时val的高32bit是position在p中的起始索引，低32bit是position的数目。
 *      由于该哈希表的实现，key的最低1bit不影响key的hash的结果。
 * */
typedef struct mm_idx_bucket_s {
	mm256_v a;
//	int32_t n;   // size of the _p_ array
//	uint64_t *p; // position array for minimizers appearing >1 times
	void *h;     // hash table indexing _p_ and minimizers appearing once
} mm_idx_bucket_t;

mm_idx_t *mm_idx_init(int w, int k, int b, int flag)
{
	mm_idx_t *mi;
	if (k*2 < b) b = k * 2;
	if (w < 1) w = 1;
	mi = (mm_idx_t*)calloc(1, sizeof(mm_idx_t));
	mi->w = w, mi->k = k, mi->b = b, mi->flag = flag;
	mi->B = (mm_idx_bucket_t*)calloc(1<<b, sizeof(mm_idx_bucket_t));
	if (!(mm_dbg_flag & 1)) mi->km = km_init();
	return mi;
}

void mm_idx_destroy(mm_idx_t *mi)
{
	uint32_t i;
	if (mi == 0) return;
	if (mi->h) kh_destroy(str, (khash_t(str)*)mi->h);
	if (mi->B) {
		for (i = 0; i < 1U<<mi->b; ++i) {
			free(mi->B[i].a.a);
			kh_destroy(idx, (idxhash_t*)mi->B[i].h);
		}
	}
	if (!mi->km) {
		for (i = 0; i < mi->n_seq; ++i)
			free(mi->seq[i].name);
		free(mi->seq);
	} else km_destroy(mi->km);
	free(mi->B); free(mi->S); free(mi);
}

/**
 *  根据minimizer查索引，返回索引中记录的minimizer的位置信息
 *
 *  @param mi           索引
 *  @param minier       minimizer，其有效位数为2k，k是kmer长度，因为minmizer是窗口中最小的kmer经过hash得到的，hash前后kmer的位数不变。
 *  @param n            记录返回的uint64_t动态数组的大小，每个uint64_t值记录minimizer的位置信息，应该是对应mm128_t中的y。
 *  @retrun             返回的指针指向动态数组
 * */
bwtintv_x_t mm_idx_get(const mm_idx_t *mi, uint64_t minier)
{
	int mask = (1<<mi->b) - 1; // minimizer的哈希值的mask
	khint_t k;
	bwtintv_x_t ret = {0,0,0};
	mm_idx_bucket_t *b = &mi->B[minier&mask]; //
	idxhash_t *h = (idxhash_t*)b->h;
	if (h == 0) return ret;
	k = kh_get(idx, h, minier>>mi->b);
	if (k == kh_end(h)) return ret;
	ret = kh_val(h, k);
	return ret;
}

/**
 * 统计索引的状态信息
 *
 * @param mi
 */
void mm_idx_stat(const mm_idx_t *mi)
{
	int n = 0, n1 = 0;
	uint32_t i;
	uint64_t sum = 0, len = 0;
	fprintf(stderr, "[M::%s] kmer size: %d; skip: %d; is_hpc: %d; #seq: %d\n", __func__, mi->k, mi->w, mi->flag&MM_I_HPC, mi->n_seq);
	for (i = 0; i < mi->n_seq; ++i)
		len += mi->seq[i].len;
	for (i = 0; i < 1U<<mi->b; ++i)
		if (mi->B[i].h) n += kh_size((idxhash_t*)mi->B[i].h);
	for (i = 0; i < 1U<<mi->b; ++i) {
		idxhash_t *h = (idxhash_t*)mi->B[i].h;
		khint_t k;
		if (h == 0) continue;
		for (k = 0; k < kh_end(h); ++k)
			if (kh_exist(h, k)) {
				++sum;
				++n1;
			}
	}
	fprintf(stderr, "[M::%s::%.3f*%.2f] distinct minimizers: %d (%.2f%% are singletons); average occurrences: %.3lf; average spacing: %.3lf\n",
			__func__, realtime() - mm_realtime0, cputime() / (realtime() - mm_realtime0), n, 100.0*n1/n, (double)sum / n, (double)len / sum);
}

int mm_idx_index_name(mm_idx_t *mi)
{
	khash_t(str) *h;
	uint32_t i;
	int has_dup = 0, absent;
	if (mi->h) return 0;
	h = kh_init(str);
	for (i = 0; i < mi->n_seq; ++i) {
		khint_t k;
		k = kh_put(str, h, mi->seq[i].name, &absent);
		if (absent) kh_val(h, k) = i;
		else has_dup = 1;
	}
	mi->h = h;
	if (has_dup && mm_verbose >= 2)
		fprintf(stderr, "[WARNING] some database sequences have identical sequence names\n");
	return has_dup;
}

int mm_idx_name2id(const mm_idx_t *mi, const char *name)
{
	khash_t(str) *h = (khash_t(str)*)mi->h;
	khint_t k;
	if (h == 0) return -2;
	k = kh_get(str, h, name);
	return k == kh_end(h)? -1 : kh_val(h, k);
}

int mm_idx_getseq(const mm_idx_t *mi, uint32_t rid, uint32_t st, uint32_t en, uint8_t *seq)
{
	uint64_t i, st1, en1;
	if (rid >= mi->n_seq || st >= mi->seq[rid].len) return -1;
	if (en > mi->seq[rid].len) en = mi->seq[rid].len;
	st1 = mi->seq[rid].offset + st;
	en1 = mi->seq[rid].offset + en;
	for (i = st1; i < en1; ++i)
		seq[i - st1] = mm_seq4_get(mi->S, i);
	return en - st;
}

int32_t mm_idx_cal_max_occ(const mm_idx_t *mi, float f)
{
	int i;
	size_t n = 0;
	uint32_t thres;
	khint_t *a, k;
	if (f <= 0.) return INT32_MAX;
	for (i = 0; i < 1<<mi->b; ++i)
		if (mi->B[i].h) n += kh_size((idxhash_t*)mi->B[i].h);
	a = (uint32_t*)malloc(n * 4);
	for (i = n = 0; i < 1<<mi->b; ++i) {
		idxhash_t *h = (idxhash_t*)mi->B[i].h;
		if (h == 0) continue;
		for (k = 0; k < kh_end(h); ++k) {
			if (!kh_exist(h, k)) continue;
			a[n++] = 1;
		}
	}
	thres = ks_ksmall_uint32_t(n, a, (uint32_t)((1. - f) * n)) + 1;
	free(a);
	return thres;
}

bwtintv_x_t get_kmer_interval(uint64_t seq, uint8_t kmer_span, const bwt_t *bwt){
    int i, j, bp, bp_c;
    bwtintv_t ik, ok[4]; // 我要存储的是kv_push的结果
    bwtintv_x_t ret = {0, 0, 0};
    bp = seq >> (2*kmer_span-2) & 3;
    if (bp > 3){
        fprintf(stderr, "WT: In getInterval: bp > 3\n");
        return ret;
    }
    bwt_set_intv(bwt, bp, ik);

    for (i = 1; i < kmer_span; ++i) { // forward search
        bp = (seq >> ((kmer_span-1-i)*2)) & 3;
        bp_c = 3 - bp;// complement of q[i], 因为是forward extention，所以用互补碱基
        if (bp < 4) { // an A/C/G/T base
            bwt_extend(bwt, &ik, ok, 0);
            ik = ok[bp_c];
            if(ik.x[2] <= 0){
                break;
            }
        }else{
            fprintf(stderr, "WT: In getInterval: bp > 3\n");
            return ret;
        }
    }
    ret.x[0] = ik.x[0]; ret.x[1] = ik.x[1]; ret.x[2] = ik.x[2];
	return ret;
}

/*********************************
 * Sort and generate hash tables *
 * 此函数仅仅处理mi中的一个bucket，bucket的所在数组的索引是i。该函数负责
 * 把bucket中的a中的（minimizer, position）信息存入到bucket的h和p中，
 * 也就是把a中的元素存入一个独立的哈希表中。
 *
 * @param g 		实际为mm_idx_t* mi,指向要处理的表示索引的数据结构，此时fasta序列已经生成为mm128_t存入mi的各个bucket中。
 * @param i 		用作mi->B[i]，本函数只处理该bucket
 * @param tid		线程id
 *********************************/

static void worker_post(void *g, long i, int tid)
{
	int n, n_keys;
	size_t j, start_a, start_p;
	idxhash_t *h;
	mm_idx_t *mi = (mm_idx_t*)((void**)g)[0];
    bwt_t* bwt = (bwt_t*)((void**)g)[1];
	mm_idx_bucket_t *b = &mi->B[i];
	if (b->a.n == 0) return;

	// sort by minimizer，升序排序，根据mm128_t的x比较大小。此时，相同的minimzer被连续放置，
	// 相同的minimizer的来源：可能是同一窗口内的最小kmer有多个，也可能是不同窗口的最小kmer相同。
	radix_sort_256x(b->a.a, b->a.a + b->a.n);

	// count and preallocate，该段代码执行后，b->n记录存在重复的minimizer的总数，n_keys则统计bucket中的minimizer有几种。
	// 比如说，bucket中的(minimizer,pos)，是(1,x) (1,y) (2,z) (3,a) (3,b)，那么，b->n的值就是2+2=4，n_keys的值是3。
	for (j = 1, n = 1, n_keys = 0; j <= b->a.n; ++j) {
		if (j == b->a.n || b->a.a[j].x>>8 != b->a.a[j-1].x>>8) {
			++n_keys;
		}
	}
	h = kh_init(idx);
	kh_resize(idx, h, n_keys);

	// create the hash table
	for (j = 1; j <= b->a.n; ++j) {
		if (j == b->a.n || b->a.a[j].x>>8 != b->a.a[j-1].x>>8) {
			khint_t itr;
			int absent;
			mm256_t *p = &b->a.a[j-1];
			itr = kh_put(idx, h, p->x>>8>>mi->b, &absent);
			assert(absent);

			uint8_t kmer_span = (uint8_t)p->x;
            bwtintv_x_t intv = get_kmer_interval(p->y[0], kmer_span, bwt);// 低kmer_span*2位存储的是kmer。其中高位表示左边的碱基，低位表示右边的碱基。
            if(intv.x[2] > 0) {
                kh_val(h, itr) = intv; // TODO：建表和查询的数据结构需要再精简一下。
            }
		}
	}
	b->h = h;

	// deallocate and clear b->a
	kfree(0, b->a.a);
	b->a.n = b->a.m = 0, b->a.a = 0;
}
 
static void mm_idx_post(mm_idx_t *mi, int n_threads, const bwt_t *bwt)
{
	const void *data[2] = {mi, bwt};
    kt_for(n_threads, worker_post, data, 1<<mi->b);
}

/******************
 * Generate index *
 ******************/

#include <string.h>
#include <zlib.h>
#include "mmp2_bseq.h"
#include "utils.h"

typedef struct {
	int mini_batch_size; 			// 一批处理的大小
	uint64_t batch_size, sum_len; 	// 已经处理过的序列的总长度（bp）
	mm_bseq_file_t *fp; 			// 序列文件
	mm_idx_t *mi; 					// 索引
} pipeline_t;

typedef struct {
    int n_seq; // 数组的大小
	mm_bseq1_t *seq; // 动态数组，每个元素表示一条序列（目前已知的是fasta的序列）
	mm256_v a; // 保存seq计算出来的kmer和interval信息
} step_t; // 用于存储每一step中处理的序列，也就是一个batch的序列


/**
 * 把数组a中的mm128_t元素添加进mi中对应的bucket。
 * bucket的确定根据minimizer的哈希值的低b位。
 * */
static void mm_idx_add(mm_idx_t *mi, int n, const mm256_t *a)
{
	int i, mask = (1<<mi->b) - 1;
	for (i = 0; i < n; ++i) {
		mm256_v *p = &mi->B[a[i].x>>8&mask].a;
		kv_push(mm256_t, 0, *p, a[i]);
	}
}


/**
 * 构建索引的函数，分为三个step，step0 读取一批序列；step1 计算sketch；step2 保存sketch。
 *
 * @param shared: 指向pipeline_t结构体
 * @param step: step编号，取值0 1 2
 * @param in: step0的时候用不到；step1和step2的时候，in指向step_t结构体，该结构体保存的是step0中读取的一批序列。
 *
 * @return step0和step1的时候，返回的是step_t*，记录了这两步中用到的序列；step2的时候，返回0。
 * */
static void *worker_pipeline(void *shared, int step, void *in)
{
	int i;
    pipeline_t *p = (pipeline_t*)shared;
    if (step == 0) { // step 0: read sequences
        step_t *s;
		if (p->sum_len > p->batch_size) return 0;
        s = (step_t*)calloc(1, sizeof(step_t));
		s->seq = mm_bseq_read(p->fp, p->mini_batch_size, 0, &s->n_seq); // read a mini-batch // 一个batch多少read会存在s->n_se	q中。
		if (s->seq) {
			uint32_t old_m, m;
			assert((uint64_t)p->mi->n_seq + s->n_seq <= UINT32_MAX); // to prevent integer overflow
			// make room for p->mi->seq
			old_m = p->mi->n_seq, m = p->mi->n_seq + s->n_seq;
			kroundup32(m); kroundup32(old_m); // kroundup32宏的作用是，把输入的m，向上二进制取整，例如，0101的数取整后就是1000；但是如果是2^31+1，那么取整后的结果是2^32，就overflow了。
			if (old_m != m)
				p->mi->seq = (mm_idx_seq_t*)krealloc(p->mi->km, p->mi->seq, m * sizeof(mm_idx_seq_t));
			// make room for p->mi->S
			if (!(p->mi->flag & MM_I_NO_SEQ)) {
				uint64_t sum_len, old_max_len, max_len;
				for (i = 0, sum_len = 0; i < s->n_seq; ++i) sum_len += s->seq[i].l_seq;
				old_max_len = (p->sum_len + 7) / 8;
				max_len = (p->sum_len + sum_len + 7) / 8;
				kroundup64(old_max_len); kroundup64(max_len);
				if (old_max_len != max_len) {
					p->mi->S = (uint32_t*)realloc(p->mi->S, max_len * 4);
					memset(&p->mi->S[old_max_len], 0, 4 * (max_len - old_max_len));
				}
			}
			// populate p->mi->seq
			for (i = 0; i < s->n_seq; ++i) {
				mm_idx_seq_t *seq = &p->mi->seq[p->mi->n_seq];
				uint32_t j;
				if (!(p->mi->flag & MM_I_NO_NAME)) {
					seq->name = (char*)kmalloc(p->mi->km, strlen(s->seq[i].name) + 1);
					strcpy(seq->name, s->seq[i].name);
				} else seq->name = 0;
				seq->len = s->seq[i].l_seq;
				seq->offset = p->sum_len;
				// copy the sequence
				if (!(p->mi->flag & MM_I_NO_SEQ)) {
					for (j = 0; j < seq->len; ++j) { // TODO: this is not the fastest way, but let's first see if speed matters here
						uint64_t o = p->sum_len + j;
						int c = seq_nt4_table[(uint8_t)s->seq[i].seq[j]]; // 把char转为int，值只有0到4
						mm_seq4_set(p->mi->S, o, c);
					}
				}
				// update p->sum_len and p->mi->n_seq
				p->sum_len += seq->len;
				s->seq[i].rid = p->mi->n_seq++;
			}
			return s;
		} else free(s);
    } else if (step == 1) { // step 1: compute sketch
        step_t *s = (step_t*)in;
		for (i = 0; i < s->n_seq; ++i) {
			mm_bseq1_t *t = &s->seq[i];
			if (t->l_seq > 0){
				mm_sketch_intv_char(0, t->seq, t->l_seq, p->mi->w, p->mi->k, p->mi->flag&MM_I_HPC, &s->a);
			}
			else if (mm_verbose >= 2)
				fprintf(stderr, "[WARNING] the length database sequence '%s' is 0\n", t->name);
			free(t->seq); free(t->name);
		}
		free(s->seq); s->seq = 0;
		return s;
    } else if (step == 2) { // step 2: dispatch sketch to buckets
        step_t *s = (step_t*)in;
		mm_idx_add(p->mi, s->a.n, s->a.a);
		kfree(0, s->a.a); free(s);
	}
    return 0;
}

/**
 * Generate Index (mmi file)
 *
 * @param fp	fasta handle
 * @param io         pointer to indexing parameters
 * @param mo         pointer to mapping parameters
 *
 * @return 0 if success; -1 if _present_ unknown
 */
mm_idx_t *mm_idx_gen(mm_bseq_file_t *fp, int w, int k, int b, int flag, int mini_batch_size, int n_threads, uint64_t batch_size, const bwt_t *bwt)
{
	pipeline_t pl;
	if (fp == 0 || mm_bseq_eof(fp)) return 0;
	memset(&pl, 0, sizeof(pipeline_t));
	pl.mini_batch_size = (uint64_t)mini_batch_size < batch_size? mini_batch_size : batch_size;
	pl.batch_size = batch_size;
	pl.fp = fp;
	pl.mi = mm_idx_init(w, k, b, flag);

	kt_pipeline(n_threads < 3? n_threads : 3, worker_pipeline, &pl, 3);
	if (mm_verbose >= 3)
		fprintf(stderr, "[M::%s::%.3f*%.2f] collected minimizers\n", __func__, realtime() - mm_realtime0, cputime() / (realtime() - mm_realtime0));

	mm_idx_post(pl.mi, n_threads, bwt);
	if (mm_verbose >= 3)
		fprintf(stderr, "[M::%s::%.3f*%.2f] sorted minimizers\n", __func__, realtime() - mm_realtime0, cputime() / (realtime() - mm_realtime0));

	return pl.mi;
}

mm_idx_t *mm_idx_build(const char *fn, int w, int k, int flag, int n_threads, const bwt_t *bwt) // a simpler interface; deprecated
{
	mm_bseq_file_t *fp;
	mm_idx_t *mi;
	fp = mm_bseq_open(fn);
	if (fp == 0) return 0;
	mi = mm_idx_gen(fp, w, k, 14, flag, 1<<18, n_threads, UINT64_MAX, bwt);
	mm_bseq_close(fp);
	return mi;
}

/**
 * index I/O，把索引保存成文件
 *
 * @param fp
 * @param mi
 */
void mm_idx_dump(FILE *fp, const mm_idx_t *mi)
{
	uint64_t sum_len = 0;
	uint32_t x[4], i;

	x[0] = mi->w, x[1] = mi->k, x[2] = mi->b,  x[3] = mi->flag;
	fwrite(MM_IDX_MAGIC, 1, 4, fp);
	fwrite(x, 4, 4, fp);

	for (i = 0; i < 1<<mi->b; ++i) {
		mm_idx_bucket_t *b = &mi->B[i];
		khint_t k;
		idxhash_t *h = (idxhash_t*)b->h;
		uint32_t size = h? h->size : 0;
		fwrite(&size, 4, 1, fp);
		if (size == 0) continue;
		for (k = 0; k < kh_end(h); ++k) {
			uint64_t x;
			bwtintv_x_t y;
			if (!kh_exist(h, k)) continue;
			x = kh_key(h, k);
			y = kh_val(h, k);
			fwrite(&x, 8, 1,  fp);
			fwrite(&y, sizeof(bwtintv_x_t), 1, fp);
			LOG(stderr, "mm_idx_dump: %lx\n", (x<<mi->b) + i);
        }
	}

	fflush(fp);
}

mm_idx_t *mm_idx_load(FILE *fp)
{
	char magic[4];
	uint32_t x[4], i;
	mm_idx_t *mi;

	if (fread(magic, 1, 4, fp) != 4) return 0;
	if (strncmp(magic, MM_IDX_MAGIC, 4) != 0) return 0;
	if (fread(x, 4, 4, fp) != 4) return 0;
	mi = mm_idx_init(x[0], x[1], x[2], x[3]);

	for (i = 0; i < 1<<mi->b; ++i) {
		mm_idx_bucket_t *b = &mi->B[i];
		uint32_t j, size;
		khint_t k;
		idxhash_t *h;
		fread(&size, 4, 1, fp);
		if (size == 0) continue;
		b->h = h = kh_init(idx);
		kh_resize(idx, h, size);
		for (j = 0; j < size; ++j) {
			uint64_t x;
			bwtintv_x_t y;
			int absent;
			fread(&x, 8, 1, fp);
			fread(&y, sizeof(bwtintv_x_t), 1, fp);
			k = kh_put(idx, h, x, &absent);
			assert(absent);
			kh_val(h, k) = y;
			LOG(stderr, "mm_idx_load: %lx\n", (x<<mi->b) + i);
		}
	}
	return mi;
}

int64_t mm_idx_is_idx(const char *fn)
{
	int fd, is_idx = 0;
	int64_t ret, off_end;
	char magic[4];

	if (strcmp(fn, "-") == 0) return 0; // read from pipe; not an index
	fd = open(fn, O_RDONLY);
	if (fd < 0) return -1; // error
#ifdef WIN32
	if ((off_end = _lseeki64(fd, 0, SEEK_END)) >= 4) {
		_lseeki64(fd, 0, SEEK_SET);
#else
	if ((off_end = lseek(fd, 0, SEEK_END)) >= 4) {
		lseek(fd, 0, SEEK_SET);
#endif // WIN32
		ret = read(fd, magic, 4);
		if (ret == 4 && strncmp(magic, MM_IDX_MAGIC, 4) == 0)
			is_idx = 1;
	}
	close(fd);
	return is_idx? off_end : 0;
}

mm_idx_reader_t *mm_idx_reader_open(const char *fn, const mm_idxopt_t *opt, const char *fn_out)
{
	int64_t is_idx;
	mm_idx_reader_t *r;
	is_idx = mm_idx_is_idx(fn);
	if (is_idx < 0) return 0; // failed to open the index
	r = (mm_idx_reader_t*)calloc(1, sizeof(mm_idx_reader_t));
	r->is_idx = is_idx;
	if (opt) r->opt = *opt;
	else mm_idxopt_init(&r->opt);
	if (r->is_idx) {
		r->fp.idx = fopen(fn, "rb");
		r->idx_size = is_idx;
	} else r->fp.seq = mm_bseq_open(fn);
	if (fn_out) r->fp_out = fopen(fn_out, "wb");
	return r;
}

void mm_idx_reader_close(mm_idx_reader_t *r)
{
	if (r->is_idx) fclose(r->fp.idx);
	else mm_bseq_close(r->fp.seq);
	if (r->fp_out) fclose(r->fp_out);
	free(r);
}

/**
 * 读取索引，返回mm_idx_t指针，该结构体表示reference的索引。
 * 如果命令行中指定的是索引文件（.mmi），那么就会直接加载索引文件；
 * 如果命令行中的是reference文件（.fasta），那么会根据fasta文件计算得到索引。
 *
 * @return 如果参数r已经执行过一次本函数，那么再次执行本函数会返回0。
 * */
mm_idx_t *mm_idx_reader_read(mm_idx_reader_t *r, int n_threads, const bwt_t *bwt)
{
	mm_idx_t *mi;
	if (r->is_idx) {
		mi = mm_idx_load(r->fp.idx);
		if (mi && mm_verbose >= 2 && (mi->k != r->opt.k || mi->w != r->opt.w || (mi->flag&MM_I_HPC) != (r->opt.flag&MM_I_HPC)))
			fprintf(stderr, "[WARNING]\033[1;31m Indexing parameters (-k, -w or -H) overridden by parameters used in the prebuilt index.\033[0m\n");
	} else
		mi = mm_idx_gen(r->fp.seq, r->opt.w, r->opt.k, r->opt.bucket_bits, r->opt.flag, r->opt.mini_batch_size, n_threads, r->opt.batch_size, bwt);
	if (mi) {
		if (r->fp_out) mm_idx_dump(r->fp_out, mi); // 把index保存成文件
		mi->index = r->n_parts++;
	}
	return mi;
}

int mm_idx_reader_eof(const mm_idx_reader_t *r) // TODO: in extremely rare cases, mm_bseq_eof() might not work
{
	return r->is_idx? (feof(r->fp.idx) || ftell(r->fp.idx) == r->idx_size) : mm_bseq_eof(r->fp.seq);
}

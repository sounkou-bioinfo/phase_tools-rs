// phase_mnv: minimal phased-variant to MNV/complex haplotype builder.
//
// This intentionally keeps only the useful construction core: read a phased
// VCF/BCF with htslib, collect alternate alleles carried on each phased
// haplotype of one sample, join adjacent (or near-adjacent with --max-gap)
// alleles on the same phase set, and emit merged records using bases fetched
// from a FASTA reference.  Pure SNV blocks are MNVs; blocks containing indels
// are emitted as TYPE=COMPLEX.

#define _GNU_SOURCE

#include <ctype.h>
#include <errno.h>
#include <getopt.h>
#include <inttypes.h>
#include <limits.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>

#include <htslib/faidx.h>
#include <htslib/hts.h>
#include <htslib/vcf.h>

#define PS_MISSING LLONG_MIN

#ifndef bcf_int32_missing
#define bcf_int32_missing INT32_MIN
#endif
#ifndef bcf_int32_vector_end
#define bcf_int32_vector_end (INT32_MIN + 1)
#endif

typedef struct {
    const char *input_path;
    const char *fasta_path;
    const char *output_path;
    const char *sample_name;
    int64_t max_gap;
    int min_variants;
    bool no_ref_check;
    bool no_header;
    bool quiet;
} Config;

typedef struct {
    int32_t rid;
    int hap;          // 0 or 1 for diploid selected sample
    int64_t ps;       // FORMAT/PS or PS_MISSING
    int64_t pos;      // 1-based VCF POS
    int64_t end;      // 1-based inclusive REF end
    char *ref;        // upper-case REF allele
    char *alt;        // upper-case selected ALT allele
    bool is_snv;      // REF and ALT are both length 1
} Obs;

typedef struct {
    Obs *data;
    size_t n;
    size_t m;
} ObsVec;

typedef struct {
    int32_t rid;
    int64_t start;    // 1-based inclusive
    int64_t end;      // 1-based inclusive
    char *ref_seq;
    char *alt_seq;
    char *positions;  // comma-separated source variant positions
    int nvars;
    int nsnps;
    const char *type; // MNV or COMPLEX
    int hap_mask;     // bit 0 = first phased haplotype, bit 1 = second
    int64_t ps;
} MnvCall;

typedef struct {
    MnvCall *data;
    size_t n;
    size_t m;
} CallVec;

typedef struct {
    uint64_t records;
    uint64_t phased_records;
    uint64_t observations;
    uint64_t skipped_no_gt;
    uint64_t skipped_not_diploid;
    uint64_t skipped_missing_gt;
    uint64_t skipped_unphased;
    uint64_t skipped_ref;
    uint64_t skipped_unsupported_alt;
    uint64_t skipped_ref_allele;
    uint64_t emitted;
} Stats;

static void die(const char *msg) {
    fprintf(stderr, "error: %s\n", msg);
    exit(EXIT_FAILURE);
}

static void *xmalloc(size_t n) {
    void *p = malloc(n ? n : 1);
    if (!p) die("out of memory");
    return p;
}

static void *xrealloc(void *p, size_t n) {
    void *q = realloc(p, n ? n : 1);
    if (!q) die("out of memory");
    return q;
}

static char *xstrdup(const char *s) {
    if (!s) return NULL;
    size_t n = strlen(s) + 1;
    char *p = (char *)xmalloc(n);
    memcpy(p, s, n);
    return p;
}

static char upbase(char c) {
    return (char)toupper((unsigned char)c);
}

static char *xstrdup_upper(const char *s) {
    char *p = xstrdup(s);
    for (char *q = p; *q; ++q) *q = upbase(*q);
    return p;
}

static bool is_symbolic_or_breakend(const char *allele) {
    if (!allele || !allele[0]) return true;
    return allele[0] == '<' || strchr(allele, '[') != NULL || strchr(allele, ']') != NULL;
}

static bool is_dna_base(char c) {
    switch (upbase(c)) {
        case 'A': case 'C': case 'G': case 'T': case 'N':
            return true;
        default:
            return false;
    }
}

static bool is_plain_dna_allele(const char *allele) {
    if (!allele || !allele[0] || is_symbolic_or_breakend(allele)) return false;
    for (const char *p = allele; *p; ++p) {
        if (!is_dna_base(*p)) return false;
    }
    return true;
}

static bool is_snv_allele_pair(const char *ref, const char *alt) {
    return ref && alt && ref[0] && alt[0] && ref[1] == '\0' && alt[1] == '\0';
}

static void obs_push(ObsVec *v, Obs x) {
    if (v->n == v->m) {
        v->m = v->m ? v->m * 2 : 1024;
        v->data = (Obs *)xrealloc(v->data, v->m * sizeof(v->data[0]));
    }
    v->data[v->n++] = x;
}

static void call_push(CallVec *v, MnvCall x) {
    if (v->n == v->m) {
        v->m = v->m ? v->m * 2 : 256;
        v->data = (MnvCall *)xrealloc(v->data, v->m * sizeof(v->data[0]));
    }
    v->data[v->n++] = x;
}

static void free_obs(ObsVec *v) {
    for (size_t i = 0; i < v->n; ++i) {
        free(v->data[i].ref);
        free(v->data[i].alt);
    }
    free(v->data);
    v->data = NULL;
    v->n = v->m = 0;
}

static void free_call(MnvCall *c) {
    free(c->ref_seq);
    free(c->alt_seq);
    free(c->positions);
    memset(c, 0, sizeof(*c));
}

static void free_calls(CallVec *v) {
    for (size_t i = 0; i < v->n; ++i) free_call(&v->data[i]);
    free(v->data);
    v->data = NULL;
    v->n = v->m = 0;
}

static void print_usage(FILE *out) {
    fprintf(out,
        "usage: phase_mnv -r ref.fa [options] input.vcf|input.bcf\n"
        "\n"
        "Build minimal merged MNV/complex records from phased variants in one sample.\n"
        "\n"
        "required:\n"
        "  -r, --reference FILE   Indexed or indexable FASTA reference\n"
        "\n"
        "options:\n"
        "  -s, --sample NAME      Sample to read (default: first sample)\n"
        "  -o, --output FILE      Output VCF path (default: stdout; plain text)\n"
        "  -g, --max-gap N        Allow up to N unchanged reference bases between\n"
        "                        phased variants when building one merged call (default: 0)\n"
        "      --min-vars N       Minimum source variants per emitted call (default: 2)\n"
        "      --min-snvs N       Alias for --min-vars\n"
        "      --no-ref-check     Do not fail when VCF REF differs from FASTA\n"
        "      --no-header        Suppress VCF header\n"
        "  -q, --quiet            Suppress summary on stderr\n"
        "  -h, --help             Show this help\n"
        "\n"
        "Notes:\n"
        "  * Only phased diploid GT (e.g. 0|1, 1|0, 1|1) is used. Unphased\n"
        "    genotypes and symbolic/breakend/non-DNA alleles are skipped.\n"
        "  * FORMAT/PS is honored when present; variants are only merged within the\n"
        "    same phase set. If PS is absent, the phase separator and proximity\n"
        "    define the merge block.\n"
        "  * With the default --max-gap 0, only adjacent phased variants are\n"
        "    merged. Pure SNV blocks are TYPE=MNV; blocks containing indels are\n"
        "    TYPE=COMPLEX.\n"
    );
}

static int64_t parse_i64(const char *s, const char *name) {
    errno = 0;
    char *end = NULL;
    long long v = strtoll(s, &end, 10);
    if (errno || end == s || *end != '\0') {
        fprintf(stderr, "error: invalid %s: %s\n", name, s);
        exit(EXIT_FAILURE);
    }
    return (int64_t)v;
}

static Config parse_args(int argc, char **argv) {
    Config cfg;
    memset(&cfg, 0, sizeof(cfg));
    cfg.max_gap = 0;
    cfg.min_variants = 2;

    static const struct option long_opts[] = {
        {"reference",    required_argument, 0, 'r'},
        {"sample",       required_argument, 0, 's'},
        {"output",       required_argument, 0, 'o'},
        {"max-gap",      required_argument, 0, 'g'},
        {"min-snvs",     required_argument, 0,  1 },
        {"min-vars",     required_argument, 0,  1 },
        {"no-ref-check", no_argument,       0,  2 },
        {"no-header",    no_argument,       0,  3 },
        {"quiet",        no_argument,       0, 'q'},
        {"help",         no_argument,       0, 'h'},
        {0, 0, 0, 0}
    };

    int c;
    while ((c = getopt_long(argc, argv, "r:s:o:g:qh", long_opts, NULL)) != -1) {
        switch (c) {
            case 'r': cfg.fasta_path = optarg; break;
            case 's': cfg.sample_name = optarg; break;
            case 'o': cfg.output_path = optarg; break;
            case 'g': cfg.max_gap = parse_i64(optarg, "--max-gap"); break;
            case 'q': cfg.quiet = true; break;
            case 'h': print_usage(stdout); exit(EXIT_SUCCESS);
            case 1: cfg.min_variants = (int)parse_i64(optarg, "--min-vars"); break;
            case 2: cfg.no_ref_check = true; break;
            case 3: cfg.no_header = true; break;
            default: print_usage(stderr); exit(EXIT_FAILURE);
        }
    }

    if (!cfg.fasta_path) {
        print_usage(stderr);
        die("--reference is required");
    }
    if (cfg.max_gap < 0) die("--max-gap must be >= 0");
    if (cfg.min_variants < 2) die("--min-vars must be >= 2");
    if (optind + 1 != argc) {
        print_usage(stderr);
        die("exactly one input VCF/BCF is required");
    }
    cfg.input_path = argv[optind];
    return cfg;
}

static int resolve_sample_index(const bcf_hdr_t *hdr, const char *sample_name) {
    int nsamples = bcf_hdr_nsamples(hdr);
    if (nsamples <= 0) die("input has no samples; phased GT is required");
    if (!sample_name) return 0;
    for (int i = 0; i < nsamples; ++i) {
        if (strcmp(hdr->samples[i], sample_name) == 0) return i;
    }
    fprintf(stderr, "error: sample '%s' not found. Available samples:", sample_name);
    for (int i = 0; i < nsamples; ++i) fprintf(stderr, " %s", hdr->samples[i]);
    fputc('\n', stderr);
    exit(EXIT_FAILURE);
}

static int64_t get_sample_ps(const bcf_hdr_t *hdr, bcf1_t *rec, int sample_idx,
                             int32_t **ps_arr, int *nps_arr) {
    int nps = bcf_get_format_int32((bcf_hdr_t *)hdr, rec, "PS", ps_arr, nps_arr);
    if (nps <= 0) return PS_MISSING;
    int nsamples = bcf_hdr_nsamples(hdr);
    if (sample_idx >= nps) return PS_MISSING;
    // FORMAT/PS is normally Number=1, so nps == nsamples. If a writer encoded
    // more than one integer per sample, use the first value for this sample.
    int values_per_sample = nps / nsamples;
    if (values_per_sample < 1) values_per_sample = 1;
    int32_t value = (*ps_arr)[sample_idx * values_per_sample];
    if (value == bcf_int32_missing || value == bcf_int32_vector_end) return PS_MISSING;
    return (int64_t)value;
}

static void read_observations(const Config *cfg, bcf_hdr_t **out_hdr, int *out_sample_idx,
                              ObsVec *obs, Stats *st) {
    htsFile *in = hts_open(cfg->input_path, "r");
    if (!in) {
        fprintf(stderr, "error: cannot open input '%s'\n", cfg->input_path);
        exit(EXIT_FAILURE);
    }

    bcf_hdr_t *hdr = bcf_hdr_read(in);
    if (!hdr) die("failed to read VCF/BCF header");
    int sample_idx = resolve_sample_index(hdr, cfg->sample_name);

    bcf1_t *rec = bcf_init();
    int32_t *gt_arr = NULL;
    int ngt_arr = 0;
    int32_t *ps_arr = NULL;
    int nps_arr = 0;
    int nsamples = bcf_hdr_nsamples(hdr);

    while (bcf_read(in, hdr, rec) == 0) {
        st->records++;
        bcf_unpack(rec, BCF_UN_STR | BCF_UN_FMT);

        int ngt = bcf_get_genotypes(hdr, rec, &gt_arr, &ngt_arr);
        if (ngt <= 0) {
            st->skipped_no_gt++;
            continue;
        }
        int ploidy = ngt / nsamples;
        if (ploidy < 2) {
            st->skipped_not_diploid++;
            continue;
        }
        int32_t *gt = gt_arr + sample_idx * ploidy;

        int actual_ploidy = 0;
        for (int k = 0; k < ploidy; ++k) {
            if (gt[k] == bcf_int32_vector_end) break;
            actual_ploidy++;
        }
        if (actual_ploidy != 2) {
            st->skipped_not_diploid++;
            continue;
        }
        if (bcf_gt_is_missing(gt[0]) || bcf_gt_is_missing(gt[1])) {
            st->skipped_missing_gt++;
            continue;
        }
        if (!bcf_gt_is_phased(gt[1])) {
            st->skipped_unphased++;
            continue;
        }
        st->phased_records++;

        if (rec->n_allele < 1 || !is_plain_dna_allele(rec->d.allele[0])) {
            st->skipped_ref++;
            continue;
        }

        int64_t ps = get_sample_ps(hdr, rec, sample_idx, &ps_arr, &nps_arr);
        int64_t pos1 = (int64_t)rec->pos + 1;
        const char *ref_allele = rec->d.allele[0];
        size_t ref_len = strlen(ref_allele);

        for (int hap = 0; hap < 2; ++hap) {
            int allele = bcf_gt_allele(gt[hap]);
            if (allele == 0) {
                st->skipped_ref_allele++;
                continue;
            }
            if (allele < 0 || allele >= rec->n_allele) {
                st->skipped_unsupported_alt++;
                continue;
            }
            const char *alt_allele = rec->d.allele[allele];
            if (!is_plain_dna_allele(alt_allele)) {
                st->skipped_unsupported_alt++;
                continue;
            }
            if (strcasecmp(ref_allele, alt_allele) == 0) {
                st->skipped_unsupported_alt++;
                continue;
            }
            Obs x;
            memset(&x, 0, sizeof(x));
            x.rid = rec->rid;
            x.hap = hap;
            x.ps = ps;
            x.pos = pos1;
            x.end = pos1 + (int64_t)ref_len - 1;
            x.ref = xstrdup_upper(ref_allele);
            x.alt = xstrdup_upper(alt_allele);
            x.is_snv = is_snv_allele_pair(x.ref, x.alt);
            obs_push(obs, x);
            st->observations++;
        }
    }

    free(gt_arr);
    free(ps_arr);
    bcf_destroy(rec);
    hts_close(in);

    *out_hdr = hdr;
    *out_sample_idx = sample_idx;
}

static int cmp_obs(const void *ap, const void *bp) {
    const Obs *a = (const Obs *)ap;
    const Obs *b = (const Obs *)bp;
    if (a->rid != b->rid) return (a->rid < b->rid) ? -1 : 1;
    if (a->hap != b->hap) return a->hap - b->hap;
    if (a->ps != b->ps) return (a->ps < b->ps) ? -1 : 1;
    if (a->pos != b->pos) return (a->pos < b->pos) ? -1 : 1;
    if (a->end != b->end) return (a->end < b->end) ? -1 : 1;
    int c = strcmp(a->alt, b->alt);
    if (c) return c;
    return strcmp(a->ref, b->ref);
}

static int cmp_calls(const void *ap, const void *bp) {
    const MnvCall *a = (const MnvCall *)ap;
    const MnvCall *b = (const MnvCall *)bp;
    if (a->rid != b->rid) return (a->rid < b->rid) ? -1 : 1;
    if (a->start != b->start) return (a->start < b->start) ? -1 : 1;
    if (a->end != b->end) return (a->end < b->end) ? -1 : 1;
    int c = strcmp(a->ref_seq, b->ref_seq);
    if (c) return c;
    c = strcmp(a->alt_seq, b->alt_seq);
    if (c) return c;
    c = strcmp(a->positions, b->positions);
    if (c) return c;
    return 0;
}

static bool can_extend(const Obs *a, const Obs *b, int64_t max_gap) {
    if (a->rid != b->rid || a->hap != b->hap || a->ps != b->ps) return false;
    // Overlapping records on the same haplotype are ambiguous to compose unless
    // first normalized/decomposed upstream. Adjacent records have gap 0.
    if (b->pos <= a->end) return false;
    return (b->pos - a->end - 1) <= max_gap;
}

static char *make_positions_string(const Obs *obs, size_t first, size_t last) {
    size_t cap = (last - first + 1) * 24 + 1;
    char *s = (char *)xmalloc(cap);
    s[0] = '\0';
    size_t used = 0;
    for (size_t i = first; i <= last; ++i) {
        char buf[64];
        int n = snprintf(buf, sizeof(buf), "%s%" PRId64,
                         i == first ? "" : ",", obs[i].pos);
        if (used + (size_t)n + 1 > cap) {
            cap = (cap + (size_t)n + 64) * 2;
            s = (char *)xrealloc(s, cap);
        }
        memcpy(s + used, buf, (size_t)n + 1);
        used += (size_t)n;
    }
    return s;
}

static char *prepend_base(char *s, char base) {
    size_t len = strlen(s);
    char *out = (char *)xmalloc(len + 2);
    out[0] = base;
    memcpy(out + 1, s, len + 1);
    free(s);
    return out;
}

static int fetch_left_base(const faidx_t *fai, const char *chrom, int64_t new_pos1, char *base) {
    hts_pos_t len = 0;
    char *seq = faidx_fetch_seq64(fai, chrom, (hts_pos_t)(new_pos1 - 1),
                                  (hts_pos_t)(new_pos1 - 1), &len);
    if (!seq || len != 1) {
        free(seq);
        return -1;
    }
    *base = upbase(seq[0]);
    free(seq);
    return 0;
}

// Normalize a biallelic VCF representation in-place using the left-aligned +
// parsimonious rules from:
//
//   Tan A, Abecasis GR, Kang HM. Unified representation of genetic variants.
//   Bioinformatics. 2015;31(13):2202-2204. doi:10.1093/bioinformatics/btv112
//
// Algorithm: repeatedly right-trim common suffixes with left-extension when an
// allele becomes empty, then left-trim common prefixes while all alleles remain
// non-empty. This makes our output left-aligned and parsimonious without an
// external `vt normalize`/`bcftools norm` pass.
static int normalize_biallelic(const faidx_t *fai, const char *chrom,
                               int64_t *pos1, char **ref, char **alt) {
    bool changed = true;
    while (changed) {
        changed = false;
        size_t rlen = strlen(*ref);
        size_t alen = strlen(*alt);

        if (rlen == 0 || alen == 0) {
            if (*pos1 <= 1) {
                fprintf(stderr, "error: cannot left-extend variant at beginning of contig %s\n", chrom);
                return -1;
            }
            char base = 'N';
            int64_t new_pos1 = *pos1 - 1;
            if (fetch_left_base(fai, chrom, new_pos1, &base) != 0) {
                fprintf(stderr, "error: failed to fetch left-extension base %s:%" PRId64 "\n",
                        chrom, new_pos1);
                return -1;
            }
            *ref = prepend_base(*ref, base);
            *alt = prepend_base(*alt, base);
            *pos1 = new_pos1;
            changed = true;
            continue;
        }

        bool would_empty_at_contig_start = (*pos1 == 1 && (rlen == 1 || alen == 1));
        if (!would_empty_at_contig_start && upbase((*ref)[rlen - 1]) == upbase((*alt)[alen - 1])) {
            (*ref)[rlen - 1] = '\0';
            (*alt)[alen - 1] = '\0';
            changed = true;
        }
    }

    while (strlen(*ref) > 1 && strlen(*alt) > 1 && upbase((*ref)[0]) == upbase((*alt)[0])) {
        size_t rlen = strlen(*ref);
        size_t alen = strlen(*alt);
        memmove(*ref, *ref + 1, rlen); // includes NUL
        memmove(*alt, *alt + 1, alen);
        ++(*pos1);
    }

    return 0;
}

static int append_to_string(char **buf, size_t *len, size_t *cap, const char *src, size_t n) {
    if (*len + n + 1 > *cap) {
        while (*len + n + 1 > *cap) *cap = *cap ? (*cap * 2) : 64;
        *buf = (char *)xrealloc(*buf, *cap);
    }
    memcpy(*buf + *len, src, n);
    *len += n;
    (*buf)[*len] = '\0';
    return 0;
}

static int add_call_from_block(const Config *cfg, const bcf_hdr_t *hdr, const faidx_t *fai,
                               const Obs *obs, size_t first, size_t last, CallVec *calls) {
    const Obs *a = &obs[first];
    const Obs *b = &obs[last];
    const char *chrom = bcf_hdr_id2name(hdr, a->rid);
    if (!chrom) {
        fprintf(stderr, "error: input record has invalid contig id %d\n", a->rid);
        return -1;
    }

    int64_t span_end = b->end;
    for (size_t i = first; i <= last; ++i) {
        if (obs[i].end > span_end) span_end = obs[i].end;
    }

    hts_pos_t ref_len = 0;
    char *ref = faidx_fetch_seq64(fai, chrom, (hts_pos_t)(a->pos - 1),
                                  (hts_pos_t)(span_end - 1), &ref_len);
    int64_t expected_len = span_end - a->pos + 1;
    if (!ref || ref_len != (hts_pos_t)expected_len) {
        fprintf(stderr,
                "error: failed to fetch reference %s:%" PRId64 "-%" PRId64
                " from '%s' (got length %" PRId64 ")\n",
                chrom, a->pos, span_end, cfg->fasta_path, (int64_t)ref_len);
        free(ref);
        return -1;
    }

    for (int64_t i = 0; i < expected_len; ++i) ref[i] = upbase(ref[i]);

    size_t alt_cap = (size_t)expected_len + 64;
    size_t alt_len = 0;
    char *alt = (char *)xmalloc(alt_cap);
    alt[0] = '\0';

    int64_t cursor = a->pos;
    int nsnps = 0;
    int nvars = (int)(last - first + 1);
    bool all_snvs = true;

    for (size_t i = first; i <= last; ++i) {
        if (obs[i].pos < cursor) {
            fprintf(stderr,
                    "error: overlapping phased records at %s:%" PRId64
                    " on one haplotype; normalize/decompose input first\n",
                    chrom, obs[i].pos);
            free(ref);
            free(alt);
            return -1;
        }

        int64_t copy_start_off = cursor - a->pos;
        int64_t copy_len = obs[i].pos - cursor;
        if (copy_len > 0) {
            append_to_string(&alt, &alt_len, &alt_cap,
                             ref + copy_start_off, (size_t)copy_len);
        }

        int64_t off = obs[i].pos - a->pos;
        size_t ref_allele_len = strlen(obs[i].ref);
        if (off < 0 || off + (int64_t)ref_allele_len > expected_len) {
            fprintf(stderr, "error: internal offset bug while building merged call\n");
            free(ref);
            free(alt);
            return -1;
        }
        if (!cfg->no_ref_check && strncasecmp(ref + off, obs[i].ref, ref_allele_len) != 0) {
            char *fasta_piece = (char *)xmalloc(ref_allele_len + 1);
            memcpy(fasta_piece, ref + off, ref_allele_len);
            fasta_piece[ref_allele_len] = '\0';
            fprintf(stderr,
                    "error: REF/FASTA mismatch at %s:%" PRId64
                    " (VCF REF=%s FASTA=%s). Use --no-ref-check to ignore.\n",
                    chrom, obs[i].pos, obs[i].ref, fasta_piece);
            free(fasta_piece);
            free(ref);
            free(alt);
            return -1;
        }

        append_to_string(&alt, &alt_len, &alt_cap, obs[i].alt, strlen(obs[i].alt));
        cursor = obs[i].end + 1;
        if (obs[i].is_snv) nsnps++;
        else all_snvs = false;
    }

    int64_t tail_off = cursor - a->pos;
    int64_t tail_len = span_end - cursor + 1;
    if (tail_len > 0) {
        append_to_string(&alt, &alt_len, &alt_cap, ref + tail_off, (size_t)tail_len);
    }

    if (strcmp(ref, alt) == 0) {
        free(ref);
        free(alt);
        return 0;
    }

    int64_t norm_pos = a->pos;
    if (normalize_biallelic(fai, chrom, &norm_pos, &ref, &alt) != 0) {
        free(ref);
        free(alt);
        return -1;
    }

    MnvCall call;
    memset(&call, 0, sizeof(call));
    call.rid = a->rid;
    call.start = norm_pos;
    call.end = norm_pos + (int64_t)strlen(ref) - 1;
    call.ref_seq = ref;
    call.alt_seq = alt;
    call.positions = make_positions_string(obs, first, last);
    call.nvars = nvars;
    call.nsnps = nsnps;
    call.type = all_snvs ? "MNV" : "COMPLEX";
    call.hap_mask = 1 << a->hap;
    call.ps = a->ps;
    call_push(calls, call);
    return 0;
}

static int build_calls(const Config *cfg, const bcf_hdr_t *hdr, const faidx_t *fai,
                       ObsVec *obs, CallVec *calls) {
    if (obs->n == 0) return 0;
    qsort(obs->data, obs->n, sizeof(obs->data[0]), cmp_obs);

    size_t i = 0;
    while (i < obs->n) {
        size_t j = i;
        while (j + 1 < obs->n && can_extend(&obs->data[j], &obs->data[j + 1], cfg->max_gap)) {
            ++j;
        }
        if ((int)(j - i + 1) >= cfg->min_variants) {
            if (add_call_from_block(cfg, hdr, fai, obs->data, i, j, calls) != 0) return -1;
        }
        i = j + 1;
    }
    return 0;
}

static void merge_duplicate_calls(CallVec *calls) {
    if (calls->n == 0) return;
    qsort(calls->data, calls->n, sizeof(calls->data[0]), cmp_calls);

    size_t out = 0;
    for (size_t i = 0; i < calls->n; ++i) {
        if (out > 0 && calls->data[out - 1].rid == calls->data[i].rid &&
            calls->data[out - 1].start == calls->data[i].start &&
            calls->data[out - 1].end == calls->data[i].end &&
            strcmp(calls->data[out - 1].ref_seq, calls->data[i].ref_seq) == 0 &&
            strcmp(calls->data[out - 1].alt_seq, calls->data[i].alt_seq) == 0 &&
            strcmp(calls->data[out - 1].positions, calls->data[i].positions) == 0) {
            calls->data[out - 1].hap_mask |= calls->data[i].hap_mask;
            if (calls->data[out - 1].ps != calls->data[i].ps) calls->data[out - 1].ps = PS_MISSING;
            free_call(&calls->data[i]);
        } else {
            if (out != i) calls->data[out] = calls->data[i];
            ++out;
        }
    }
    calls->n = out;
}

static const char *gt_for_mask(int mask) {
    switch (mask & 3) {
        case 1: return "1|0";
        case 2: return "0|1";
        case 3: return "1|1";
        default: return "./.";
    }
}

static FILE *open_output(const char *path) {
    if (!path || strcmp(path, "-") == 0) return stdout;
    FILE *fp = fopen(path, "w");
    if (!fp) {
        fprintf(stderr, "error: cannot open output '%s': %s\n", path, strerror(errno));
        exit(EXIT_FAILURE);
    }
    return fp;
}

static void write_header(FILE *out, const Config *cfg, const bcf_hdr_t *hdr, int sample_idx) {
    fprintf(out, "##fileformat=VCFv4.3\n");
    fprintf(out, "##source=phase_mnv\n");
    fprintf(out, "##phase_mnv_normalization=Tan2015_left_aligned_parsimonious\n");
    fprintf(out, "##phase_mnv_normalization_citation=Tan_A_Abecasis_GR_Kang_HM_Bioinformatics_2015_31_13_2202_2204_doi_10.1093/bioinformatics/btv112\n");
    fprintf(out, "##phase_mnv_input=%s\n", cfg->input_path);
    fprintf(out, "##reference=%s\n", cfg->fasta_path);
    fprintf(out, "##INFO=<ID=TYPE,Number=1,Type=String,Description=\"Merged call type: MNV for pure SNV blocks, COMPLEX when indels are included\">\n");
    fprintf(out, "##INFO=<ID=NVAR,Number=1,Type=Integer,Description=\"Number of phased source variants merged into this call\">\n");
    fprintf(out, "##INFO=<ID=NSNPS,Number=1,Type=Integer,Description=\"Number of source SNVs in this merged call\">\n");
    fprintf(out, "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End coordinate of merged reference span\">\n");
    fprintf(out, "##INFO=<ID=SOURCE_POS,Number=.,Type=Integer,Description=\"Original source variant positions merged into this call\">\n");
    fprintf(out, "##INFO=<ID=HAPS,Number=.,Type=Integer,Description=\"One-based phased haplotypes carrying this merged call\">\n");
    fprintf(out, "##INFO=<ID=PS,Number=1,Type=Integer,Description=\"Phase set shared by merged variants, when present in input FORMAT/PS\">\n");
    fprintf(out, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Phased genotype for the constructed call in the selected sample\">\n");
    fprintf(out, "##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"Phase set for the constructed call, or missing if absent/ambiguous\">\n");

    int nseq = 0;
    const char **seqnames = bcf_hdr_seqnames(hdr, &nseq);
    if (seqnames) {
        for (int i = 0; i < nseq; ++i) {
            if (seqnames[i]) fprintf(out, "##contig=<ID=%s>\n", seqnames[i]);
        }
        free(seqnames);
    }

    fprintf(out, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n",
            hdr->samples[sample_idx]);
}

static void write_calls(FILE *out, const bcf_hdr_t *hdr, const CallVec *calls, Stats *st) {
    for (size_t i = 0; i < calls->n; ++i) {
        const MnvCall *c = &calls->data[i];
        const char *chrom = bcf_hdr_id2name(hdr, c->rid);
        if (!chrom) chrom = ".";
        const char *gt = gt_for_mask(c->hap_mask);
        char psbuf[64];
        const char *ps = ".";
        if (c->ps != PS_MISSING) {
            snprintf(psbuf, sizeof(psbuf), "%" PRId64, c->ps);
            ps = psbuf;
        }
        char haps[8];
        if ((c->hap_mask & 3) == 3) strcpy(haps, "1,2");
        else if (c->hap_mask & 1) strcpy(haps, "1");
        else if (c->hap_mask & 2) strcpy(haps, "2");
        else strcpy(haps, ".");

        fprintf(out,
                "%s\t%" PRId64 "\t.\t%s\t%s\t.\tPASS\t"
                "TYPE=%s;NVAR=%d;NSNPS=%d;END=%" PRId64 ";SOURCE_POS=%s;HAPS=%s",
                chrom, c->start, c->ref_seq, c->alt_seq, c->type,
                c->nvars, c->nsnps, c->end, c->positions, haps);
        if (c->ps != PS_MISSING) fprintf(out, ";PS=%" PRId64, c->ps);
        fprintf(out, "\tGT:PS\t%s:%s\n", gt, ps);
        st->emitted++;
    }
}

static void print_summary(const Config *cfg, const Stats *st, const char *sample) {
    if (cfg->quiet) return;
    fprintf(stderr,
            "phase_mnv: sample=%s records=%" PRIu64 " phased_records=%" PRIu64
            " haplotype_variant_observations=%" PRIu64 " emitted_calls=%" PRIu64 "\n",
            sample, st->records, st->phased_records, st->observations, st->emitted);
    fprintf(stderr,
            "phase_mnv: skipped no_gt=%" PRIu64 " non_diploid=%" PRIu64
            " missing_gt=%" PRIu64 " unphased=%" PRIu64 " unsupported_ref=%" PRIu64
            " unsupported_alt=%" PRIu64 " ref_hap_alleles=%" PRIu64 "\n",
            st->skipped_no_gt, st->skipped_not_diploid, st->skipped_missing_gt,
            st->skipped_unphased, st->skipped_ref, st->skipped_unsupported_alt,
            st->skipped_ref_allele);
}

int main(int argc, char **argv) {
    Config cfg = parse_args(argc, argv);

    ObsVec obs = {0};
    CallVec calls = {0};
    Stats st;
    memset(&st, 0, sizeof(st));

    bcf_hdr_t *hdr = NULL;
    int sample_idx = -1;
    read_observations(&cfg, &hdr, &sample_idx, &obs, &st);

    faidx_t *fai = fai_load(cfg.fasta_path);
    if (!fai) {
        fprintf(stderr, "error: cannot load or create FASTA index for '%s'\n", cfg.fasta_path);
        bcf_hdr_destroy(hdr);
        free_obs(&obs);
        return EXIT_FAILURE;
    }

    if (build_calls(&cfg, hdr, fai, &obs, &calls) != 0) {
        fai_destroy(fai);
        bcf_hdr_destroy(hdr);
        free_obs(&obs);
        free_calls(&calls);
        return EXIT_FAILURE;
    }
    merge_duplicate_calls(&calls);

    FILE *out = open_output(cfg.output_path);
    if (!cfg.no_header) write_header(out, &cfg, hdr, sample_idx);
    write_calls(out, hdr, &calls, &st);
    if (out != stdout) fclose(out);

    print_summary(&cfg, &st, hdr->samples[sample_idx]);

    fai_destroy(fai);
    bcf_hdr_destroy(hdr);
    free_obs(&obs);
    free_calls(&calls);
    return EXIT_SUCCESS;
}

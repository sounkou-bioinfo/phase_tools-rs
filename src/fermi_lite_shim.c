// Minimal fermi-lite helpers normally supplied by bseq.c.
//
// These two functions are copied from fermi-lite's MIT-licensed bseq.c
// (see vendor/fermi-lite/LICENSE.txt). We do not compile bseq.c itself because
// this crate feeds reads from Rust and does not need zlib-backed FASTA/FASTQ
// input helpers.

void seq_reverse(int l, unsigned char *s)
{
    int i;
    for (i = 0; i < l>>1; ++i) {
        int tmp = s[l-1-i];
        s[l-1-i] = s[i]; s[i] = tmp;
    }
}

void seq_revcomp6(int l, unsigned char *s)
{
    int i;
    for (i = 0; i < l>>1; ++i) {
        int tmp = s[l-1-i];
        tmp = (tmp >= 1 && tmp <= 4)? 5 - tmp : tmp;
        s[l-1-i] = (s[i] >= 1 && s[i] <= 4)? 5 - s[i] : s[i];
        s[i] = tmp;
    }
    if (l&1) s[i] = (s[i] >= 1 && s[i] <= 4)? 5 - s[i] : s[i];
}

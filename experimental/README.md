Experiments
-----------

### Decycling-set based minimizer orders

I adapted the code from https://github.com/OrensteinLab/DecyclingSetBasedMinimizerOrder and described in the paper
*"Efficient minimizer orders for large values of k using minimum decycling sets"* by Pellow et al.

My adaptation is in `include/util.hpp`.

These new orders can be used as drop-in replacements for the random minimizer order currently used in SSHash. If an order induces lower density, then SSHash maintains less minimizers and longer (on average) super-kmers, thus improving space usage.

A small example for E. Coli (seed is 1 by default):

	./sshash build -i ../data/unitigs_stitched/with_weights/ecoli_sakai.ust.k31.fa.gz -k 31 -m 13

gives the following results.

- Random order: `4.77532 [bits/kmer]` (`488490` minimizers and `525171` super-kmers)
- Decycling: `4.30659 [bits/kmer]` (`397127` minimizers and `435437` super-kmers)
- Double-Decycling: `4.24312 [bits/kmer]` (`384739` minimizers and `423233` super-kmers)

So, yes, we see a good improvement in space. Beware that, however, query time becomes problematic since checking the membership of a m-mer a decycling set is theoretically O(m), hence we compute a minimizer for a k-mer in O(m(k-m+1)).

Using seed `1234567890`:

	./sshash build -i ../data/unitigs_stitched/with_weights/ecoli_sakai.ust.k31.fa.gz -k 31 -m 13 -s 1234567890

gives the following results.

- Random order: `4.83112 [bits/kmer]` (`489522` minimizers and `526050` super-kmers)
- Decycling: `4.66423 [bits/kmer]` (`461927` minimizers and `504579` super-kmers)
- Double-Decycling: `4.66423 [bits/kmer]` (`461927` minimizers and `504579` super-kmers)

So, a minor improvement. It seems the decycling-set based method is very susceptible to the seed, which is rather strange.

For a whole human genome:

	./sshash build -k 31 -m 21 -i /data2/DNA/lphash_datasets/human.k31.unitigs.fa.ust.fa.gz --verbose

we've got the following results.

- Random order: `8.69547 [bits/kmer]` (`389482029` minimizers and `428428536` super-kmers)
- Decycling: `8.4507 [bits/kmer]` (`372475238` minimizers and `412169244` super-kmers)
- Double-Decycling: `8.38363 [bits/kmer]` (`369875186` minimizers and `408404053` super-kmers)


For a whole yeast genome:

	./sshash build -k 31 -m 15 -i /data2/DNA/lphash_datasets/yeast.k31.unitigs.fa.ust.fa.gz --verbose -s 1234567890

we have the following numbers.

- Random order: `5.26637 [bits/kmer]` (`1257282` minimizers and `1288199` super-kmers)
- Decycling: `5.37412 [bits/kmer]` (`1326252` minimizers and `1358643` super-kmers)
- Double-Decycling: `5.37412 [bits/kmer]` (`1326252` minimizers and `1358643` super-kmers)

Using seed `0`:

	./sshash build -k 31 -m 15 -i /data2/DNA/lphash_datasets/yeast.k31.unitigs.fa.ust.fa.gz --verbose -s 0

we've got the following results instead.

- Random order: `5.19991 [bits/kmer]` (`1256895` minimizers and `1288072` super-kmers)
- Decycling: `4.74185 [bits/kmer]` (`1062025` minimizers and `1094399` super-kmers)
- Double-Decycling: `4.56576 [bits/kmer]` (`989783` minimizers and `1022506` super-kmers)

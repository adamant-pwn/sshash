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

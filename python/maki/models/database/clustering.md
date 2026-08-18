# Clustered Database

The shape of the colour array relates to the dataset:

  1. It's length is determined by the number of unique $k$-mers in the dataset.
  2. It's width is determined by the number of unique $k$-mer distributions throughout the dataset.
  3. It also takes space to store the distributions.

We need to derive:

  1. The expected number of unique $k$-mers in a collection of $N$ genomes at pairwise mutation rate $r$.
  2. The expected number of unique distributions of a $k$-mer in a collection of genomes at pairwise mutation rate $r$.

## Number of unique $k$-mers

Say any pair of genomes has an average mutation rate $r$ and average length $L$. Consider an alignment of $N$ such genomes. When considering the number of $k$-mers conserved in this alignment, it can be treated as a pair of sequences with a higher effective mutation rate.

Each position in one of the genomes has probability $r$ of being mutated. For the first non-seed genome, the set of mutation positions is found by taking $rL$ draws from the set $\left[L\right]$ with replacement. The next genome takes a further $rL$ draws. For a set of $N$ genomes, the mutated positions are found by taking $(N-1)rL$ draws from $\left[L\right]$ with replacement. The effective mutation rate therefore becomes

$$
\begin{align}
  \hat{r}_N &= 1 - \hat{i}_N, \\
  \hat{i}_N &= \left( 1 - 1/L \right)^{\left( N-1 \right)rL},
\end{align}
$$

where $\hat{i}_N$ denotes the effective identity rate.

To estimate the total number of unique $k$-mers

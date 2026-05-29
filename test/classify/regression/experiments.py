
from generate_templates import generate_templates
from simulate_community import ReferenceToken, mutate, simulate_community
from kmer_count import kmer_counts, match_kmers
from inference import formulate_matrices, estimate_abundances


def main():
    k = 7
    
    templates = list(generate_templates(1, 1, min_len=50, max_len=100))
    for i in range(2):
        templates.append(mutate(templates[0], 10, i+1))
        
    community = [
      ReferenceToken(0, 0, 6),
      ReferenceToken(1, 0, 5),
    ]
    reads = simulate_community(templates, community, 10, read_len=10, insert_size=10)
    
    counts = kmer_counts(reads, k)
    matches = match_kmers(counts, templates)
    
    X, y = formulate_matrices(matches, counts)
    results = estimate_abundances(X, y)
    print(results)

        
if __name__ == "__main__":
    main()

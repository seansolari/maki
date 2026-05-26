import random

DNA = "ACGT"
COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}

def random_seq(length):
    return "".join(random.choice(DNA) for _ in range(length))

def generate_templates(num_templates: int, seed: int, min_len: int = 500, max_len = 1000):
    random.seed(seed)
    
    for _ in range(num_templates):
        length = random.randint(min_len, max_len)
        yield random_seq(length)


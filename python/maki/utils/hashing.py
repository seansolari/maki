import hashlib


def hash_strings(values):
    h = hashlib.sha256()
    for v in sorted(values):
        h.update(v.encode())
    return h.hexdigest()

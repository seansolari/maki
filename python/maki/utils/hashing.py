import hashlib
from typing import Iterable


def hash_strings(values: Iterable[str]):
    h = hashlib.sha256()
    for v in sorted(values):
        h.update(v.encode())
    return h.hexdigest()

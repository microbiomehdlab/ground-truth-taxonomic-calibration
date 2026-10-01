"""Stable nested biological allocations; no profile or model execution."""
from __future__ import annotations
import hashlib
import itertools
import random


def allocations(samples, cohort, background, repetitions=1000, seed='da_v1_20261001',
                sizes=(5,10,15,20), exhaustive=False):
    samples=sorted(samples)
    if len(samples)!=len(set(samples)) or not samples:
        raise ValueError('Unique nonempty biological pool required')
    pool_hash=hashlib.sha256('\n'.join(samples).encode()).hexdigest()
    namespace='|'.join((seed,cohort,background,'subset' if exhaustive else 'full',pool_hash))
    if exhaustive:
        if len(samples)!=10: raise ValueError('Enumeration requires exactly ten people')
        for index,cases in enumerate(itertools.combinations(samples,5)):
            controls=tuple(s for s in samples if s not in cases)
            yield dict(allocation_id=f'subset_{index:03d}',n=5,cases=cases,controls=controls,
                       pool_hash=pool_hash,seed_namespace=namespace)
        return
    feasible=sorted(n for n in sizes if 2*n<=len(samples))
    if not feasible: return
    maximum=max(feasible)
    if repetitions<1: raise ValueError('Positive repetition count required')
    for index in range(repetitions):
        integer=int.from_bytes(hashlib.sha256((namespace+'|'+str(index)).encode()).digest(),'big')
        shuffled=samples.copy()
        random.Random(integer).shuffle(shuffled)
        for n in feasible:
            yield dict(allocation_id=f'full_{index:04d}',n=n,
                       cases=tuple(shuffled[:n]),controls=tuple(shuffled[maximum:maximum+n]),
                       pool_hash=pool_hash,seed_namespace=namespace)

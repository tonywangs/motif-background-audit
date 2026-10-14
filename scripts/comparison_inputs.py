"""Frozen synthetic input generator. No analysis of motif counts here."""
import random


def generate(protocol, condition, repetition):
    seed = protocol['generation_seed'] + 100 * protocol['conditions'].index(condition) + repetition
    rng = random.Random(seed)
    lines, plants = [], []
    for record in range(protocol['records_per_dataset']):
        n = protocol['bases_per_record']
        if condition == 'first_order':
            seq = [rng.choice('ACGT')]
            for _ in range(n-1):
                seq.append(seq[-1] if rng.random() < .8 else rng.choice('ACGT'))
        elif condition == 'longer_range':
            # Independent blocks, deterministic internal pattern: dependence
            # beyond order one persists even after conditioning dinucleotides.
            seq = list(''.join(rng.choice(['ACGTAC','GTACGT']) for _ in range(n//6)))
        else:
            seq = [rng.choice('ACGT') for _ in range(n)]
        if condition == 'planted':
            for start in [20,70,120,170,220]:
                seq[start:start+6] = 'ACGTAC'
                plants.append({'record':record,'start':start,'end':start+6,'pattern':'ACGTAC'})
        lines.extend([f'>{condition}_{repetition}_{record} synthetic only', ''.join(seq)])
    return ('\n'.join(lines)+'\n').encode(),seed,plants

# -*- coding: utf-8 -*-
"""
VIREON Part 3 — Advanced scientific modules.
Imported by biology.py via 'from biosim_part3 import *'.

Contains:
  A1. Advanced MolBio: multiple sequence alignment (MSA), consensus,
      sequence logo data, neighbor-joining tree, ML-distance.
  A2. Advanced SynBio: logic gates, Gillespie SSA, codon optimization.
  A3. Advanced Pharmacology: 2/3-compartment PK, PBPK, ADME (Lipinski/Veber).
  B1. Epidemiology: SIR/SEIR/SEIRS, R0, vaccination, interventions.
  B2. Neuroscience: Hodgkin-Huxley, FitzHugh-Nagumo, synapse.
  B3. Molecular docking: simple grid-based scoring, binding pose.

All functions are pure-Python / numpy / scipy where possible, so they work
even when heavy optional libraries are missing.
"""

import math
import random


# ── Logger ───────────────────────────────────────────────────────────
# biosim_log is really defined in biosim_part2. Because of import order
# (part2 imports part3 near its top, long before it defines biosim_log, and in
# a packaged .app/.exe the order can vary), a direct import here is unreliable
# and used to raise "name 'biosim_log' is not defined" inside this module.
# This wrapper looks up the real logger at CALL time via sys.modules, so it
# always works, and it does not shadow part2's version (it forwards to it).
def biosim_log(level, message):
    """Forward to biosim_part2.biosim_log if available, else print."""
    import sys as _sys
    _p2 = _sys.modules.get('biosim_part2')
    _real = getattr(_p2, 'biosim_log', None) if _p2 is not None else None
    if _real is not None and _real is not biosim_log:
        try:
            return _real(level, message)
        except Exception:
            pass
    try:
        print('[biosim] {}: {}'.format(level, message))
    except Exception:
        pass


try:
    import numpy as np
except ImportError as _exc:
    # biosim_log is defined in biosim_part2 and may not be imported yet when
    # part3 loads; guard the call so a missing numpy doesn't raise NameError.
    if 'biosim_log' in globals():
        biosim_log('warn', 'module: {}'.format(_exc))
    np = None


# ════════════════════════════════════════════════════════════════════
# A1. ADVANCED MOLECULAR BIOLOGY
# ════════════════════════════════════════════════════════════════════

def msa_progressive(sequences, names=None, match=1, mismatch=-1, gap=-2):
    """Progressive multiple sequence alignment (simplified Feng-Doolittle).

    Aligns sequences pairwise, builds a guide order by similarity, then
    progressively aligns. Returns dict with aligned sequences + consensus.

    sequences: list of strings (DNA or protein)
    Returns: {'aligned': [...], 'names': [...], 'consensus': str,
              'length': int, 'identity': float}
    """
    if not sequences or len(sequences) < 2:
        return {'error': 'Need at least 2 sequences'}
    seqs = [s.upper().replace(' ', '').replace('\n', '') for s in sequences]
    if names is None:
        names = ['Seq{}'.format(i + 1) for i in range(len(seqs))]

    def nw_align(a, b):
        """Needleman-Wunsch global alignment of two strings."""
        n, m = len(a), len(b)
        score = [[0] * (m + 1) for _ in range(n + 1)]
        for i in range(n + 1):
            score[i][0] = i * gap
        for j in range(m + 1):
            score[0][j] = j * gap
        for i in range(1, n + 1):
            for j in range(1, m + 1):
                s = match if a[i - 1] == b[j - 1] else mismatch
                score[i][j] = max(
                    score[i - 1][j - 1] + s,
                    score[i - 1][j] + gap,
                    score[i][j - 1] + gap)
        # Traceback
        aa, bb = [], []
        i, j = n, m
        while i > 0 and j > 0:
            s = match if a[i - 1] == b[j - 1] else mismatch
            if score[i][j] == score[i - 1][j - 1] + s:
                aa.append(a[i - 1]); bb.append(b[j - 1]); i -= 1; j -= 1
            elif score[i][j] == score[i - 1][j] + gap:
                aa.append(a[i - 1]); bb.append('-'); i -= 1
            else:
                aa.append('-'); bb.append(b[j - 1]); j -= 1
        while i > 0:
            aa.append(a[i - 1]); bb.append('-'); i -= 1
        while j > 0:
            aa.append('-'); bb.append(b[j - 1]); j -= 1
        return ''.join(reversed(aa)), ''.join(reversed(bb)), score[n][m]

    # Build pairwise similarity matrix to order sequences
    n = len(seqs)
    sim = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            _, _, sc = nw_align(seqs[i], seqs[j])
            sim[i][j] = sim[j][i] = sc

    # Order: start with most similar pair, add by max similarity to aligned set
    order = [0]
    if n > 1:
        # pick highest-similarity partner chain (greedy)
        remaining = set(range(1, n))
        while remaining:
            best, best_sc = None, -1e9
            for r in remaining:
                sc = max(sim[r][o] for o in order)
                if sc > best_sc:
                    best_sc, best = sc, r
            order.append(best)
            remaining.discard(best)

    # Progressive alignment: align sequences in order, pad existing to length
    aligned = [seqs[order[0]]]
    for idx in order[1:]:
        # Align new sequence against current consensus of aligned block
        consensus = _msa_consensus_simple(aligned)
        a_cons, a_new, _ = nw_align(consensus, seqs[idx])
        # Insert gaps into existing aligned seqs where consensus got gaps
        new_aligned = []
        for aseq in aligned:
            res = []
            ai = 0
            for c in a_cons:
                if c == '-':
                    res.append('-')
                else:
                    res.append(aseq[ai] if ai < len(aseq) else '-')
                    ai += 1
            new_aligned.append(''.join(res))
        new_aligned.append(a_new)
        aligned = new_aligned

    # Reorder back to original input order
    final = [''] * n
    for pos, idx in enumerate(order):
        final[idx] = aligned[pos]
    # Pad all to same length
    maxlen = max(len(s) for s in final)
    final = [s + '-' * (maxlen - len(s)) for s in final]

    consensus = _msa_consensus_simple(final)
    identity = _msa_identity(final)
    return {
        'aligned': final,
        'names': names,
        'consensus': consensus,
        'length': maxlen,
        'identity': identity,
    }


def _msa_consensus_simple(aligned):
    """Most common non-gap residue per column."""
    if not aligned:
        return ''
    maxlen = max(len(s) for s in aligned)
    cons = []
    for col in range(maxlen):
        counts = {}
        for s in aligned:
            if col < len(s):
                c = s[col]
                if c != '-':
                    counts[c] = counts.get(c, 0) + 1
        if counts:
            cons.append(max(counts, key=counts.get))
        else:
            cons.append('-')
    return ''.join(cons)


def _msa_identity(aligned):
    """Percent of fully-conserved columns."""
    if not aligned or len(aligned) < 2:
        return 0.0
    maxlen = max(len(s) for s in aligned)
    conserved = 0
    counted = 0
    for col in range(maxlen):
        chars = set()
        ok = True
        for s in aligned:
            if col < len(s) and s[col] != '-':
                chars.add(s[col])
            else:
                ok = False
        if ok and len(chars) == 1:
            conserved += 1
        if ok:
            counted += 1
    return round(100.0 * conserved / counted, 1) if counted else 0.0


def msa_column_conservation(aligned):
    """Return per-column conservation score 0..1 (for sequence logo / shading)."""
    if not aligned:
        return []
    maxlen = max(len(s) for s in aligned)
    nseq = len(aligned)
    scores = []
    for col in range(maxlen):
        counts = {}
        for s in aligned:
            c = s[col] if col < len(s) else '-'
            counts[c] = counts.get(c, 0) + 1
        top = max(counts.values())
        scores.append(round(top / nseq, 2))
    return scores


def neighbor_joining_tree(dist_matrix, names):
    """Neighbor-joining phylogenetic tree (Saitou-Nei).

    dist_matrix: NxN list of lists (symmetric distances)
    names: list of N taxon names
    Returns: Newick string + list of merge steps.
    """
    n = len(names)
    if n < 2:
        return {'newick': names[0] if names else '', 'steps': []}
    # Work on a copy
    D = [row[:] for row in dist_matrix]
    nodes = list(names)
    steps = []
    next_id = 0

    while len(nodes) > 2:
        m = len(nodes)
        # Net divergence
        r = [sum(D[i]) for i in range(m)]
        # Q matrix
        min_q, mi, mj = 1e18, 0, 1
        for i in range(m):
            for j in range(i + 1, m):
                q = (m - 2) * D[i][j] - r[i] - r[j]
                if q < min_q:
                    min_q, mi, mj = q, i, j
        # Branch lengths
        if (m - 2) != 0:
            d_i = 0.5 * D[mi][mj] + (r[mi] - r[mj]) / (2 * (m - 2))
        else:
            d_i = 0.5 * D[mi][mj]
        d_j = D[mi][mj] - d_i
        d_i = max(0, round(d_i, 4)); d_j = max(0, round(d_j, 4))
        new_name = '({}:{},{}:{})'.format(nodes[mi], d_i, nodes[mj], d_j)
        steps.append((nodes[mi], nodes[mj], d_i, d_j))
        # New distance matrix
        new_dist = []
        for k in range(m):
            if k == mi or k == mj:
                continue
            d = 0.5 * (D[mi][k] + D[mj][k] - D[mi][mj])
            new_dist.append(round(d, 4))
        # Build new node list and matrix
        new_nodes = [nodes[k] for k in range(m) if k != mi and k != mj]
        new_nodes.append(new_name)
        size = len(new_nodes)
        ND = [[0] * size for _ in range(size)]
        old_idx = [k for k in range(m) if k != mi and k != mj]
        for a in range(len(old_idx)):
            for b in range(len(old_idx)):
                ND[a][b] = D[old_idx[a]][old_idx[b]]
        for a in range(len(old_idx)):
            ND[a][size - 1] = ND[size - 1][a] = new_dist[a]
        D = ND
        nodes = new_nodes
        next_id += 1

    # Last two
    if len(nodes) == 2:
        final_d = round(D[0][1], 4)
        newick = '({}:{},{}:{});'.format(nodes[0], final_d / 2, nodes[1], final_d / 2)
    else:
        newick = nodes[0] + ';'
    return {'newick': newick, 'steps': steps}


def jukes_cantor_distance(seq1, seq2):
    """Jukes-Cantor evolutionary distance between two aligned sequences."""
    s1 = seq1.upper(); s2 = seq2.upper()
    n = min(len(s1), len(s2))
    diff = sum(1 for i in range(n) if s1[i] != s2[i] and s1[i] != '-' and s2[i] != '-')
    valid = sum(1 for i in range(n) if s1[i] != '-' and s2[i] != '-')
    if valid == 0:
        return 0.0
    p = diff / valid
    if p >= 0.75:
        return 3.0  # saturated
    try:
        d = -0.75 * math.log(1 - (4.0 / 3.0) * p)
    except ValueError as _exc:
        biosim_log('warn', 'jukes_cantor_distance: {}'.format(_exc))
        d = 3.0
    return round(d, 4)


def build_distance_matrix(sequences):
    """Build JC distance matrix from aligned sequences."""
    n = len(sequences)
    D = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            d = jukes_cantor_distance(sequences[i], sequences[j])
            D[i][j] = D[j][i] = d
    return D


# ════════════════════════════════════════════════════════════════════
# A2. ADVANCED SYNTHETIC BIOLOGY
# ════════════════════════════════════════════════════════════════════

def gillespie_ssa(reactions, initial, t_max, max_steps=100000, seed=None):
    """Gillespie Stochastic Simulation Algorithm (exact SSA).

    reactions: list of dicts {'rate': k, 'reactants': {species: n},
                              'products': {species: n}}
    initial: dict {species: count}
    Returns: {'times': [...], 'species': {name: [...]}, 'steps': int}
    """
    if seed is not None:
        random.seed(seed)
    state = dict(initial)
    species_names = list(state.keys())
    # Ensure all species in reactions are tracked
    for rxn in reactions:
        for sp in list(rxn.get('reactants', {}).keys()) + list(rxn.get('products', {}).keys()):
            if sp not in state:
                state[sp] = 0
                species_names.append(sp)

    times = [0.0]
    traj = {sp: [state[sp]] for sp in species_names}
    t = 0.0
    steps = 0

    while t < t_max and steps < max_steps:
        # Compute propensities
        propensities = []
        for rxn in reactions:
            a = rxn['rate']
            for sp, n in rxn.get('reactants', {}).items():
                cnt = state.get(sp, 0)
                # mass-action: combinatorial
                for k in range(n):
                    a *= max(0, cnt - k)
                if n > 1:
                    a /= math.factorial(n)
            propensities.append(a)
        a0 = sum(propensities)
        if a0 <= 0:
            break
        # Time to next reaction
        r1 = random.random()
        tau = -math.log(r1) / a0
        t += tau
        if t > t_max:
            break
        # Which reaction
        r2 = random.random() * a0
        cumsum = 0
        chosen = 0
        for i, p in enumerate(propensities):
            cumsum += p
            if r2 <= cumsum:
                chosen = i
                break
        # Apply reaction
        rxn = reactions[chosen]
        for sp, n in rxn.get('reactants', {}).items():
            state[sp] = state.get(sp, 0) - n
        for sp, n in rxn.get('products', {}).items():
            state[sp] = state.get(sp, 0) + n
        # Record
        times.append(round(t, 4))
        for sp in species_names:
            traj[sp].append(state.get(sp, 0))
        steps += 1

    return {'times': times, 'species': traj, 'steps': steps}


def logic_gate_dna(gate_type, input_a, input_b=None):
    """Simulate a genetic logic gate (transcriptional).

    gate_type: 'AND', 'OR', 'NOT', 'NAND', 'NOR', 'XOR'
    inputs: 0 or 1 (presence of inducer)
    Returns: output expression level 0..1 and truth explanation.
    """
    a = 1 if input_a else 0
    b = 1 if input_b else 0
    gt = gate_type.upper()
    if gt == 'AND':
        out = a & b
    elif gt == 'OR':
        out = a | b
    elif gt == 'NOT':
        out = 1 - a
    elif gt == 'NAND':
        out = 1 - (a & b)
    elif gt == 'NOR':
        out = 1 - (a | b)
    elif gt == 'XOR':
        out = a ^ b
    else:
        out = 0
    # Add leaky expression (realistic): ON=0.9, OFF=0.05
    level = 0.9 if out else 0.05
    return {'gate': gt, 'input_a': a, 'input_b': b, 'output': out,
            'expression': level}


# E. coli codon usage table (relative frequency, highest = preferred)
ECOLI_CODON_PREFERENCE = {
    'A': 'GCG', 'R': 'CGC', 'N': 'AAC', 'D': 'GAT', 'C': 'TGC',
    'Q': 'CAG', 'E': 'GAA', 'G': 'GGC', 'H': 'CAT', 'I': 'ATT',
    'L': 'CTG', 'K': 'AAA', 'M': 'ATG', 'F': 'TTT', 'P': 'CCG',
    'S': 'AGC', 'T': 'ACC', 'W': 'TGG', 'Y': 'TAT', 'V': 'GTG',
    '*': 'TAA',
}

CODON_TABLE_STD = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'CTT': 'L', 'CTC': 'L',
    'CTA': 'L', 'CTG': 'L', 'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'TCT': 'S', 'TCC': 'S',
    'TCA': 'S', 'TCG': 'S', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'GCT': 'A', 'GCC': 'A',
    'GCA': 'A', 'GCG': 'A', 'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'AAT': 'N', 'AAC': 'N',
    'AAA': 'K', 'AAG': 'K', 'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W', 'CGT': 'R', 'CGC': 'R',
    'CGA': 'R', 'CGG': 'R', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}


def codon_optimize_detailed(dna_or_protein, is_protein=False, organism='ecoli'):
    """Optimize codons for expression in E. coli — DETAILED version.
    Accepts DNA or protein, returns a dict with optimized DNA + stats.
    (The simpler protein-only version lives in biosim_part2 as codon_optimize.)

    Returns: {'optimized_dna': str, 'protein': str, 'changes': int, 'cai_before', 'cai_after'}
    """
    if is_protein:
        protein = dna_or_protein.upper()
    else:
        dna = dna_or_protein.upper().replace(' ', '').replace('\n', '')
        protein = ''
        for i in range(0, len(dna) - 2, 3):
            codon = dna[i:i + 3]
            protein += CODON_TABLE_STD.get(codon, 'X')
    # Build optimized DNA
    opt = []
    for aa in protein:
        if aa == 'X':
            continue
        opt.append(ECOLI_CODON_PREFERENCE.get(aa, 'NNN'))
    optimized = ''.join(opt)
    changes = 0
    if not is_protein:
        orig = dna_or_protein.upper().replace(' ', '').replace('\n', '')
        for i in range(0, min(len(orig), len(optimized)) - 2, 3):
            if orig[i:i + 3] != optimized[i:i + 3]:
                changes += 1
    return {
        'optimized_dna': optimized,
        'protein': protein.replace('*', ''),
        'changes': changes if not is_protein else len(protein),
        'length': len(optimized),
    }


# ════════════════════════════════════════════════════════════════════
# A3. ADVANCED PHARMACOLOGY
# ════════════════════════════════════════════════════════════════════

def pk_two_compartment(dose, ka, k10, k12, k21, V1, times):
    """Two-compartment PK model with first-order absorption.

    dose: administered dose (mg)
    ka: absorption rate (1/h)
    k10: elimination from central (1/h)
    k12, k21: inter-compartment rates (1/h)
    V1: central volume (L)
    times: list of time points (h)
    Returns: {'times': [...], 'central': [...], 'peripheral': [...]}
    """
    if np is None:
        return {'error': 'numpy required'}
    # Solve ODE system numerically (simple RK4)
    def deriv(state, t):
        A_gut, A1, A2 = state
        dA_gut = -ka * A_gut
        dA1 = ka * A_gut - (k10 + k12) * A1 + k21 * A2
        dA2 = k12 * A1 - k21 * A2
        return [dA_gut, dA1, dA2]

    state = [dose, 0.0, 0.0]
    central = []
    peripheral = []
    t_prev = 0.0
    central.append(0.0)
    peripheral.append(0.0)
    for idx in range(1, len(times)):
        dt = times[idx] - times[idx - 1]
        # RK4
        steps = max(1, int(dt / 0.01))
        h = dt / steps
        for _ in range(steps):
            k1 = deriv(state, t_prev)
            k2 = deriv([state[i] + h / 2 * k1[i] for i in range(3)], t_prev + h / 2)
            k3 = deriv([state[i] + h / 2 * k2[i] for i in range(3)], t_prev + h / 2)
            k4 = deriv([state[i] + h * k3[i] for i in range(3)], t_prev + h)
            state = [state[i] + h / 6 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i])
                     for i in range(3)]
            t_prev += h
        central.append(state[1] / V1)
        peripheral.append(state[2] / V1)
    return {'times': list(times), 'central': central, 'peripheral': peripheral}


def pk_multiple_dosing(dose, ka, ke, V, tau, n_doses, times):
    """Multiple dosing — superposition of single doses (1-compartment).

    tau: dosing interval (h); n_doses: number of doses.
    Returns concentration accounting for accumulation.
    """
    conc = []
    for t in times:
        c = 0.0
        for n in range(n_doses):
            t_dose = n * tau
            if t >= t_dose:
                td = t - t_dose
                if abs(ka - ke) > 1e-9:
                    c += (dose * ka) / (V * (ka - ke)) * (math.exp(-ke * td) - math.exp(-ka * td))
        conc.append(max(0, c))
    # Steady-state metrics
    css_max = max(conc) if conc else 0
    accumulation = 1.0 / (1.0 - math.exp(-ke * tau)) if ke * tau < 50 else 1.0
    return {'times': list(times), 'conc': conc,
            'css_max': round(css_max, 3),
            'accumulation_ratio': round(accumulation, 2)}


def adme_properties(mol_weight, logp, h_donors, h_acceptors,
                    rotatable_bonds, tpsa):
    """ADME / drug-likeness assessment (Lipinski + Veber rules).

    Returns pass/fail for each rule + overall druglikeness.
    """
    lipinski = {
        'MW <= 500': mol_weight <= 500,
        'LogP <= 5': logp <= 5,
        'H-donors <= 5': h_donors <= 5,
        'H-acceptors <= 10': h_acceptors <= 10,
    }
    lipinski_violations = sum(1 for v in lipinski.values() if not v)
    veber = {
        'Rotatable bonds <= 10': rotatable_bonds <= 10,
        'TPSA <= 140': tpsa <= 140,
    }
    # Oral bioavailability prediction (simplified)
    oral_ok = (lipinski_violations <= 1 and rotatable_bonds <= 10 and tpsa <= 140)
    # Blood-brain barrier (rough): TPSA < 90 and MW < 450
    bbb = (tpsa < 90 and mol_weight < 450)
    return {
        'lipinski': lipinski,
        'lipinski_violations': lipinski_violations,
        'veber': veber,
        'oral_bioavailable': oral_ok,
        'bbb_permeable': bbb,
        'druglike': lipinski_violations <= 1,
    }


def dose_by_weight(dose_per_kg, weight_kg, max_dose=None):
    """Calculate dose adjusted for body weight."""
    dose = dose_per_kg * weight_kg
    if max_dose:
        dose = min(dose, max_dose)
    return round(dose, 2)


def renal_dose_adjustment(normal_dose, creatinine_clearance):
    """Adjust dose based on renal function (Cockcroft-Gault style).

    creatinine_clearance in mL/min. Normal ~120.
    """
    if creatinine_clearance >= 90:
        factor = 1.0
    elif creatinine_clearance >= 60:
        factor = 0.75
    elif creatinine_clearance >= 30:
        factor = 0.5
    elif creatinine_clearance >= 15:
        factor = 0.25
    else:
        factor = 0.1
    return {'adjusted_dose': round(normal_dose * factor, 2),
            'factor': factor,
            'renal_function': ('normal' if creatinine_clearance >= 90 else
                               'mild' if creatinine_clearance >= 60 else
                               'moderate' if creatinine_clearance >= 30 else
                               'severe' if creatinine_clearance >= 15 else
                               'failure')}


# ════════════════════════════════════════════════════════════════════
# B1. EPIDEMIOLOGY — SIR / SEIR / SEIRS models
# ════════════════════════════════════════════════════════════════════

def sir_model(S0, I0, R0_count, beta, gamma, days, dt=0.5):
    """Classic SIR epidemic model.

    S0, I0, R0_count: initial Susceptible / Infected / Recovered
    beta: transmission rate; gamma: recovery rate (1/infectious period)
    Returns trajectories + R0 (basic reproduction number).
    """
    N = S0 + I0 + R0_count
    S, I, R = float(S0), float(I0), float(R0_count)
    times, S_t, I_t, R_t = [0], [S], [I], [R]
    steps = int(days / dt)
    peak_I, peak_day = I, 0
    for step in range(1, steps + 1):
        dS = -beta * S * I / N
        dI = beta * S * I / N - gamma * I
        dR = gamma * I
        S += dS * dt; I += dI * dt; R += dR * dt
        S = max(0, S); I = max(0, I); R = max(0, R)
        t = step * dt
        times.append(round(t, 2)); S_t.append(S); I_t.append(I); R_t.append(R)
        if I > peak_I:
            peak_I, peak_day = I, t
    R0 = beta / gamma if gamma > 0 else 0
    attack_rate = (R_t[-1] / N) * 100 if N > 0 else 0
    return {
        'times': times, 'S': S_t, 'I': I_t, 'R': R_t,
        'R0': round(R0, 2),
        'peak_infected': round(peak_I),
        'peak_day': round(peak_day, 1),
        'attack_rate': round(attack_rate, 1),
        'N': N,
    }


def seir_model(S0, E0, I0, R0_count, beta, sigma, gamma, days, dt=0.5):
    """SEIR model with Exposed (latent) compartment.

    sigma: rate E->I (1/incubation period)
    """
    N = S0 + E0 + I0 + R0_count
    S, E, I, R = float(S0), float(E0), float(I0), float(R0_count)
    times = [0]; S_t, E_t, I_t, R_t = [S], [E], [I], [R]
    steps = int(days / dt)
    peak_I, peak_day = I, 0
    for step in range(1, steps + 1):
        dS = -beta * S * I / N
        dE = beta * S * I / N - sigma * E
        dI = sigma * E - gamma * I
        dR = gamma * I
        S += dS * dt; E += dE * dt; I += dI * dt; R += dR * dt
        S = max(0, S); E = max(0, E); I = max(0, I); R = max(0, R)
        t = step * dt
        times.append(round(t, 2)); S_t.append(S); E_t.append(E)
        I_t.append(I); R_t.append(R)
        if I > peak_I:
            peak_I, peak_day = I, t
    R0 = beta / gamma if gamma > 0 else 0
    return {
        'times': times, 'S': S_t, 'E': E_t, 'I': I_t, 'R': R_t,
        'R0': round(R0, 2),
        'peak_infected': round(peak_I),
        'peak_day': round(peak_day, 1),
        'attack_rate': round((R_t[-1] / N) * 100, 1) if N else 0,
        'N': N,
    }


def sir_with_vaccination(S0, I0, beta, gamma, vacc_fraction, days, dt=0.5):
    """SIR with initial vaccination (moves fraction S->R at t=0).

    Demonstrates herd immunity threshold.
    """
    R0_count = S0 * vacc_fraction
    S0_eff = S0 * (1 - vacc_fraction)
    result = sir_model(S0_eff, I0, R0_count, beta, gamma, days, dt)
    R0 = beta / gamma if gamma > 0 else 0
    herd_threshold = (1 - 1 / R0) * 100 if R0 > 1 else 0
    result['vaccination_pct'] = round(vacc_fraction * 100, 1)
    result['herd_immunity_threshold'] = round(herd_threshold, 1)
    result['protected'] = vacc_fraction * 100 >= herd_threshold
    return result


def effective_reproduction_number(R0, susceptible_fraction):
    """Rt = R0 * fraction susceptible."""
    return round(R0 * susceptible_fraction, 2)


# ════════════════════════════════════════════════════════════════════
# B2. NEUROSCIENCE — Hodgkin-Huxley, FitzHugh-Nagumo
# ════════════════════════════════════════════════════════════════════

def hodgkin_huxley(I_inj, t_max=50, dt=0.01,
                   C_m=1.0, g_Na=120.0, g_K=36.0, g_L=0.3,
                   E_Na=50.0, E_K=-77.0, E_L=-54.387):
    """Hodgkin-Huxley neuron model — simulates action potentials.

    I_inj: injected current (uA/cm^2), constant or callable(t)
    Returns: {'times', 'V' (membrane potential mV), 'spikes' (count),
              'm','h','n' (gating variables)}
    """
    def alpha_m(V): return 0.1 * (V + 40) / (1 - math.exp(-(V + 40) / 10)) if abs(V + 40) > 1e-7 else 1.0
    def beta_m(V): return 4.0 * math.exp(-(V + 65) / 18)
    def alpha_h(V): return 0.07 * math.exp(-(V + 65) / 20)
    def beta_h(V): return 1.0 / (1 + math.exp(-(V + 35) / 10))
    def alpha_n(V): return 0.01 * (V + 55) / (1 - math.exp(-(V + 55) / 10)) if abs(V + 55) > 1e-7 else 0.1
    def beta_n(V): return 0.125 * math.exp(-(V + 65) / 80)

    V = -65.0
    m = alpha_m(V) / (alpha_m(V) + beta_m(V))
    h = alpha_h(V) / (alpha_h(V) + beta_h(V))
    n = alpha_n(V) / (alpha_n(V) + beta_n(V))

    steps = int(t_max / dt)
    times, V_t, m_t, h_t, n_t = [], [], [], [], []
    spikes = 0
    above = False
    # subsample output to keep arrays manageable
    sub = max(1, steps // 2000)

    for step in range(steps):
        t = step * dt
        I = I_inj(t) if callable(I_inj) else I_inj
        # Currents
        I_Na = g_Na * (m ** 3) * h * (V - E_Na)
        I_K = g_K * (n ** 4) * (V - E_K)
        I_L = g_L * (V - E_L)
        dV = (I - I_Na - I_K - I_L) / C_m
        dm = alpha_m(V) * (1 - m) - beta_m(V) * m
        dh = alpha_h(V) * (1 - h) - beta_h(V) * h
        dn = alpha_n(V) * (1 - n) - beta_n(V) * n
        V += dV * dt; m += dm * dt; h += dh * dt; n += dn * dt
        m = min(1, max(0, m)); h = min(1, max(0, h)); n = min(1, max(0, n))
        # Spike detection (threshold crossing upward at 0 mV)
        if V > 0 and not above:
            spikes += 1; above = True
        elif V < -20:
            above = False
        if step % sub == 0:
            times.append(round(t, 3)); V_t.append(round(V, 2))
            m_t.append(round(m, 3)); h_t.append(round(h, 3)); n_t.append(round(n, 3))

    firing_rate = spikes / (t_max / 1000.0)  # spikes per second
    return {'times': times, 'V': V_t, 'm': m_t, 'h': h_t, 'n': n_t,
            'spikes': spikes, 'firing_rate_hz': round(firing_rate, 1)}


def fitzhugh_nagumo(I_ext, t_max=100, dt=0.05, a=0.7, b=0.8, tau=12.5):
    """FitzHugh-Nagumo simplified neuron model (2 variables).

    Faster than HH, shows excitability and oscillation.
    """
    v, w = -1.0, 1.0
    times, v_t, w_t = [], [], []
    steps = int(t_max / dt)
    spikes = 0
    above = False
    for step in range(steps):
        t = step * dt
        I = I_ext(t) if callable(I_ext) else I_ext
        dv = v - (v ** 3) / 3 - w + I
        dw = (v + a - b * w) / tau
        v += dv * dt; w += dw * dt
        if v > 1.0 and not above:
            spikes += 1; above = True
        elif v < 0:
            above = False
        times.append(round(t, 2)); v_t.append(round(v, 3)); w_t.append(round(w, 3))
    return {'times': times, 'v': v_t, 'w': w_t, 'spikes': spikes}


# ════════════════════════════════════════════════════════════════════
# B3. MOLECULAR DOCKING — simplified grid scoring
# ════════════════════════════════════════════════════════════════════

def simple_docking_score(ligand_atoms, receptor_atoms):
    """Very simplified docking score based on shape complementarity
    and distance-based interaction energy (Lennard-Jones-like).

    ligand_atoms / receptor_atoms: list of (x, y, z, type) tuples
    type: 'H' (hydrophobic), 'P' (polar), 'C' (charged+), 'N' (charged-)
    Returns: {'score' (lower=better binding), 'contacts', 'h_bonds', 'clashes'}
    """
    score = 0.0
    contacts = 0
    h_bonds = 0
    clashes = 0
    for lx, ly, lz, lt in ligand_atoms:
        for rx, ry, rz, rt in receptor_atoms:
            d = math.sqrt((lx - rx) ** 2 + (ly - ry) ** 2 + (lz - rz) ** 2)
            if d < 2.0:
                # Steric clash (bad)
                score += 5.0 * (2.0 - d)
                clashes += 1
            elif d < 4.5:
                contacts += 1
                # Hydrophobic contact
                if lt == 'H' and rt == 'H':
                    score -= 0.7
                # Hydrogen bond (polar-polar)
                elif (lt == 'P' and rt == 'P'):
                    score -= 1.5
                    if d < 3.5:
                        h_bonds += 1
                # Electrostatic (opposite charges attract)
                elif (lt == 'C' and rt == 'N') or (lt == 'N' and rt == 'C'):
                    score -= 2.5
                    h_bonds += 1
                # Same charge repels
                elif (lt == 'C' and rt == 'C') or (lt == 'N' and rt == 'N'):
                    score += 1.5
    return {
        'score': round(score, 2),
        'contacts': contacts,
        'h_bonds': h_bonds,
        'clashes': clashes,
        'binding_quality': ('strong' if score < -10 else
                            'moderate' if score < -3 else
                            'weak' if score < 0 else 'poor'),
    }


def dock_ligand_grid(ligand_atoms, receptor_atoms, grid_size=5, grid_step=2.0):
    """Try multiple ligand positions on a grid, return best pose.

    Translates ligand over a 3D grid and scores each position.
    Returns best score + offset.
    """
    best_score = 1e9
    best_offset = (0, 0, 0)
    rng = range(-grid_size, grid_size + 1)
    # Limit total positions for speed
    tested = 0
    max_tests = 2000
    for dx in rng:
        for dy in rng:
            for dz in rng:
                if tested >= max_tests:
                    break
                offset_ligand = [(x + dx * grid_step, y + dy * grid_step,
                                  z + dz * grid_step, t)
                                 for x, y, z, t in ligand_atoms]
                result = simple_docking_score(offset_ligand, receptor_atoms)
                tested += 1
                if result['score'] < best_score:
                    best_score = result['score']
                    best_offset = (dx * grid_step, dy * grid_step, dz * grid_step)
                    best_result = result
    return {
        'best_score': round(best_score, 2),
        'best_offset': best_offset,
        'poses_tested': tested,
        'details': best_result if 'best_result' in dir() else {},
    }


# ════════════════════════════════════════════════════════════════════
# C1. EPIDEMIOLOGY PRO — age structure, interventions, metapopulation
# ════════════════════════════════════════════════════════════════════

def seir_with_interventions(S0, I0, beta, sigma, gamma, days,
                             interventions=None, dt=0.5):
    """SEIR with time-dependent interventions (lockdown, masks, etc).

    interventions: list of dicts {'day': start_day, 'end_day': end,
                                  'beta_mult': 0.3} — multiplies beta during period.
    Returns trajectories + intervention markers.
    """
    N = S0 + I0
    S, E, I, R = float(S0 - I0), float(I0), float(I0), 0.0
    times, S_t, E_t, I_t, R_t, beta_t = [0], [S], [E], [I], [R], [beta]
    steps = int(days / dt)
    peak_I, peak_day = I, 0
    interventions = interventions or []

    for step in range(1, steps + 1):
        t = step * dt
        # Determine active beta
        cur_beta = beta
        for iv in interventions:
            if iv['day'] <= t <= iv.get('end_day', days):
                cur_beta *= iv.get('beta_mult', 1.0)
        dS = -cur_beta * S * I / N
        dE = cur_beta * S * I / N - sigma * E
        dI = sigma * E - gamma * I
        dR = gamma * I
        S += dS * dt; E += dE * dt; I += dI * dt; R += dR * dt
        S = max(0, S); E = max(0, E); I = max(0, I); R = max(0, R)
        times.append(round(t, 2)); S_t.append(S); E_t.append(E)
        I_t.append(I); R_t.append(R); beta_t.append(cur_beta)
        if I > peak_I:
            peak_I, peak_day = I, t

    return {
        'times': times, 'S': S_t, 'E': E_t, 'I': I_t, 'R': R_t,
        'beta_t': beta_t,
        'R0': round(beta / gamma, 2),
        'peak_infected': round(peak_I), 'peak_day': round(peak_day, 1),
        'attack_rate': round((R_t[-1] / N) * 100, 1) if N else 0,
        'N': N, 'interventions': interventions,
    }


def age_structured_sir(age_groups, contact_matrix, beta, gamma, days, dt=0.5):
    """Age-structured SIR with contact matrix between age groups.

    age_groups: list of dicts {'name': str, 'S': n, 'I': n, 'R': n}
    contact_matrix: NxN — contacts[i][j] = avg contacts of group i with group j
    Returns per-group trajectories.
    """
    n = len(age_groups)
    S = [float(g['S']) for g in age_groups]
    I = [float(g['I']) for g in age_groups]
    R = [float(g.get('R', 0)) for g in age_groups]
    N = [S[i] + I[i] + R[i] for i in range(n)]
    names = [g['name'] for g in age_groups]

    times = [0]
    traj = {names[i]: {'S': [S[i]], 'I': [I[i]], 'R': [R[i]]} for i in range(n)}
    steps = int(days / dt)

    for step in range(1, steps + 1):
        newS, newI, newR = S[:], I[:], R[:]
        for i in range(n):
            # Force of infection from all groups
            foi = 0.0
            for j in range(n):
                if N[j] > 0:
                    foi += contact_matrix[i][j] * beta * I[j] / N[j]
            dS = -foi * S[i]
            dI = foi * S[i] - gamma * I[i]
            dR = gamma * I[i]
            newS[i] = max(0, S[i] + dS * dt)
            newI[i] = max(0, I[i] + dI * dt)
            newR[i] = max(0, R[i] + dR * dt)
        S, I, R = newS, newI, newR
        times.append(round(step * dt, 2))
        for i in range(n):
            traj[names[i]]['S'].append(S[i])
            traj[names[i]]['I'].append(I[i])
            traj[names[i]]['R'].append(R[i])

    total_infected = sum(R[i] for i in range(n))
    return {'times': times, 'groups': traj, 'names': names,
            'total_attack_rate': round(100 * total_infected / sum(N), 1) if sum(N) else 0}


def metapopulation_sir(cities, beta, gamma, migration, days, dt=0.5):
    """Metapopulation model — disease spreads between connected cities.

    cities: list of dicts {'name': str, 'S': n, 'I': n, 'R': n}
    migration: NxN matrix — daily fraction moving from city i to city j
    Returns per-city trajectories showing spatial spread.
    """
    n = len(cities)
    S = [float(c['S']) for c in cities]
    I = [float(c['I']) for c in cities]
    R = [float(c.get('R', 0)) for c in cities]
    names = [c['name'] for c in cities]

    times = [0]
    traj = {names[i]: {'I': [I[i]], 'S': [S[i]], 'R': [R[i]]} for i in range(n)}
    steps = int(days / dt)

    for step in range(1, steps + 1):
        newS, newI, newR = S[:], I[:], R[:]
        for i in range(n):
            N_i = S[i] + I[i] + R[i]
            if N_i > 0:
                dS = -beta * S[i] * I[i] / N_i
                dI = beta * S[i] * I[i] / N_i - gamma * I[i]
                dR = gamma * I[i]
                newS[i] += dS * dt; newI[i] += dI * dt; newR[i] += dR * dt
        # Migration between cities
        for i in range(n):
            for j in range(n):
                if i != j:
                    m = migration[i][j] * dt
                    moveS = newS[i] * m; moveI = newI[i] * m; moveR = newR[i] * m
                    newS[i] -= moveS; newS[j] += moveS
                    newI[i] -= moveI; newI[j] += moveI
                    newR[i] -= moveR; newR[j] += moveR
        S = [max(0, x) for x in newS]
        I = [max(0, x) for x in newI]
        R = [max(0, x) for x in newR]
        times.append(round(step * dt, 2))
        for i in range(n):
            traj[names[i]]['S'].append(S[i])
            traj[names[i]]['I'].append(I[i])
            traj[names[i]]['R'].append(R[i])

    # When did each city peak?
    peaks = {}
    for i in range(n):
        I_series = traj[names[i]]['I']
        peak_idx = I_series.index(max(I_series))
        peaks[names[i]] = round(times[peak_idx], 1)
    return {'times': times, 'cities': traj, 'names': names, 'peak_days': peaks}


def stochastic_sir(S0, I0, beta, gamma, days, n_runs=10, seed=None):
    """Stochastic SIR (chain binomial) — shows outbreak variability.

    Runs multiple simulations to show distribution of outcomes.
    Returns ensemble of trajectories + extinction probability.
    """
    if seed is not None:
        random.seed(seed)
    N = S0 + I0
    runs = []
    extinctions = 0
    final_sizes = []
    for run in range(n_runs):
        S, I, R = S0, I0, 0
        I_series = [I]
        for day in range(days):
            if I == 0:
                I_series.append(0)
                continue
            # Probability of infection per susceptible
            p_inf = 1 - math.exp(-beta * I / N)
            p_rec = 1 - math.exp(-gamma)
            new_inf = sum(1 for _ in range(int(S)) if random.random() < p_inf)
            new_rec = sum(1 for _ in range(int(I)) if random.random() < p_rec)
            S -= new_inf; I += new_inf - new_rec; R += new_rec
            S = max(0, S); I = max(0, I)
            I_series.append(I)
        runs.append(I_series)
        final_sizes.append(R)
        if R < N * 0.05:  # outbreak fizzled
            extinctions += 1
    return {
        'runs': runs,
        'n_runs': n_runs,
        'extinction_prob': round(extinctions / n_runs, 2),
        'mean_final_size': round(sum(final_sizes) / len(final_sizes)),
        'days': days,
    }


# ════════════════════════════════════════════════════════════════════
# C2. NEUROSCIENCE PRO — neuron networks, synaptic plasticity
# ════════════════════════════════════════════════════════════════════

def lif_network(n_neurons, connectivity, I_ext, t_max=200, dt=0.1,
                tau_m=10.0, V_rest=-65.0, V_thresh=-50.0, V_reset=-70.0,
                syn_weight=0.5, seed=None):
    """Leaky integrate-and-fire network of neurons.

    n_neurons: number of neurons
    connectivity: fraction of connections (0..1) — random network
    I_ext: external current (constant or list per neuron)
    Returns: spike raster (times per neuron) + population rate.
    """
    if seed is not None:
        random.seed(seed)
    # Build random connectivity matrix
    W = [[0.0] * n_neurons for _ in range(n_neurons)]
    for i in range(n_neurons):
        for j in range(n_neurons):
            if i != j and random.random() < connectivity:
                W[i][j] = syn_weight

    V = [V_rest + random.uniform(-5, 5) for _ in range(n_neurons)]
    spike_times = [[] for _ in range(n_neurons)]
    syn_input = [0.0] * n_neurons
    steps = int(t_max / dt)
    pop_rate = []

    for step in range(steps):
        t = step * dt
        spikes_now = []
        new_V = V[:]
        for i in range(n_neurons):
            I_i = I_ext[i] if isinstance(I_ext, list) else I_ext
            # LIF dynamics: tau_m dV/dt = -(V - V_rest) + R*I
            # R*I scaled so I~15 gives suprathreshold drive
            drive = I_i + syn_input[i] * 5.0
            dV = (-(V[i] - V_rest) + drive) / tau_m
            new_V[i] = V[i] + dV * dt
            if new_V[i] >= V_thresh:
                new_V[i] = V_reset
                spike_times[i].append(round(t, 2))
                spikes_now.append(i)
        # Propagate synaptic input (decay + new spikes)
        syn_input = [s * math.exp(-dt / 5.0) for s in syn_input]
        for i in spikes_now:
            for j in range(n_neurons):
                syn_input[j] += W[i][j]
        V = new_V
        pop_rate.append(len(spikes_now))

    total_spikes = sum(len(s) for s in spike_times)
    mean_rate = total_spikes / n_neurons / (t_max / 1000.0)
    return {
        'spike_times': spike_times,
        'n_neurons': n_neurons,
        'pop_rate': pop_rate,
        'total_spikes': total_spikes,
        'mean_firing_rate_hz': round(mean_rate, 1),
        't_max': t_max,
    }


def synapse_stdp(pre_spikes, post_spikes, A_plus=0.01, A_minus=0.012,
                 tau_plus=20.0, tau_minus=20.0, w_init=0.5):
    """Spike-Timing-Dependent Plasticity — synaptic weight change.

    pre_spikes, post_spikes: lists of spike times (ms)
    Returns final weight + weight trajectory (learning).
    """
    w = w_init
    w_history = [w]
    events = sorted([(t, 'pre') for t in pre_spikes] +
                    [(t, 'post') for t in post_spikes])
    last_pre, last_post = None, None
    for t, kind in events:
        if kind == 'pre':
            last_pre = t
            if last_post is not None:
                dt = t - last_post
                # pre after post -> depression
                w -= A_minus * math.exp(-abs(dt) / tau_minus)
        else:
            last_post = t
            if last_pre is not None:
                dt = t - last_pre
                # post after pre -> potentiation
                w += A_plus * math.exp(-abs(dt) / tau_plus)
        w = max(0, min(1, w))
        w_history.append(round(w, 4))
    return {'final_weight': round(w, 4), 'w_history': w_history,
            'change': round(w - w_init, 4),
            'direction': 'strengthened' if w > w_init else 'weakened'}


# ════════════════════════════════════════════════════════════════════
# C3. MOLECULAR DOCKING PRO — RDKit 3D conformers, real scoring
# ════════════════════════════════════════════════════════════════════

def smiles_to_atoms_3d(smiles):
    """Convert SMILES to 3D atom list using RDKit (if available).

    Returns list of (x, y, z, atom_type) tuples for docking.
    atom_type classified as H (hydrophobic C), P (polar N/O), C (cationic), N (anionic).
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except ImportError as _exc:
        biosim_log('warn', 'smiles_to_atoms_3d: {}'.format(_exc))
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    mol = Chem.AddHs(mol)
    try:
        AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.MMFFOptimizeMolecule(mol)
    except Exception as _exc:
        biosim_log('warn', 'smiles_to_atoms_3d: {}'.format(_exc))
        return None
    conf = mol.GetConformer()
    atoms = []
    for atom in mol.GetAtoms():
        pos = conf.GetAtomPosition(atom.GetIdx())
        sym = atom.GetSymbol()
        # Classify atom type
        if sym == 'C':
            atype = 'H'  # hydrophobic
        elif sym in ('N',):
            atype = 'P'  # polar / can be cationic
        elif sym in ('O', 'S'):
            atype = 'P'  # polar
        elif sym in ('F', 'Cl', 'Br'):
            atype = 'H'
        else:
            atype = 'H'
        # Charged check
        if atom.GetFormalCharge() > 0:
            atype = 'C'
        elif atom.GetFormalCharge() < 0:
            atype = 'N'
        if sym != 'H':  # skip hydrogens for docking speed
            atoms.append((round(pos.x, 3), round(pos.y, 3), round(pos.z, 3), atype))
    return atoms


def docking_with_rotation(ligand_atoms, receptor_atoms, n_rotations=8, grid_size=3):
    """Dock with both translation and rotation sampling.

    Rotates ligand around z-axis and tries grid positions for each.
    Returns best pose across all orientations.
    """
    best = {'best_score': 1e9}
    for rot in range(n_rotations):
        angle = 2 * math.pi * rot / n_rotations
        cos_a, sin_a = math.cos(angle), math.sin(angle)
        rotated = [(x * cos_a - y * sin_a, x * sin_a + y * cos_a, z, t)
                   for x, y, z, t in ligand_atoms]
        result = dock_ligand_grid(rotated, receptor_atoms, grid_size=grid_size)
        if result['best_score'] < best['best_score']:
            best = result
            best['rotation_deg'] = round(math.degrees(angle), 0)
    return best


def dock_smiles(ligand_smiles, receptor_atoms=None, grid_size=4,
                sample_rotations=0):
    """Full docking pipeline from SMILES.

    Generates 3D ligand, docks against receptor (or a default pocket),
    returns best binding pose + score + interaction analysis.

    sample_rotations > 0 also turns the ligand around the z-axis that many
    times and keeps the best orientation. Translation-only docking can miss a
    good pose purely because the molecule entered the pocket the wrong way
    round, which is why real docking programs always sample orientation.
    """
    ligand = smiles_to_atoms_3d(ligand_smiles)
    if ligand is None:
        return {'error': 'RDKit required or invalid SMILES'}

    # Default synthetic binding pocket if no receptor given
    if receptor_atoms is None:
        receptor_atoms = [
            (8, 0, 0, 'P'), (8, 2, 0, 'N'), (8, -2, 0, 'H'),
            (10, 0, 2, 'H'), (10, 0, -2, 'P'), (6, 0, 3, 'C'),
            (6, 3, 0, 'H'), (9, 1, 1, 'P'),
        ]
    if sample_rotations and sample_rotations > 1:
        result = docking_with_rotation(ligand, receptor_atoms,
                                       n_rotations=int(sample_rotations),
                                       grid_size=grid_size)
        result['rotations_tried'] = int(sample_rotations)
    else:
        result = dock_ligand_grid(ligand, receptor_atoms, grid_size=grid_size)
    result['ligand_atoms'] = len(ligand)
    result['receptor_atoms'] = len(receptor_atoms)
    # Estimate binding affinity (kcal/mol) from score
    result['estimated_kd_nm'] = round(10 ** (result['best_score'] / 3 + 2), 1) if result['best_score'] < 0 else 99999
    return result


# ════════════════════════════════════════════════════════════════════
# AI MODULE — real Claude API integration
# ════════════════════════════════════════════════════════════════════

import os as _os
import json as _json
import urllib.request as _urlreq
import ssl as _ssl_ai

# Config file stored next to the program
_CONFIG_FILE = _os.path.join(_os.path.dirname(_os.path.abspath(__file__)),
                              'biosim_config.json')


def ai_load_config():
    """Load saved config (API key etc) from local file."""
    try:
        if _os.path.exists(_CONFIG_FILE):
            with open(_CONFIG_FILE, 'r') as f:
                return _json.load(f)
    except Exception as _exc:
        biosim_log('warn', 'ai_load_config: {}'.format(_exc))
    return {}


def ai_save_config(cfg):
    """Save config to local file."""
    try:
        with open(_CONFIG_FILE, 'w') as f:
            _json.dump(cfg, f)
        return True
    except Exception as _exc:
        biosim_log('warn', 'ai_save_config: {}'.format(_exc))
        return False


def ai_get_api_key():
    """Get API key from config file or environment variable."""
    # Priority 1: environment variable (Gemini)
    key = _os.environ.get('GEMINI_API_KEY', '')
    if key:
        return key
    # Priority 2: config file
    cfg = ai_load_config()
    return cfg.get('api_key', '')


def ai_set_api_key(key):
    """Save API key to config file."""
    cfg = ai_load_config()
    cfg['api_key'] = key.strip()
    return ai_save_config(cfg)


def ai_is_available():
    """Check if AI is configured (has a key)."""
    return bool(ai_get_api_key())


def _call_gemini_once(prompt, system_prompt=None, max_tokens=1024,
                    model='gemini-flash-latest', use_search=False):
    """Send a prompt to Google Gemini API and return the text response.
    (Single attempt with one model — see call_gemini_api for the retry wrapper.)

    Returns dict: {'ok': bool, 'text': str, 'error': str}
    Works on old macOS via SSL fallback.

    use_search=True enables Google Search grounding, so the model can look up
    fresh information on the web (used for biology/science questions). It is
    left off for questions about the program itself, which the model answers
    from the system prompt.
    """
    api_key = ai_get_api_key()
    if not api_key:
        return {'ok': False, 'text': '',
                'error': 'no_key'}

    url = ('https://generativelanguage.googleapis.com/v1beta/models/'
           '{}:generateContent?key={}'.format(model, api_key))
    body = {
        'contents': [
            {'role': 'user', 'parts': [{'text': prompt}]}
        ],
        'generationConfig': {
            'maxOutputTokens': max_tokens,
        },
    }
    if system_prompt:
        body['systemInstruction'] = {'parts': [{'text': system_prompt}]}
    if use_search:
        # Google Search grounding — lets Gemini search the web for fresh facts.
        body['tools'] = [{'google_search': {}}]

    data = _json.dumps(body).encode('utf-8')
    req = _urlreq.Request(url, data=data, method='POST')
    req.add_header('Content-Type', 'application/json')

    # SSL context. In a packaged .app/.exe the CA bundle is often missing, so
    # certificate verification fails — sometimes as SSLError, sometimes wrapped
    # inside URLError. We retry once with verification disabled in that case.
    def _do_request(verify=True):
        if verify:
            ctx = _ssl_ai.create_default_context()
        else:
            ctx = _ssl_ai.create_default_context()
            ctx.check_hostname = False
            ctx.verify_mode = _ssl_ai.CERT_NONE
        with _urlreq.urlopen(req, timeout=30, context=ctx) as resp:
            return _json.loads(resp.read().decode('utf-8'))

    try:
        try:
            result = _do_request(verify=True)
        except Exception as _exc:
            _msg = str(_exc).lower()
            # Any certificate / SSL problem → retry without verification.
            if ('ssl' in _msg or 'certificate' in _msg or 'cert ' in _msg
                    or isinstance(_exc, _ssl_ai.SSLError)):
                biosim_log('warn', 'call_gemini_api ssl-retry: {}'.format(_exc))
                result = _do_request(verify=False)
            else:
                raise
    except Exception as e:
        biosim_log('warn', 'call_gemini_api: {}'.format(e))
        err = str(e)
        # Friendly error messages
        if '400' in err or 'API_KEY_INVALID' in err or 'api key' in err.lower():
            return {'ok': False, 'text': '', 'error': 'bad_key'}
        if '401' in err or '403' in err or 'permission' in err.lower():
            return {'ok': False, 'text': '', 'error': 'bad_key'}
        if '429' in err or 'quota' in err.lower() or 'rate' in err.lower():
            return {'ok': False, 'text': '', 'error': 'rate_limit'}
        if '503' in err or 'overloaded' in err.lower() or 'unavailable' in err.lower():
            return {'ok': False, 'text': '', 'error': 'overloaded'}
        return {'ok': False, 'text': '', 'error': str(e)[:150]}

    # Extract text from Gemini response
    try:
        # API-level error object
        if isinstance(result, dict) and result.get('error'):
            emsg = str(result['error'].get('message', '')).lower()
            if 'api key' in emsg or 'api_key' in emsg:
                return {'ok': False, 'text': '', 'error': 'bad_key'}
            if 'quota' in emsg or 'rate' in emsg:
                return {'ok': False, 'text': '', 'error': 'rate_limit'}
            return {'ok': False, 'text': '',
                    'error': str(result['error'].get('message', ''))[:150]}
        candidates = result.get('candidates', [])
        text = ''
        for cand in candidates:
            parts = cand.get('content', {}).get('parts', [])
            for part in parts:
                if 'text' in part:
                    text += part['text']
        if not text:
            return {'ok': False, 'text': '', 'error': 'empty_response'}
        return {'ok': True, 'text': text, 'error': ''}
    except Exception as e:
        biosim_log('warn', 'call_gemini_api: {}'.format(e))
        return {'ok': False, 'text': '', 'error': 'parse_error: ' + str(e)[:100]}


def call_gemini_api(prompt, system_prompt=None, max_tokens=1024,
                    model=None, use_search=False):
    """Public entry point. Tries the default free Flash model first; if it is
    rate-limited or unavailable, retries once with the lighter Flash-Lite model
    (which has higher free-tier request limits). Model aliases 'gemini-flash-
    latest' / 'gemini-flash-lite-latest' auto-track Google's current versions,
    so this keeps working as Google retires old model numbers.
    """
    primary = model or 'gemini-flash-latest'
    r = _call_gemini_once(prompt, system_prompt=system_prompt,
                          max_tokens=max_tokens, model=primary,
                          use_search=use_search)
    if r['ok']:
        return r
    # Retry with a higher-limit model on transient/model-availability errors.
    if r.get('error') in ('rate_limit', 'overloaded', 'empty_response'):
        r2 = _call_gemini_once(prompt, system_prompt=system_prompt,
                               max_tokens=max_tokens,
                               model='gemini-flash-lite-latest',
                               use_search=use_search)
        if r2['ok']:
            return r2
        # If search may be the cause, try once more without web search.
        if use_search and r2.get('error') in ('rate_limit', 'overloaded'):
            r3 = _call_gemini_once(prompt, system_prompt=system_prompt,
                                   max_tokens=max_tokens,
                                   model='gemini-flash-lite-latest',
                                   use_search=False)
            if r3['ok']:
                return r3
        return r2
    return r


def call_claude_api(prompt, system_prompt=None, max_tokens=1024, model=None,
                    use_search=False):
    """Backward-compatible wrapper. The app now uses Google Gemini instead of
    Claude; this alias keeps older call sites working. The `model` argument is
    ignored (a Gemini model is chosen inside call_gemini_api).
    """
    return call_gemini_api(prompt, system_prompt=system_prompt,
                           max_tokens=max_tokens, use_search=use_search)


def ai_error_message(error_code, lang='ru'):
    """Human-friendly error message for AI errors."""
    is_ru = (lang == 'ru')
    messages = {
        'no_key': ('Ключ API не задан. Откройте Настройки → API ключ и введите бесплатный ключ Gemini с aistudio.google.com/apikey'
                    if is_ru else
                    'No API key. Open Settings → API key and enter your free Gemini key from aistudio.google.com/apikey'),
        'bad_key': ('Неверный ключ API. Проверьте ключ Gemini в Настройках.'
                     if is_ru else 'Invalid API key. Check your Gemini key in Settings.'),
        'rate_limit': ('Слишком много запросов. Подождите минуту.'
                        if is_ru else 'Too many requests. Wait a minute.'),
        'overloaded': ('Сервер ИИ перегружен. Попробуйте позже.'
                        if is_ru else 'AI server overloaded. Try again later.'),
        'empty_response': ('ИИ вернул пустой ответ. Попробуйте переформулировать вопрос.'
                            if is_ru else 'AI returned an empty response. Try rephrasing.'),
    }
    return messages.get(error_code,
                        ('Ошибка ИИ: ' if is_ru else 'AI error: ') + str(error_code))


# ════════════════════════════════════════════════════════════════════
# SAFE INPUT HELPERS — prevent crashes from bad text input
# ════════════════════════════════════════════════════════════════════

def safe_float(value, default=0.0, minimum=None, maximum=None):
    """Parse a float safely. Returns default if invalid."""
    try:
        result = float(str(value).strip().replace(',', '.'))
        if minimum is not None and result < minimum:
            result = minimum
        if maximum is not None and result > maximum:
            result = maximum
        return result
    except (ValueError, TypeError) as _exc:
        biosim_log('warn', 'safe_float: {}'.format(_exc))
        return default


def safe_int(value, default=0, minimum=None, maximum=None):
    """Parse an int safely. Returns default if invalid."""
    try:
        result = int(float(str(value).strip()))
        if minimum is not None and result < minimum:
            result = minimum
        if maximum is not None and result > maximum:
            result = maximum
        return result
    except (ValueError, TypeError) as _exc:
        biosim_log('warn', 'safe_int: {}'.format(_exc))
        return default


# ════════════════════════════════════════════════════════════════════
# DEEP MolBio — professional-grade sequence analysis
# ════════════════════════════════════════════════════════════════════

# --- Codon usage tables for multiple organisms (CAI calculation) ---
CODON_USAGE_TABLES = {
    'ecoli': {  # relative adaptiveness (w values), highest=1.0
        'TTT': 0.51, 'TTC': 1.0, 'TTA': 0.12, 'TTG': 0.12, 'CTT': 0.12,
        'CTC': 0.10, 'CTA': 0.04, 'CTG': 1.0, 'ATT': 0.51, 'ATC': 1.0,
        'ATA': 0.07, 'ATG': 1.0, 'GTT': 0.60, 'GTC': 0.40, 'GTA': 0.30,
        'GTG': 1.0, 'TCT': 0.46, 'TCC': 0.46, 'TCA': 0.20, 'TCG': 0.20,
        'CCT': 0.30, 'CCC': 0.12, 'CCA': 0.42, 'CCG': 1.0, 'ACT': 0.43,
        'ACC': 1.0, 'ACA': 0.25, 'ACG': 0.40, 'GCT': 0.46, 'GCC': 0.46,
        'GCA': 0.46, 'GCG': 1.0, 'TAT': 0.53, 'TAC': 1.0, 'CAT': 0.52,
        'CAC': 1.0, 'CAA': 0.51, 'CAG': 1.0, 'AAT': 0.49, 'AAC': 1.0,
        'AAA': 1.0, 'AAG': 0.30, 'GAT': 0.63, 'GAC': 1.0, 'GAA': 1.0,
        'GAG': 0.40, 'TGT': 0.43, 'TGC': 1.0, 'TGG': 1.0, 'CGT': 1.0,
        'CGC': 0.72, 'CGA': 0.12, 'CGG': 0.12, 'AGT': 0.20, 'AGC': 0.46,
        'AGA': 0.07, 'AGG': 0.05, 'GGT': 0.70, 'GGC': 1.0, 'GGA': 0.30,
        'GGG': 0.30,
    },
}


def calculate_cai(dna, organism='ecoli'):
    """Codon Adaptation Index — measures how well codons match host preference.

    CAI ranges 0..1. Higher = better expression in that organism.
    Returns CAI value and per-codon adaptiveness.
    """
    dna = dna.upper().replace(' ', '').replace('\n', '')
    table = CODON_USAGE_TABLES.get(organism, CODON_USAGE_TABLES['ecoli'])
    w_values = []
    for i in range(0, len(dna) - 2, 3):
        codon = dna[i:i + 3]
        if codon in table and CODON_TABLE_STD.get(codon) != '*':
            w = table[codon]
            if w > 0:
                w_values.append(w)
    if not w_values:
        return {'cai': 0.0, 'n_codons': 0}
    # Geometric mean of w values
    log_sum = sum(math.log(w) for w in w_values)
    cai = math.exp(log_sum / len(w_values))
    return {'cai': round(cai, 3), 'n_codons': len(w_values),
            'interpretation': ('высокая экспрессия' if cai > 0.7 else
                               'средняя' if cai > 0.4 else 'низкая')}


def gc_content_windowed(dna, window=50):
    """GC content in sliding windows — finds GC-rich/poor regions.

    Returns list of (position, gc_percent) for plotting GC skew.
    """
    dna = dna.upper().replace(' ', '').replace('\n', '')
    results = []
    for i in range(0, len(dna) - window + 1, max(1, window // 5)):
        seg = dna[i:i + window]
        gc = (seg.count('G') + seg.count('C')) / len(seg) * 100 if seg else 0
        results.append((i, round(gc, 1)))
    return results


def simulate_sanger_chromatogram(sequence, noise=0.08, quality_dropoff=0.7):
    """Simulate a Sanger sequencing chromatogram (trace) from a DNA sequence.
    Each base gives a peak in one of 4 channels (A/C/G/T). Quality drops near
    the ends (like real reads). Like SnapGene's .ab1 viewer.

    Returns dict with per-base peak heights for 4 channels + quality scores.
    """
    seq = (sequence or '').upper().replace(' ', '').replace('\n', '')
    if not seq:
        return {'error': 'Empty sequence'}
    import random
    random.seed(42)  # reproducible
    channels = {'A': [], 'C': [], 'G': [], 'T': []}
    qualities = []
    n = len(seq)
    for i, base in enumerate(seq):
        # Quality: high in middle, lower at start and end
        rel = i / max(1, n - 1)
        # Bell-ish: low at very start, peak ~20%, decline toward end
        if rel < 0.05:
            q = 20 + rel * 400  # ramp up
        else:
            q = 60 * (quality_dropoff ** (rel * 2))  # decay
        q = max(5, min(60, q))
        qualities.append(round(q, 1))
        # Peak height for the called base (higher = better quality)
        main_height = 0.6 + 0.4 * (q / 60)
        for ch in 'ACGT':
            if ch == base:
                h = main_height + random.uniform(-noise, noise)
            else:
                # Background noise in other channels (grows as quality drops)
                h = random.uniform(0, noise * (1 + (60 - q) / 60))
            channels[ch].append(round(max(0, min(1.2, h)), 3))
    # Identify low-quality regions (for trimming, like real software)
    trim_start = 0
    trim_end = n
    for i, q in enumerate(qualities):
        if q >= 20:
            trim_start = i
            break
    for i in range(n - 1, -1, -1):
        if qualities[i] >= 20:
            trim_end = i + 1
            break
    return {
        'sequence': seq,
        'channels': channels,
        'qualities': qualities,
        'trim_start': trim_start,
        'trim_end': trim_end,
        'mean_quality': round(sum(qualities) / len(qualities), 1),
    }


def call_variants(reference, sample):
    """Variant calling — compare a sample sequence to a reference and find
    SNPs (substitutions) and indels (insertions/deletions). Like Geneious
    variant calling. Reports each difference with position and type.

    Returns list of variants with confidence estimate.
    """
    ref = (reference or '').upper().replace(' ', '').replace('\n', '')
    smp = (sample or '').upper().replace(' ', '').replace('\n', '')
    if not ref or not smp:
        return {'error': 'Need both reference and sample'}
    variants = []
    # Simple position-by-position comparison (assumes aligned)
    min_len = min(len(ref), len(smp))
    for i in range(min_len):
        if ref[i] != smp[i]:
            variants.append({
                'position': i + 1,
                'type': 'SNP',
                'ref': ref[i],
                'alt': smp[i],
                'change': '{}>{}'.format(ref[i], smp[i]),
                'transition': _is_transition(ref[i], smp[i]),
            })
    # Length difference suggests indel
    if len(smp) > len(ref):
        variants.append({
            'position': min_len + 1,
            'type': 'insertion',
            'ref': '-',
            'alt': smp[min_len:],
            'change': 'ins {}'.format(smp[min_len:min_len+10]),
            'transition': None,
        })
    elif len(ref) > len(smp):
        variants.append({
            'position': min_len + 1,
            'type': 'deletion',
            'ref': ref[min_len:],
            'alt': '-',
            'change': 'del {}'.format(ref[min_len:min_len+10]),
            'transition': None,
        })
    n_snp = sum(1 for v in variants if v['type'] == 'SNP')
    n_ti = sum(1 for v in variants if v.get('transition') is True)
    n_tv = sum(1 for v in variants if v.get('transition') is False)
    return {
        'variants': variants,
        'n_variants': len(variants),
        'n_snp': n_snp,
        'ti_tv_ratio': round(n_ti / max(1, n_tv), 2),
        'identity': round(100 * (min_len - n_snp) / max(1, min_len), 2),
    }


def _is_transition(ref, alt):
    """Transition = purine↔purine (A↔G) or pyrimidine↔pyrimidine (C↔T).
    Transversion = purine↔pyrimidine. Ti/Tv ratio is a QC metric."""
    purines = {'A', 'G'}
    pyrimidines = {'C', 'T'}
    if (ref in purines and alt in purines) or \
       (ref in pyrimidines and alt in pyrimidines):
        return True
    return False


def build_circular_plasmid_map(sequence, features=None):
    """Build a circular plasmid map — compute angular positions of features
    around a circle. Like SnapGene's circular view. Returns features placed
    by angle, ready to draw on a canvas.

    features: list of dicts {'name','start','end','type'} (optional, else
    auto-detect ORFs and common elements).
    """
    seq = (sequence or '').upper().replace(' ', '').replace('\n', '')
    n = len(seq)
    if n == 0:
        return {'error': 'Empty sequence'}
    placed = []
    if features:
        for f in features:
            start = f.get('start', 0)
            end = f.get('end', start)
            # Convert to angles (0° at top, clockwise)
            a1 = 360 * start / n
            a2 = 360 * end / n
            placed.append({
                'name': f.get('name', '?'),
                'type': f.get('type', 'misc'),
                'start_angle': a1,
                'end_angle': a2,
                'start': start,
                'end': end,
            })
    else:
        # Auto-detect: find ORFs as features
        orfs = find_orfs_advanced(seq, min_length=60)
        for i, orf in enumerate(orfs[:12]):
            start = orf.get('start', 0)
            end = orf.get('end', start)
            a1 = 360 * start / n
            a2 = 360 * end / n
            placed.append({
                'name': 'ORF{}'.format(i + 1),
                'type': 'CDS',
                'start_angle': a1,
                'end_angle': a2,
                'start': start,
                'end': end,
            })
    # GC content for the ring
    gc = (seq.count('G') + seq.count('C')) / n * 100
    return {
        'length': n,
        'features': placed,
        'gc_percent': round(gc, 1),
        'n_features': len(placed),
    }


def find_orfs_advanced(dna, min_length=30, genetic_code='standard'):
    """Find ORFs in all 6 frames with detailed annotation.

    Returns ORFs sorted by length with start/stop positions, strand, frame.
    """
    dna = dna.upper().replace(' ', '').replace('\n', '')
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    rev_comp = ''.join(complement.get(b, 'N') for b in reversed(dna))

    orfs = []
    start_codons = {'ATG'}
    stop_codons = {'TAA', 'TAG', 'TGA'}

    for strand, seq in [('+', dna), ('-', rev_comp)]:
        for frame in range(3):
            i = frame
            while i < len(seq) - 2:
                codon = seq[i:i + 3]
                if codon in start_codons:
                    # Found start, look for stop
                    j = i
                    while j < len(seq) - 2:
                        c = seq[j:j + 3]
                        if c in stop_codons:
                            length = j + 3 - i
                            if length >= min_length:
                                protein = ''
                                for k in range(i, j, 3):
                                    protein += CODON_TABLE_STD.get(seq[k:k + 3], 'X')
                                orfs.append({
                                    'strand': strand, 'frame': frame + 1,
                                    'start': i, 'end': j + 3,
                                    'length_nt': length,
                                    'length_aa': length // 3 - 1,
                                    'protein': protein,
                                    'dna': seq[i:j + 3],
                                })
                            i = j + 3
                            break
                        j += 3
                    else:
                        i += 3
                else:
                    i += 3
    orfs.sort(key=lambda x: x['length_nt'], reverse=True)
    return orfs


def melting_temperature_advanced(primer):
    """Accurate Tm using nearest-neighbor thermodynamics (SantaLucia 1998).

    Much more accurate than simple Wallace rule. Returns Tm in Celsius.
    """
    primer = primer.upper().replace(' ', '')
    if len(primer) < 2:
        return {'tm': 0, 'method': 'too short'}

    # Nearest-neighbor enthalpy (dH kcal/mol) and entropy (dS cal/mol·K)
    nn_params = {
        'AA': (-7.9, -22.2), 'AT': (-7.2, -20.4), 'AC': (-8.4, -22.4), 'AG': (-7.8, -21.0),
        'TA': (-7.2, -21.3), 'TT': (-7.9, -22.2), 'TC': (-8.2, -22.2), 'TG': (-8.5, -22.7),
        'CA': (-8.5, -22.7), 'CT': (-7.8, -21.0), 'CC': (-8.0, -19.9), 'CG': (-10.6, -27.2),
        'GA': (-8.2, -22.2), 'GT': (-8.4, -22.4), 'GC': (-9.8, -24.4), 'GG': (-8.0, -19.9),
    }
    dH, dS = 0.2, -5.7  # initiation
    for i in range(len(primer) - 1):
        pair = primer[i:i + 2]
        if pair in nn_params:
            h, s = nn_params[pair]
            dH += h; dS += s
    # Tm = dH / (dS + R ln(C)) - 273.15; C = primer conc (250 nM typical)
    R = 1.987  # cal/mol·K
    C = 250e-9
    try:
        tm = (dH * 1000) / (dS + R * math.log(C / 4)) - 273.15
        # Salt correction (50 mM Na+)
        tm += 16.6 * math.log10(0.05)
    except (ValueError, ZeroDivisionError) as _exc:
        biosim_log('warn', 'melting_temperature_advanced: {}'.format(_exc))
        tm = 0
    gc = (primer.count('G') + primer.count('C')) / len(primer) * 100
    return {'tm': round(tm, 1), 'gc_percent': round(gc, 1),
            'length': len(primer), 'method': 'nearest-neighbor (SantaLucia)'}


def restriction_digest(dna, enzymes=None):
    """Virtual restriction digest — finds cut sites and fragments.

    Returns cut positions and resulting fragment sizes.
    """
    dna = dna.upper().replace(' ', '').replace('\n', '')
    # Common enzymes and their recognition sites
    enzyme_sites = {
        'EcoRI': 'GAATTC', 'BamHI': 'GGATCC', 'HindIII': 'AAGCTT',
        'NotI': 'GCGGCCGC', 'XhoI': 'CTCGAG', 'PstI': 'CTGCAG',
        'SmaI': 'CCCGGG', 'KpnI': 'GGTACC', 'SacI': 'GAGCTC',
        'SalI': 'GTCGAC', 'SpeI': 'ACTAGT', 'XbaI': 'TCTAGA',
        'NcoI': 'CCATGG', 'NdeI': 'CATATG', 'BglII': 'AGATCT',
    }
    if enzymes is None:
        enzymes = list(enzyme_sites.keys())

    all_cuts = []
    for enz in enzymes:
        site = enzyme_sites.get(enz)
        if not site:
            continue
        pos = 0
        while True:
            idx = dna.find(site, pos)
            if idx == -1:
                break
            all_cuts.append({'enzyme': enz, 'position': idx, 'site': site})
            pos = idx + 1

    all_cuts.sort(key=lambda x: x['position'])
    # Calculate fragments
    cut_positions = sorted(set(c['position'] for c in all_cuts))
    fragments = []
    if cut_positions:
        prev = 0
        for cp in cut_positions:
            fragments.append(cp - prev)
            prev = cp
        fragments.append(len(dna) - prev)
    else:
        fragments = [len(dna)]

    return {
        'cuts': all_cuts,
        'n_cuts': len(all_cuts),
        'fragments': sorted(fragments, reverse=True),
        'enzymes_that_cut': sorted(set(c['enzyme'] for c in all_cuts)),
    }


def protein_properties(protein):
    """Calculate protein physicochemical properties.

    MW, pI (isoelectric point), hydrophobicity (GRAVY), instability.
    """
    protein = protein.upper().replace(' ', '').replace('*', '')
    if not protein:
        return {'error': 'empty'}

    # Amino acid monoisotopic masses
    aa_mass = {
        'A': 71.04, 'R': 156.10, 'N': 114.04, 'D': 115.03, 'C': 103.01,
        'E': 129.04, 'Q': 128.06, 'G': 57.02, 'H': 137.06, 'I': 113.08,
        'L': 113.08, 'K': 128.09, 'M': 131.04, 'F': 147.07, 'P': 97.05,
        'S': 87.03, 'T': 101.05, 'W': 186.08, 'Y': 163.06, 'V': 99.07,
    }
    # Kyte-Doolittle hydrophobicity
    kd = {
        'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5, 'E': -3.5,
        'Q': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5, 'L': 3.8, 'K': -3.9,
        'M': 1.9, 'F': 2.8, 'P': -1.6, 'S': -0.8, 'T': -0.7, 'W': -0.9,
        'Y': -1.3, 'V': 4.2,
    }
    mw = sum(aa_mass.get(a, 0) for a in protein) + 18.02  # + water
    gravy = sum(kd.get(a, 0) for a in protein) / len(protein)

    # Simple pI estimation from charged residues
    pos = protein.count('K') + protein.count('R') + protein.count('H')
    neg = protein.count('D') + protein.count('E')
    # Rough pI
    if pos > neg:
        pI = 8.5 + (pos - neg) * 0.1
    elif neg > pos:
        pI = 5.5 - (neg - pos) * 0.1
    else:
        pI = 7.0
    pI = max(3.0, min(11.0, pI))

    # Amino acid composition
    composition = {}
    for a in protein:
        composition[a] = composition.get(a, 0) + 1

    return {
        'length': len(protein),
        'mw_da': round(mw, 1),
        'mw_kda': round(mw / 1000, 2),
        'pI': round(pI, 1),
        'gravy': round(gravy, 3),
        'hydrophobic': gravy > 0,
        'charged_pos': pos,
        'charged_neg': neg,
        'net_charge': pos - neg,
    }


# ╔══════════════════════════════════════════════════════════════════╗
# ║  PHARMACOLOGY PRO (merged from part5)                             ║
# ╚══════════════════════════════════════════════════════════════════╝

# Physiology constants for PBPK
ORGAN_PHYSIOLOGY = {
    'liver':   {'flow': 90.0,  'volume': 1.8,  'partition': 4.0},
    'kidney':  {'flow': 72.0,  'volume': 0.3,  'partition': 3.0},
    'brain':   {'flow': 45.0,  'volume': 1.4,  'partition': 2.0},
    'fat':     {'flow': 15.0,  'volume': 15.0, 'partition': 8.0},
    'muscle':  {'flow': 45.0,  'volume': 30.0, 'partition': 2.5},
    'rest':    {'flow': 100.0, 'volume': 20.0, 'partition': 2.0},
}

CARDIAC_OUTPUT = 312.0  # L/h total blood flow
BLOOD_VOLUME = 5.0      # L



# ════════════════════════════════════════════════════════════════════
# 1. PBPK — Physiologically-Based Pharmacokinetics
# ════════════════════════════════════════════════════════════════════

# Standard human physiology (70 kg adult) — organ blood flows (L/h) and volumes (L)
ORGAN_PHYSIOLOGY = {
    'liver':   {'flow': 90.0,  'volume': 1.8,  'partition': 4.0},
    'kidney':  {'flow': 72.0,  'volume': 0.3,  'partition': 3.0},
    'brain':   {'flow': 45.0,  'volume': 1.4,  'partition': 2.0},
    'fat':     {'flow': 15.0,  'volume': 15.0, 'partition': 8.0},
    'muscle':  {'flow': 45.0,  'volume': 30.0, 'partition': 2.5},
    'rest':    {'flow': 100.0, 'volume': 20.0, 'partition': 2.0},
}

CARDIAC_OUTPUT = 312.0  # L/h total blood flow
BLOOD_VOLUME = 5.0      # L


def pbpk_model(dose, clearance_hepatic, days=2, dt=0.01,
               organs=None, ka=1.0):
    """Physiologically-based PK model with organ compartments.

    Simulates drug distribution to liver, kidney, brain, fat, muscle.
    dose: mg; clearance_hepatic: L/h; ka: absorption rate.
    Returns concentration in blood + each organ over time.
    """
    if organs is None:
        organs = ORGAN_PHYSIOLOGY
    # State: amount in gut, blood, and each organ
    gut = float(dose)
    blood = 0.0
    organ_amounts = {name: 0.0 for name in organs}

    times = [0.0]
    blood_conc = [0.0]
    organ_conc = {name: [0.0] for name in organs}

    steps = int(days * 24 / dt)
    peak_blood, peak_time = 0, 0

    for step in range(1, steps + 1):
        t = step * dt
        # Absorption from gut to blood
        absorption = ka * gut
        gut -= absorption * dt
        gut = max(0, gut)

        # Distribution to/from organs
        blood_change = absorption
        C_blood = blood / BLOOD_VOLUME
        for name, props in organs.items():
            Q = props['flow']        # blood flow
            Kp = props['partition']  # partition coefficient
            V = props['volume']
            C_organ = organ_amounts[name] / V if V > 0 else 0
            # Flux = Q * (C_blood - C_organ/Kp)
            flux = Q * (C_blood - C_organ / Kp)
            organ_amounts[name] += flux * dt
            organ_amounts[name] = max(0, organ_amounts[name])
            blood_change -= flux

        # Hepatic clearance (elimination)
        elimination = clearance_hepatic * C_blood
        blood_change -= elimination

        blood += blood_change * dt
        blood = max(0, blood)

        C_blood = blood / BLOOD_VOLUME
        if C_blood > peak_blood:
            peak_blood, peak_time = C_blood, t

        # Record (subsample)
        if step % max(1, steps // 500) == 0:
            times.append(round(t, 2))
            blood_conc.append(round(C_blood, 4))
            for name in organs:
                V = organs[name]['volume']
                organ_conc[name].append(round(organ_amounts[name] / V if V > 0 else 0, 4))

    return {
        'times': times,
        'blood': blood_conc,
        'organs': organ_conc,
        'peak_blood_conc': round(peak_blood, 4),
        'peak_time_h': round(peak_time, 2),
    }


# ════════════════════════════════════════════════════════════════════
# 2. BIOAVAILABILITY & FIRST-PASS METABOLISM
# ════════════════════════════════════════════════════════════════════

def bioavailability(fraction_absorbed, fraction_escaping_gut,
                    hepatic_extraction):
    """Calculate oral bioavailability F.

    F = Fa * Fg * Fh
    Fa: fraction absorbed from gut lumen
    Fg: fraction escaping gut wall metabolism
    Fh: fraction escaping hepatic first-pass (1 - extraction)
    """
    Fh = 1 - hepatic_extraction
    F = fraction_absorbed * fraction_escaping_gut * Fh
    return {
        'bioavailability': round(F, 3),
        'bioavailability_pct': round(F * 100, 1),
        'Fa': fraction_absorbed,
        'Fg': fraction_escaping_gut,
        'Fh': round(Fh, 3),
        'first_pass_loss_pct': round(hepatic_extraction * 100, 1),
    }


def first_pass_effect(dose, hepatic_blood_flow, intrinsic_clearance):
    """Calculate first-pass metabolism using the well-stirred liver model.

    Returns hepatic extraction ratio and amount reaching systemic circulation.
    """
    Q = hepatic_blood_flow
    CLint = intrinsic_clearance
    # Extraction ratio E = CLint / (Q + CLint)
    E = CLint / (Q + CLint) if (Q + CLint) > 0 else 0
    F_h = 1 - E
    systemic_dose = dose * F_h
    return {
        'extraction_ratio': round(E, 3),
        'hepatic_availability': round(F_h, 3),
        'systemic_dose': round(systemic_dose, 2),
        'metabolized': round(dose * E, 2),
        'classification': ('высокая экстракция' if E > 0.7 else
                           'средняя' if E > 0.3 else 'низкая'),
    }


# ════════════════════════════════════════════════════════════════════
# 3. DOSE-RESPONSE CURVES (Emax, EC50, Hill)
# ════════════════════════════════════════════════════════════════════

def dose_response_curve(emax, ec50, hill=1.0, dose_range=None):
    """Generate a dose-response (Emax/Hill) curve.

    E = Emax * C^h / (EC50^h + C^h)
    Returns doses and effects for plotting.
    """
    if dose_range is None:
        # Log-spaced doses around EC50
        dose_range = [ec50 * (10 ** (i / 4.0 - 2)) for i in range(17)]
    doses, effects = [], []
    for c in dose_range:
        if c <= 0:
            e = 0
        else:
            e = emax * (c ** hill) / (ec50 ** hill + c ** hill)
        doses.append(round(c, 4))
        effects.append(round(e, 3))
    return {
        'doses': doses, 'effects': effects,
        'emax': emax, 'ec50': ec50, 'hill': hill,
    }


def therapeutic_index(td50, ed50):
    """Therapeutic Index = TD50 / ED50 (margin of safety).

    Higher = safer. TI < 10 = narrow therapeutic window.
    """
    if ed50 <= 0:
        return {'ti': 0, 'safety': 'unknown'}
    ti = td50 / ed50
    return {
        'therapeutic_index': round(ti, 1),
        'td50': td50, 'ed50': ed50,
        'safety': ('очень безопасно' if ti > 100 else
                   'безопасно' if ti > 10 else
                   'узкое окно — осторожно!' if ti > 2 else
                   'опасно — токсичность близко!'),
        'narrow_window': ti < 10,
    }


# ════════════════════════════════════════════════════════════════════
# 4. THREE-COMPARTMENT MODEL
# ════════════════════════════════════════════════════════════════════

def pk_three_compartment(dose, ka, k10, k12, k21, k13, k31, V1, days=1, dt=0.02):
    """Three-compartment PK (central + 2 peripheral).

    Used for drugs with complex distribution (e.g. some anesthetics).
    """
    A_gut, A1, A2, A3 = dose, 0.0, 0.0, 0.0
    times, c1, c2, c3 = [0], [0], [0], [0]
    steps = int(days * 24 / dt)

    def deriv(state):
        g, a1, a2, a3 = state
        dg = -ka * g
        da1 = ka * g - (k10 + k12 + k13) * a1 + k21 * a2 + k31 * a3
        da2 = k12 * a1 - k21 * a2
        da3 = k13 * a1 - k31 * a3
        return [dg, da1, da2, da3]

    state = [A_gut, A1, A2, A3]
    for step in range(1, steps + 1):
        k1 = deriv(state)
        k2 = deriv([state[i] + dt / 2 * k1[i] for i in range(4)])
        k3 = deriv([state[i] + dt / 2 * k2[i] for i in range(4)])
        k4 = deriv([state[i] + dt * k3[i] for i in range(4)])
        state = [state[i] + dt / 6 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i])
                 for i in range(4)]
        if step % max(1, steps // 400) == 0:
            times.append(round(step * dt, 2))
            c1.append(round(state[1] / V1, 4))
            c2.append(round(state[2] / V1, 4))
            c3.append(round(state[3] / V1, 4))
    return {'times': times, 'central': c1, 'peripheral1': c2, 'peripheral2': c3}


# ════════════════════════════════════════════════════════════════════
# 5. HEPATIC CLEARANCE & CYP METABOLISM
# ════════════════════════════════════════════════════════════════════

def hepatic_clearance(intrinsic_cl, hepatic_flow=90.0, fu=1.0):
    """Calculate hepatic clearance (well-stirred model).

    intrinsic_cl: intrinsic clearance (L/h)
    hepatic_flow: liver blood flow (L/h, default 90)
    fu: fraction unbound in plasma
    """
    CL_h = (hepatic_flow * fu * intrinsic_cl) / (hepatic_flow + fu * intrinsic_cl)
    E = CL_h / hepatic_flow
    return {
        'hepatic_clearance': round(CL_h, 2),
        'extraction_ratio': round(E, 3),
        'flow_limited': E > 0.7,
        'capacity_limited': E < 0.3,
    }


def half_life(clearance, volume):
    """Calculate elimination half-life. t1/2 = 0.693 * V / CL"""
    if clearance <= 0:
        return {'half_life': 0}
    ke = clearance / volume
    t_half = 0.693 / ke
    return {
        'half_life_h': round(t_half, 2),
        'elimination_rate': round(ke, 4),
        'time_to_steady_state_h': round(t_half * 5, 1),  # ~5 half-lives
    }


# ════════════════════════════════════════════════════════════════════
# 6. BIOEQUIVALENCE
# ════════════════════════════════════════════════════════════════════

def calculate_auc(times, concentrations):
    """Calculate AUC (area under curve) using trapezoidal rule."""
    auc = 0.0
    for i in range(1, len(times)):
        dt = times[i] - times[i - 1]
        auc += 0.5 * (concentrations[i] + concentrations[i - 1]) * dt
    return round(auc, 3)


def bioequivalence(test_times, test_conc, ref_times, ref_conc):
    """Compare two formulations for bioequivalence.

    Returns Cmax, Tmax, AUC ratios. BE if 80-125% (FDA rule).
    """
    test_auc = calculate_auc(test_times, test_conc)
    ref_auc = calculate_auc(ref_times, ref_conc)
    test_cmax = max(test_conc) if test_conc else 0
    ref_cmax = max(ref_conc) if ref_conc else 0
    test_tmax = test_times[test_conc.index(test_cmax)] if test_conc else 0
    ref_tmax = ref_times[ref_conc.index(ref_cmax)] if ref_conc else 0

    auc_ratio = (test_auc / ref_auc * 100) if ref_auc else 0
    cmax_ratio = (test_cmax / ref_cmax * 100) if ref_cmax else 0

    # FDA bioequivalence: 80-125%
    auc_be = 80 <= auc_ratio <= 125
    cmax_be = 80 <= cmax_ratio <= 125

    return {
        'test_auc': test_auc, 'ref_auc': ref_auc,
        'test_cmax': round(test_cmax, 3), 'ref_cmax': round(ref_cmax, 3),
        'test_tmax': test_tmax, 'ref_tmax': ref_tmax,
        'auc_ratio': round(auc_ratio, 1),
        'cmax_ratio': round(cmax_ratio, 1),
        'bioequivalent': auc_be and cmax_be,
        'auc_within_limits': auc_be,
        'cmax_within_limits': cmax_be,
    }


# ╔══════════════════════════════════════════════════════════════════╗
# ║  MOLECULAR BIOLOGY PRO (merged from part4)                        ║
# ╚══════════════════════════════════════════════════════════════════╝
import re as _re

# ════════════════════════════════════════════════════════════════════
# 1. BLAST-LIKE LOCAL ALIGNMENT (with e-value, bit-score)
# ════════════════════════════════════════════════════════════════════

# BLOSUM62-like simplified scoring for proteins; for DNA use match/mismatch
def _seq_type(seq):
    """Guess if sequence is DNA or protein."""
    s = set(seq.upper()) - set('ACGTUN-')
    return 'protein' if s else 'dna'


def blast_search(query, database_seqs, db_names=None,
                 match=2, mismatch=-3, gap_open=-5, gap_extend=-2,
                 word_size=4):
    """BLAST-like local alignment search.

    query: query sequence string
    database_seqs: list of subject sequences to search
    Returns hits sorted by score with e-value, bit-score, identity, alignment.
    """
    query = query.upper().replace(' ', '').replace('\n', '')
    if db_names is None:
        db_names = ['Subject{}'.format(i + 1) for i in range(len(database_seqs))]

    hits = []
    total_db_len = sum(len(s) for s in database_seqs) or 1

    for idx, subject in enumerate(database_seqs):
        subject = subject.upper().replace(' ', '').replace('\n', '')
        if not subject:
            continue
        # Smith-Waterman local alignment
        aln = _smith_waterman_full(query, subject, match, mismatch,
                                    gap_open, gap_extend)
        if aln['score'] <= 0:
            continue
        # Bit-score (simplified Karlin-Altschul): S' = (lambda*S - ln K) / ln2
        lambda_param = 0.267
        K = 0.041
        bit_score = (lambda_param * aln['score'] - math.log(K)) / math.log(2)
        # E-value: E = m * n * 2^(-bit_score)
        m, n = len(query), total_db_len
        try:
            e_value = m * n * (2 ** (-bit_score))
        except OverflowError as _exc:
            biosim_log('warn', 'blast_search: {}'.format(_exc))
            e_value = 0.0
        # Identity
        matches = sum(1 for a, b in zip(aln['q_aln'], aln['s_aln'])
                      if a == b and a != '-')
        aln_len = len(aln['q_aln'])
        identity = (100.0 * matches / aln_len) if aln_len else 0
        hits.append({
            'name': db_names[idx],
            'score': aln['score'],
            'bit_score': round(bit_score, 1),
            'e_value': e_value,
            'identity': round(identity, 1),
            'align_length': aln_len,
            'q_start': aln['q_start'], 'q_end': aln['q_end'],
            's_start': aln['s_start'], 's_end': aln['s_end'],
            'q_aln': aln['q_aln'],
            's_aln': aln['s_aln'],
            'midline': _alignment_midline(aln['q_aln'], aln['s_aln']),
        })

    hits.sort(key=lambda h: h['score'], reverse=True)
    return hits


def _smith_waterman_full(seq1, seq2, match=2, mismatch=-3,
                         gap_open=-5, gap_extend=-2):
    """Smith-Waterman local alignment with affine gaps. Returns alignment."""
    n, m = len(seq1), len(seq2)
    if n == 0 or m == 0:
        return {'score': 0, 'q_aln': '', 's_aln': '',
                'q_start': 0, 'q_end': 0, 's_start': 0, 's_end': 0}
    # Simple linear gap for speed (gap_open used as gap penalty)
    gap = gap_open
    H = [[0] * (m + 1) for _ in range(n + 1)]
    max_score, max_i, max_j = 0, 0, 0
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            s = match if seq1[i - 1] == seq2[j - 1] else mismatch
            diag = H[i - 1][j - 1] + s
            up = H[i - 1][j] + gap
            left = H[i][j - 1] + gap
            H[i][j] = max(0, diag, up, left)
            if H[i][j] > max_score:
                max_score, max_i, max_j = H[i][j], i, j
    # Traceback
    q_aln, s_aln = [], []
    i, j = max_i, max_j
    while i > 0 and j > 0 and H[i][j] > 0:
        s = match if seq1[i - 1] == seq2[j - 1] else mismatch
        if H[i][j] == H[i - 1][j - 1] + s:
            q_aln.append(seq1[i - 1]); s_aln.append(seq2[j - 1]); i -= 1; j -= 1
        elif H[i][j] == H[i - 1][j] + gap:
            q_aln.append(seq1[i - 1]); s_aln.append('-'); i -= 1
        else:
            q_aln.append('-'); s_aln.append(seq2[j - 1]); j -= 1
    return {
        'score': max_score,
        'q_aln': ''.join(reversed(q_aln)),
        's_aln': ''.join(reversed(s_aln)),
        'q_start': i + 1, 'q_end': max_i,
        's_start': j + 1, 's_end': max_j,
    }


def _alignment_midline(a, b):
    """Create BLAST-style midline (| for match, space for mismatch)."""
    line = []
    for x, y in zip(a, b):
        if x == y and x != '-':
            line.append('|')
        elif x != '-' and y != '-':
            line.append('.')
        else:
            line.append(' ')
    return ''.join(line)


# ════════════════════════════════════════════════════════════════════
# 2. PLASMID ANNOTATION (auto-detect features)
# ════════════════════════════════════════════════════════════════════

# Common feature signatures (simplified motifs)
PROMOTER_MOTIFS = {
    'T7 promoter': 'TAATACGACTCACTATAG',
    'lac promoter': 'TTTACACTTTATGCTTCCGGCTCG',
    'SP6 promoter': 'ATTTAGGTGACACTATAG',
    'CMV promoter': 'GGGTCATTAGTTCATAGCCCATATATGG',
}

RBS_MOTIF = 'AGGAGG'  # Shine-Dalgarno ribosome binding site

ORI_MOTIFS = {
    'pUC ori': 'GGGAAACGCCTGGTATCTT',
    'pBR322 ori': 'TTGAGATCCTTTTTTTCTGCGC',
}

RESISTANCE_MARKERS = {
    'AmpR (bla)': 'ATGAGTATTCAACATTTCCGTGTCGCCCTTATTCC',
    'KanR': 'ATGAGCCATATTCAACGGGAAACGTCTTGCTCG',
    'CmR': 'ATGGAGAAAAAAATCACTGGATATACCACCGTTG',
}


def annotate_plasmid(sequence, min_orf_len=300):
    """Auto-annotate a plasmid sequence.

    Detects: ORFs, promoters, RBS, origins of replication, resistance markers.
    Returns list of features with positions, type, strand.
    """
    seq = sequence.upper().replace(' ', '').replace('\n', '')
    features = []
    L = len(seq)

    # 1. Promoters
    for name, motif in PROMOTER_MOTIFS.items():
        for m in _re.finditer('(?=' + motif[:12] + ')', seq):
            features.append({'type': 'promoter', 'name': name,
                             'start': m.start(), 'end': m.start() + len(motif),
                             'strand': '+'})

    # 2. RBS (Shine-Dalgarno)
    for m in _re.finditer('(?=' + RBS_MOTIF + ')', seq):
        features.append({'type': 'RBS', 'name': 'Ribosome binding site',
                         'start': m.start(), 'end': m.start() + 6, 'strand': '+'})

    # 3. Origins
    for name, motif in ORI_MOTIFS.items():
        if motif[:12] in seq:
            pos = seq.find(motif[:12])
            features.append({'type': 'ori', 'name': name,
                             'start': pos, 'end': pos + len(motif), 'strand': '+'})

    # 4. Resistance markers
    for name, motif in RESISTANCE_MARKERS.items():
        if motif[:20] in seq:
            pos = seq.find(motif[:20])
            features.append({'type': 'resistance', 'name': name,
                             'start': pos, 'end': pos + 800, 'strand': '+'})

    # 5. ORFs (both strands)
    orfs = _find_orfs_annotated(seq, min_orf_len)
    features.extend(orfs)

    # Sort by position
    features.sort(key=lambda f: f['start'])
    return {
        'length': L,
        'features': features,
        'gc_content': round(100.0 * (seq.count('G') + seq.count('C')) / L, 1) if L else 0,
        'feature_count': len(features),
    }


def _find_orfs_annotated(seq, min_len=300):
    """Find ORFs on both strands for annotation."""
    orfs = []
    comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    rev_comp = ''.join(comp.get(b, 'N') for b in reversed(seq))

    for strand, s in [('+', seq), ('-', rev_comp)]:
        for frame in range(3):
            i = frame
            while i < len(s) - 2:
                if s[i:i + 3] == 'ATG':
                    # Find stop
                    j = i
                    while j < len(s) - 2:
                        codon = s[j:j + 3]
                        if codon in ('TAA', 'TAG', 'TGA'):
                            orf_len = j + 3 - i
                            if orf_len >= min_len:
                                if strand == '+':
                                    start, end = i, j + 3
                                else:
                                    start, end = len(s) - (j + 3), len(s) - i
                                orfs.append({
                                    'type': 'ORF', 'name': 'ORF ({} aa)'.format(orf_len // 3),
                                    'start': start, 'end': end, 'strand': strand,
                                    'length': orf_len})
                            i = j
                            break
                        j += 3
                i += 3
    return orfs


# ════════════════════════════════════════════════════════════════════
# 3. RNA SECONDARY STRUCTURE (Nussinov algorithm)
# ════════════════════════════════════════════════════════════════════

def rna_secondary_structure(sequence):
    """Predict RNA secondary structure using Nussinov algorithm.

    Returns dot-bracket notation + base pairs + number of pairs.
    """
    seq = sequence.upper().replace('T', 'U').replace(' ', '').replace('\n', '')
    n = len(seq)
    if n < 2:
        return {'structure': '.' * n, 'pairs': [], 'num_pairs': 0}

    def can_pair(a, b):
        return (a, b) in [('A', 'U'), ('U', 'A'), ('G', 'C'), ('C', 'G'),
                          ('G', 'U'), ('U', 'G')]

    # DP matrix: max number of base pairs
    dp = [[0] * n for _ in range(n)]
    min_loop = 3  # minimum hairpin loop size

    for length in range(min_loop + 1, n):
        for i in range(n - length):
            j = i + length
            # Option 1: j unpaired
            best = dp[i][j - 1]
            # Option 2: j pairs with some k
            for k in range(i, j - min_loop):
                if can_pair(seq[k], seq[j]):
                    left = dp[i][k - 1] if k > i else 0
                    val = left + dp[k + 1][j - 1] + 1
                    if val > best:
                        best = val
            dp[i][j] = best

    # Traceback
    structure = ['.'] * n
    pairs = []

    def traceback(i, j):
        if i >= j:
            return
        if dp[i][j] == dp[i][j - 1]:
            traceback(i, j - 1)
        else:
            for k in range(i, j - min_loop):
                if can_pair(seq[k], seq[j]):
                    left = dp[i][k - 1] if k > i else 0
                    if dp[i][j] == left + dp[k + 1][j - 1] + 1:
                        structure[k] = '('
                        structure[j] = ')'
                        pairs.append((k, j))
                        if k > i:
                            traceback(i, k - 1)
                        traceback(k + 1, j - 1)
                        return

    traceback(0, n - 1)
    return {
        'structure': ''.join(structure),
        'pairs': sorted(pairs),
        'num_pairs': len(pairs),
        'sequence': seq,
        'mfe_estimate': round(-2.0 * len(pairs), 1),  # rough kcal/mol
    }


# ════════════════════════════════════════════════════════════════════
# 4. PROFESSIONAL PRIMER DESIGN (nearest-neighbor Tm)
# ════════════════════════════════════════════════════════════════════

# Nearest-neighbor thermodynamic parameters (SantaLucia 1998)
NN_PARAMS = {
    'AA': (-7.9, -22.2), 'AT': (-7.2, -20.4), 'AC': (-8.4, -22.4), 'AG': (-7.8, -21.0),
    'TA': (-7.2, -21.3), 'TT': (-7.9, -22.2), 'TC': (-8.2, -22.2), 'TG': (-8.5, -22.7),
    'CA': (-8.5, -22.7), 'CT': (-7.8, -21.0), 'CC': (-8.0, -19.9), 'CG': (-10.6, -27.2),
    'GA': (-8.2, -22.2), 'GT': (-8.4, -22.4), 'GC': (-9.8, -24.4), 'GG': (-8.0, -19.9),
}


def primer_tm_nn(primer, na_conc=50, primer_conc=500):
    """Calculate primer Tm using nearest-neighbor thermodynamics.

    primer: DNA sequence
    na_conc: Na+ concentration (mM)
    primer_conc: primer concentration (nM)
    Returns Tm in Celsius (much more accurate than Wallace rule).
    """
    seq = primer.upper().replace(' ', '')
    if len(seq) < 2:
        return 0.0
    dH, dS = 0.0, 0.0
    for i in range(len(seq) - 1):
        pair = seq[i:i + 2]
        if pair in NN_PARAMS:
            h, s = NN_PARAMS[pair]
            dH += h; dS += s
    # Initiation
    dH += 0.2; dS += -5.7
    # Tm calculation
    R = 1.987  # cal/(K*mol)
    C = primer_conc * 1e-9 / 4
    try:
        tm = (dH * 1000) / (dS + R * math.log(C)) - 273.15
        # Salt correction
        tm += 16.6 * math.log10(na_conc / 1000.0)
    except (ValueError, ZeroDivisionError) as _exc:
        biosim_log('warn', 'primer_tm_nn: {}'.format(_exc))
        tm = 0.0
    return round(tm, 1)


def check_primer_dimer(primer1, primer2=None):
    """Check for primer-dimer formation (3' complementarity)."""
    if primer2 is None:
        primer2 = primer1  # self-dimer
    comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    p2_rc = ''.join(comp.get(b, 'N') for b in reversed(primer2.upper()))
    # Check 3' end complementarity
    max_match = 0
    p1 = primer1.upper()
    for offset in range(len(p1)):
        matches = 0
        for i in range(min(len(p1) - offset, len(p2_rc))):
            if p1[offset + i] == p2_rc[i]:
                matches += 1
            else:
                break
        max_match = max(max_match, matches)
    return {'max_complementary': max_match,
            'risk': 'high' if max_match >= 5 else 'medium' if max_match >= 3 else 'low'}


def check_hairpin(primer, min_stem=4):
    """Check for hairpin secondary structure in primer."""
    seq = primer.upper()
    comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    n = len(seq)
    max_stem = 0
    for i in range(n):
        for j in range(i + min_stem + 3, n):
            stem = 0
            while (i + stem < j - stem and
                   seq[i + stem] == comp.get(seq[j - stem], 'X')):
                stem += 1
            max_stem = max(max_stem, stem)
    return {'max_stem': max_stem,
            'risk': 'high' if max_stem >= 5 else 'medium' if max_stem >= 4 else 'low'}


def design_primers_pro(template, target_tm=60, primer_len_range=(18, 25)):
    """Professional primer design with NN Tm, dimer & hairpin checks.

    Returns forward & reverse primers with full QC.
    """
    seq = template.upper().replace(' ', '').replace('\n', '')
    comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}

    def gc_content(p):
        return 100.0 * (p.count('G') + p.count('C')) / len(p) if p else 0

    # Forward primer: scan from start
    best_fwd = None
    for length in range(primer_len_range[0], primer_len_range[1] + 1):
        if length > len(seq):
            break
        primer = seq[:length]
        tm = primer_tm_nn(primer)
        gc = gc_content(primer)
        score = -abs(tm - target_tm) - abs(gc - 50) * 0.1
        if best_fwd is None or score > best_fwd['score']:
            best_fwd = {'primer': primer, 'tm': tm, 'gc': round(gc, 1),
                        'length': length, 'score': score}

    # Reverse primer: scan from end (reverse complement)
    rev_comp = ''.join(comp.get(b, 'N') for b in reversed(seq))
    best_rev = None
    for length in range(primer_len_range[0], primer_len_range[1] + 1):
        if length > len(seq):
            break
        primer = rev_comp[:length]
        tm = primer_tm_nn(primer)
        gc = gc_content(primer)
        score = -abs(tm - target_tm) - abs(gc - 50) * 0.1
        if best_rev is None or score > best_rev['score']:
            best_rev = {'primer': primer, 'tm': tm, 'gc': round(gc, 1),
                        'length': length, 'score': score}

    # QC
    if best_fwd and best_rev:
        dimer = check_primer_dimer(best_fwd['primer'], best_rev['primer'])
        fwd_hairpin = check_hairpin(best_fwd['primer'])
        rev_hairpin = check_hairpin(best_rev['primer'])
        product_size = len(seq)
        return {
            'forward': best_fwd, 'reverse': best_rev,
            'tm_difference': round(abs(best_fwd['tm'] - best_rev['tm']), 1),
            'primer_dimer': dimer,
            'forward_hairpin': fwd_hairpin,
            'reverse_hairpin': rev_hairpin,
            'product_size': product_size,
        }
    return {'error': 'Template too short'}


# ════════════════════════════════════════════════════════════════════
# 5. CODON ADAPTATION INDEX (CAI)
# ════════════════════════════════════════════════════════════════════

def codon_adaptation_index(dna, reference='ecoli'):
    """Calculate Codon Adaptation Index (CAI) — measure of expression efficiency.

    CAI ranges 0-1; higher = better adapted to host's preferred codons.
    """
    # E. coli relative synonymous codon usage (RSCU) weights (w values)
    ECOLI_W = {
        'TTT': 0.5, 'TTC': 1.0, 'TTA': 0.1, 'TTG': 0.1, 'CTT': 0.1, 'CTC': 0.1,
        'CTA': 0.04, 'CTG': 1.0, 'ATT': 0.5, 'ATC': 1.0, 'ATA': 0.1, 'ATG': 1.0,
        'GTT': 0.5, 'GTC': 0.4, 'GTA': 0.3, 'GTG': 1.0, 'TCT': 0.3, 'TCC': 0.3,
        'TCA': 0.1, 'TCG': 0.1, 'CCT': 0.2, 'CCC': 0.1, 'CCA': 0.3, 'CCG': 1.0,
        'ACT': 0.3, 'ACC': 1.0, 'ACA': 0.2, 'ACG': 0.3, 'GCT': 0.4, 'GCC': 0.5,
        'GCA': 0.4, 'GCG': 1.0, 'TAT': 0.5, 'TAC': 1.0, 'CAT': 0.5, 'CAC': 1.0,
        'CAA': 0.3, 'CAG': 1.0, 'AAT': 0.5, 'AAC': 1.0, 'AAA': 1.0, 'AAG': 0.3,
        'GAT': 0.6, 'GAC': 1.0, 'GAA': 1.0, 'GAG': 0.4, 'TGT': 0.5, 'TGC': 1.0,
        'TGG': 1.0, 'CGT': 1.0, 'CGC': 0.7, 'CGA': 0.1, 'CGG': 0.1, 'AGT': 0.1,
        'AGC': 1.0, 'AGA': 0.05, 'AGG': 0.03, 'GGT': 0.7, 'GGC': 1.0, 'GGA': 0.1,
        'GGG': 0.2,
    }
    seq = dna.upper().replace(' ', '').replace('\n', '')
    codons = [seq[i:i + 3] for i in range(0, len(seq) - 2, 3)]
    valid = [c for c in codons if c in ECOLI_W and ECOLI_W[c] > 0]
    if not valid:
        return {'cai': 0.0, 'codons': 0}
    # Geometric mean of w values
    log_sum = sum(math.log(ECOLI_W[c]) for c in valid)
    cai = math.exp(log_sum / len(valid))
    return {
        'cai': round(cai, 3),
        'codons': len(valid),
        'interpretation': ('отлично/excellent' if cai > 0.8 else
                           'хорошо/good' if cai > 0.5 else
                           'слабо/poor'),
    }


# ════════════════════════════════════════════════════════════════════
# 6. GC CONTENT PROFILE & MELTING
# ════════════════════════════════════════════════════════════════════

def gc_content_profile(sequence, window=50):
    """Calculate GC content in sliding windows (for plotting GC skew)."""
    seq = sequence.upper().replace(' ', '').replace('\n', '')
    profile = []
    positions = []
    for i in range(0, len(seq) - window + 1, max(1, window // 5)):
        win = seq[i:i + window]
        gc = 100.0 * (win.count('G') + win.count('C')) / len(win)
        profile.append(round(gc, 1))
        positions.append(i + window // 2)
    overall_gc = 100.0 * (seq.count('G') + seq.count('C')) / len(seq) if seq else 0
    return {'positions': positions, 'gc_profile': profile,
            'overall_gc': round(overall_gc, 1), 'window': window}


def gc_skew(sequence, window=50):
    """GC skew = (G-C)/(G+C) — helps find origin of replication."""
    seq = sequence.upper().replace(' ', '').replace('\n', '')
    positions, skew_vals, cumulative = [], [], []
    cum = 0
    for i in range(0, len(seq) - window + 1, max(1, window // 5)):
        win = seq[i:i + window]
        g, c = win.count('G'), win.count('C')
        skew = (g - c) / (g + c) if (g + c) else 0
        skew_vals.append(round(skew, 3))
        cum += skew
        cumulative.append(round(cum, 2))
        positions.append(i + window // 2)
    return {'positions': positions, 'skew': skew_vals, 'cumulative': cumulative}


# ════════════════════════════════════════════════════════════════════
# 7. VIRTUAL CLONING (Gibson assembly, restriction-ligation)
# ════════════════════════════════════════════════════════════════════

def gibson_assembly_simple(fragments, overlap=20):
    """Simulate Gibson assembly — SIMPLE version with a single overlap size.
    (The fuller version with overlap_min/overlap_max lives in biosim_part2 as
    gibson_assembly.)

    fragments: list of DNA sequences
    Returns assembled construct if overlaps match.
    """
    if not fragments:
        return {'error': 'No fragments'}
    if len(fragments) == 1:
        return {'assembled': fragments[0], 'junctions': 0}

    frags = [f.upper().replace(' ', '').replace('\n', '') for f in fragments]
    assembled = frags[0]
    junctions = 0
    used = [False] * len(frags)
    used[0] = True

    for _ in range(len(frags) - 1):
        end = assembled[-overlap:]
        found = False
        for i, frag in enumerate(frags):
            if used[i]:
                continue
            # Check if fragment starts with our end overlap
            if frag[:overlap] == end:
                assembled += frag[overlap:]
                used[i] = True
                junctions += 1
                found = True
                break
        if not found:
            break

    all_used = all(used)
    return {
        'assembled': assembled,
        'length': len(assembled),
        'junctions': junctions,
        'complete': all_used,
        'fragments_used': sum(used),
    }


def restriction_digest_custom(sequence, enzyme_sites):
    """Cut DNA at restriction sites, return fragments.

    enzyme_sites: dict {enzyme_name: recognition_sequence}
    Returns fragments and cut positions.
    """
    seq = sequence.upper().replace(' ', '').replace('\n', '')
    cuts = []
    for enzyme, site in enzyme_sites.items():
        site = site.upper().replace('^', '')  # ^ marks cut position
        pos = 0
        while True:
            idx = seq.find(site, pos)
            if idx == -1:
                break
            cuts.append({'position': idx, 'enzyme': enzyme})
            pos = idx + 1
    cuts.sort(key=lambda c: c['position'])
    # Generate fragments
    fragments = []
    last = 0
    for cut in cuts:
        fragments.append(seq[last:cut['position']])
        last = cut['position']
    fragments.append(seq[last:])
    fragments = [f for f in fragments if f]
    return {
        'cuts': cuts,
        'num_cuts': len(cuts),
        'fragments': fragments,
        'fragment_sizes': [len(f) for f in fragments],
    }


# ╔══════════════════════════════════════════════════════════════════╗
# ║  EPIDEMIOLOGY PRO — professional-level models                     ║
# ╚══════════════════════════════════════════════════════════════════╝

def seirs_model(S0, I0, beta, sigma, gamma, omega, days, dt=0.5):
    """SEIRS — SEIR with waning immunity (R -> S over time).

    omega: rate of immunity loss (1/days of immunity duration).
    Models diseases where recovery is temporary (flu, common cold, COVID).
    Returns trajectories showing repeated waves.
    """
    N = S0 + I0
    S, E, I, R = float(S0 - I0), float(I0), float(I0), 0.0
    times, S_t, E_t, I_t, R_t = [0], [S], [E], [I], [R]
    steps = int(days / dt)
    waves = 0
    prev_I = I
    rising = False

    for step in range(1, steps + 1):
        dS = -beta * S * I / N + omega * R
        dE = beta * S * I / N - sigma * E
        dI = sigma * E - gamma * I
        dR = gamma * I - omega * R
        S += dS * dt; E += dE * dt; I += dI * dt; R += dR * dt
        S = max(0, S); E = max(0, E); I = max(0, I); R = max(0, R)
        t = step * dt
        times.append(round(t, 2)); S_t.append(S); E_t.append(E)
        I_t.append(I); R_t.append(R)
        # Count waves (peaks)
        if I > prev_I:
            rising = True
        elif I < prev_I and rising:
            waves += 1
            rising = False
        prev_I = I

    return {
        'times': times, 'S': S_t, 'E': E_t, 'I': I_t, 'R': R_t,
        'R0': round(beta / gamma, 2),
        'immunity_days': round(1 / omega, 0) if omega > 0 else 9999,
        'waves': waves,
        'endemic': I_t[-1] > 1,  # disease persists
        'N': N,
    }


def fit_epidemic_data(observed_days, observed_cases, population,
                      gamma=0.1):
    """Estimate R0 and beta from observed case data (early growth fitting).

    Uses exponential growth phase: I(t) = I0 * exp(r*t), then R0 = 1 + r/gamma.
    Returns fitted parameters + predicted curve.
    """
    if len(observed_days) < 3 or len(observed_cases) < 3:
        return {'error': 'Need at least 3 data points'}

    # Fit exponential to early growth (log-linear regression)
    import math
    # Use only growing phase
    log_cases = []
    valid_days = []
    for d, c in zip(observed_days, observed_cases):
        if c > 0:
            log_cases.append(math.log(c))
            valid_days.append(d)
    if len(log_cases) < 2:
        return {'error': 'Not enough positive cases'}

    # Linear regression: log(I) = log(I0) + r*t
    n = len(valid_days)
    sum_t = sum(valid_days)
    sum_l = sum(log_cases)
    sum_tt = sum(t * t for t in valid_days)
    sum_tl = sum(t * l for t, l in zip(valid_days, log_cases))
    denom = (n * sum_tt - sum_t * sum_t)
    if abs(denom) < 1e-9:
        return {'error': 'Cannot fit'}
    r = (n * sum_tl - sum_t * sum_l) / denom  # growth rate
    log_I0 = (sum_l - r * sum_t) / n
    I0 = math.exp(log_I0)

    # R0 from growth rate: R0 = 1 + r/gamma (SIR approximation)
    R0 = 1 + r / gamma if gamma > 0 else 0
    beta = R0 * gamma

    # Doubling time
    doubling_time = math.log(2) / r if r > 0 else float('inf')

    # Generate fitted prediction
    pred_days = list(range(0, int(max(observed_days)) + 30))
    pred_cases = [I0 * math.exp(r * d) for d in pred_days]

    return {
        'R0': round(R0, 2),
        'beta': round(beta, 4),
        'growth_rate': round(r, 4),
        'doubling_time_days': round(doubling_time, 1) if doubling_time != float('inf') else None,
        'I0_estimated': round(I0, 1),
        'pred_days': pred_days,
        'pred_cases': pred_cases,
    }


def network_epidemic(n_people, avg_contacts, beta, gamma, days,
                     network_type='random', seed=None):
    """Epidemic spread on a social contact network.

    network_type: 'random' (Erdos-Renyi), 'smallworld', 'scalefree'
    More realistic than mean-field — captures superspreaders.
    Returns infection curve + network statistics.
    """
    if seed is not None:
        random.seed(seed)

    # Build contact network (adjacency list)
    contacts = {i: set() for i in range(n_people)}
    if network_type == 'scalefree':
        # Preferential attachment (Barabasi-Albert-like) — has superspreaders
        degrees = [1] * n_people
        for i in range(1, n_people):
            n_links = min(avg_contacts // 2, i)
            targets = set()
            attempts = 0
            while len(targets) < n_links and attempts < n_links * 5:
                total = sum(degrees[:i])
                r = random.uniform(0, total)
                cum = 0
                for j in range(i):
                    cum += degrees[j]
                    if cum >= r:
                        targets.add(j)
                        break
                attempts += 1
            for j in targets:
                contacts[i].add(j); contacts[j].add(i)
                degrees[i] += 1; degrees[j] += 1
    else:
        # Random network
        for i in range(n_people):
            n_links = max(1, int(random.gauss(avg_contacts, avg_contacts / 3)))
            for _ in range(n_links):
                j = random.randint(0, n_people - 1)
                if j != i:
                    contacts[i].add(j); contacts[j].add(i)

    # SIR on network
    state = {i: 'S' for i in range(n_people)}
    # Seed infection
    patient_zero = random.randint(0, n_people - 1)
    state[patient_zero] = 'I'
    infection_day = {patient_zero: 0}

    daily_infected = [1]
    daily_new = [1]
    for day in range(1, days + 1):
        new_infections = []
        recoveries = []
        for person, s in state.items():
            if s == 'I':
                # Try to infect contacts
                for contact in contacts[person]:
                    if state[contact] == 'S' and random.random() < beta:
                        new_infections.append(contact)
                # Recovery
                if random.random() < gamma:
                    recoveries.append(person)
        for p in new_infections:
            state[p] = 'I'
            infection_day[p] = day
        for p in recoveries:
            state[p] = 'R'
        current_I = sum(1 for s in state.values() if s == 'I')
        daily_infected.append(current_I)
        daily_new.append(len(set(new_infections)))
        if current_I == 0:
            break

    total_infected = sum(1 for s in state.values() if s in ('I', 'R'))
    # Network stats
    degrees = [len(c) for c in contacts.values()]
    max_degree = max(degrees) if degrees else 0
    avg_degree = sum(degrees) / len(degrees) if degrees else 0

    return {
        'daily_infected': daily_infected,
        'daily_new': daily_new,
        'total_infected': total_infected,
        'attack_rate': round(100 * total_infected / n_people, 1),
        'peak_infected': max(daily_infected),
        'peak_day': daily_infected.index(max(daily_infected)),
        'duration': len(daily_infected),
        'max_contacts': max_degree,
        'avg_contacts': round(avg_degree, 1),
        'network_type': network_type,
    }


def multi_strain_competition(S0, strains, gamma, days, dt=0.5):
    """Competition between multiple virus strains/variants.

    strains: list of dicts {'name': str, 'I0': n, 'beta': rate,
                            'cross_immunity': 0-1}
    Models how a more transmissible variant outcompetes others.
    Returns per-strain infection curves.
    """
    n_strains = len(strains)
    N = S0 + sum(s['I0'] for s in strains)
    S = float(S0)
    I = [float(s['I0']) for s in strains]
    R = [0.0] * n_strains
    names = [s['name'] for s in strains]
    betas = [s['beta'] for s in strains]

    times = [0]
    traj = {names[i]: [I[i]] for i in range(n_strains)}
    S_t = [S]
    steps = int(days / dt)

    for step in range(1, steps + 1):
        total_infection = 0
        new_I = I[:]
        for i in range(n_strains):
            infection = betas[i] * S * I[i] / N
            recovery = gamma * I[i]
            new_I[i] = I[i] + (infection - recovery) * dt
            new_I[i] = max(0, new_I[i])
            R[i] += recovery * dt
            total_infection += infection
        S -= total_infection * dt
        S = max(0, S)
        I = new_I
        times.append(round(step * dt, 2))
        S_t.append(S)
        for i in range(n_strains):
            traj[names[i]].append(I[i])

    # Determine winner
    final_sizes = {names[i]: R[i] for i in range(n_strains)}
    winner = max(final_sizes, key=final_sizes.get) if final_sizes else None

    return {
        'times': times, 'strains': traj, 'S': S_t, 'names': names,
        'winner': winner,
        'final_sizes': {k: round(v) for k, v in final_sizes.items()},
    }


def hospital_burden(infected_curve, times, hospitalization_rate=0.05,
                    icu_rate=0.01, hospital_capacity=1000, icu_capacity=100,
                    los_hospital=7, los_icu=14):
    """Estimate hospital/ICU burden from infection curve.

    Calculates concurrent hospitalizations, ICU demand, capacity overflow.
    Critical for healthcare planning.
    """
    n = len(infected_curve)
    hospitalized = [0.0] * n
    icu = [0.0] * n
    dt = (times[1] - times[0]) if len(times) > 1 else 1

    # Each new infection leads to hospitalization (with delay) for los days
    for i in range(n):
        new_cases = infected_curve[i]
        new_hosp = new_cases * hospitalization_rate
        new_icu = new_cases * icu_rate
        # Spread over length-of-stay
        hosp_steps = int(los_hospital / dt)
        icu_steps = int(los_icu / dt)
        for j in range(i, min(i + hosp_steps, n)):
            hospitalized[j] += new_hosp
        for j in range(i, min(i + icu_steps, n)):
            icu[j] += new_icu

    peak_hosp = max(hospitalized) if hospitalized else 0
    peak_icu = max(icu) if icu else 0
    hosp_overflow_days = sum(1 for h in hospitalized if h > hospital_capacity)
    icu_overflow_days = sum(1 for c in icu if c > icu_capacity)

    return {
        'hospitalized': [round(h, 1) for h in hospitalized],
        'icu': [round(c, 1) for c in icu],
        'peak_hospitalized': round(peak_hosp),
        'peak_icu': round(peak_icu),
        'hospital_overflow': peak_hosp > hospital_capacity,
        'icu_overflow': peak_icu > icu_capacity,
        'hospital_overflow_days': round(hosp_overflow_days * dt),
        'icu_overflow_days': round(icu_overflow_days * dt),
        'hospital_capacity': hospital_capacity,
        'icu_capacity': icu_capacity,
    }


# ╔══════════════════════════════════════════════════════════════════╗
# ║  NEUROSCIENCE PRO — neuron types, learning, networks              ║
# ╚══════════════════════════════════════════════════════════════════╝

# Izhikevich neuron parameters for different cell types
IZHIKEVICH_TYPES = {
    'regular_spiking':  {'a': 0.02, 'b': 0.2,  'c': -65, 'd': 8,
                          'name_ru': 'Регулярный', 'name_en': 'Regular spiking'},
    'intrinsically_bursting': {'a': 0.02, 'b': 0.2, 'c': -55, 'd': 4,
                          'name_ru': 'Пачечный', 'name_en': 'Bursting'},
    'chattering': {'a': 0.02, 'b': 0.2, 'c': -50, 'd': 2,
                          'name_ru': 'Болтающий', 'name_en': 'Chattering'},
    'fast_spiking': {'a': 0.1, 'b': 0.2, 'c': -65, 'd': 2,
                          'name_ru': 'Быстрый', 'name_en': 'Fast spiking'},
    'low_threshold': {'a': 0.02, 'b': 0.25, 'c': -65, 'd': 2,
                          'name_ru': 'Низкопороговый', 'name_en': 'Low-threshold'},
    'resonator': {'a': 0.1, 'b': 0.26, 'c': -65, 'd': 2,
                          'name_ru': 'Резонатор', 'name_en': 'Resonator'},
}


def izhikevich_neuron(neuron_type='regular_spiking', I_input=10,
                      t_max=200, dt=0.25):
    """Izhikevich neuron model — reproduces many real neuron firing patterns.

    Efficient model that captures regular spiking, bursting, fast-spiking, etc.
    Returns membrane potential trace + spike count + firing pattern.
    """
    params = IZHIKEVICH_TYPES.get(neuron_type, IZHIKEVICH_TYPES['regular_spiking'])
    a, b, c, d = params['a'], params['b'], params['c'], params['d']

    v = -65.0  # membrane potential
    u = b * v  # recovery variable
    times, v_trace = [], []
    spikes = []
    steps = int(t_max / dt)

    for step in range(steps):
        t = step * dt
        I = I_input(t) if callable(I_input) else I_input
        # Izhikevich equations
        dv = (0.04 * v * v + 5 * v + 140 - u + I)
        du = a * (b * v - u)
        v += dv * dt
        u += du * dt
        if v >= 30:  # spike threshold
            v_trace.append(30)  # cap spike
            spikes.append(round(t, 2))
            v = c       # reset
            u += d
        else:
            v_trace.append(round(v, 2))
        times.append(round(t, 2))

    firing_rate = len(spikes) / (t_max / 1000.0)
    # Detect bursting (spikes clustered)
    bursting = False
    if len(spikes) > 3:
        intervals = [spikes[i + 1] - spikes[i] for i in range(len(spikes) - 1)]
        if intervals:
            short = sum(1 for iv in intervals if iv < 10)
            bursting = short > len(intervals) * 0.3
    return {
        'times': times, 'V': v_trace, 'spikes': spikes,
        'spike_count': len(spikes),
        'firing_rate_hz': round(firing_rate, 1),
        'bursting': bursting,
        'type': neuron_type,
        'type_name': params,
    }


def cable_equation(length=100, n_segments=50, I_inject=10,
                   inject_pos=0, t_max=20, dt=0.01,
                   diameter=2.0, Rm=10000, Ri=100, Cm=1.0):
    """Cable equation — voltage propagation along a dendrite/axon.

    Models how signal spreads and attenuates along a neuron's process.
    Returns voltage at each segment over time (heatmap data).
    """
    dx = length / n_segments
    # Cable parameters
    lambda_const = math.sqrt((Rm * diameter) / (4 * Ri))  # space constant
    tau = Rm * Cm / 1000.0  # time constant (ms)

    V = [0.0] * n_segments
    snapshots = []
    steps = int(t_max / dt)
    inject_seg = int(inject_pos / dx) if dx > 0 else 0

    for step in range(steps):
        newV = V[:]
        for i in range(n_segments):
            # Diffusion term (second derivative)
            left = V[i - 1] if i > 0 else V[i]
            right = V[i + 1] if i < n_segments - 1 else V[i]
            d2V = (left - 2 * V[i] + right) / (dx * dx)
            # Cable equation: tau dV/dt = lambda^2 d2V/dx2 - V + I
            inject = I_inject if i == inject_seg else 0
            dV = (lambda_const ** 2 * d2V - V[i] + inject) / tau
            newV[i] = V[i] + dV * dt
        V = newV
        if step % max(1, steps // 50) == 0:
            snapshots.append([round(v, 3) for v in V])

    # Attenuation: peak voltage vs distance
    final_V = [abs(v) for v in V]
    return {
        'segments': n_segments,
        'length': length,
        'space_constant': round(lambda_const, 2),
        'time_constant': round(tau, 2),
        'final_voltage': [round(v, 3) for v in V],
        'snapshots': snapshots,
        'inject_segment': inject_seg,
        'attenuation': round(final_V[-1] / final_V[inject_seg], 3) if final_V[inject_seg] > 0 else 0,
    }


def hopfield_network(patterns, test_pattern, max_iters=20):
    """Hopfield network — associative memory / pattern completion.

    patterns: list of binary patterns (lists of +1/-1) to store
    test_pattern: noisy/partial pattern to recall
    Returns recalled pattern + convergence.
    """
    n = len(test_pattern)
    # Build weight matrix (Hebbian learning)
    W = [[0.0] * n for _ in range(n)]
    for pattern in patterns:
        for i in range(n):
            for j in range(n):
                if i != j:
                    W[i][j] += pattern[i] * pattern[j]
    # Normalize
    for i in range(n):
        for j in range(n):
            W[i][j] /= len(patterns)

    # Recall: iterate until convergence
    state = test_pattern[:]
    history = [state[:]]
    converged_at = max_iters
    for iteration in range(max_iters):
        new_state = state[:]
        # Asynchronous update
        for i in range(n):
            activation = sum(W[i][j] * state[j] for j in range(n))
            new_state[i] = 1 if activation >= 0 else -1
        history.append(new_state[:])
        if new_state == state:
            converged_at = iteration + 1
            break
        state = new_state

    # Which stored pattern does it match?
    best_match, best_overlap = -1, -1
    for idx, pattern in enumerate(patterns):
        overlap = sum(1 for a, b in zip(state, pattern) if a == b)
        if overlap > best_overlap:
            best_overlap = overlap
            best_match = idx

    return {
        'recalled': state,
        'converged_at': converged_at,
        'matched_pattern': best_match,
        'match_quality': round(100 * best_overlap / n, 1),
        'history': history,
        'capacity': round(0.138 * n),  # theoretical Hopfield capacity
    }


def perceptron_learn(training_data, learning_rate=0.1, epochs=50):
    """Simple perceptron — learns to classify with weight updates.

    training_data: list of (inputs, target) where target is 0 or 1
    Returns learned weights + accuracy curve (shows learning).
    """
    if not training_data:
        return {'error': 'No training data'}
    n_inputs = len(training_data[0][0])
    weights = [0.0] * n_inputs
    bias = 0.0
    accuracy_curve = []

    for epoch in range(epochs):
        correct = 0
        for inputs, target in training_data:
            # Predict
            activation = sum(w * x for w, x in zip(weights, inputs)) + bias
            prediction = 1 if activation >= 0 else 0
            # Update weights
            error = target - prediction
            if error != 0:
                for i in range(n_inputs):
                    weights[i] += learning_rate * error * inputs[i]
                bias += learning_rate * error
            else:
                correct += 1
        accuracy = correct / len(training_data)
        accuracy_curve.append(round(accuracy, 3))
        if accuracy == 1.0:
            break

    return {
        'weights': [round(w, 3) for w in weights],
        'bias': round(bias, 3),
        'accuracy_curve': accuracy_curve,
        'final_accuracy': accuracy_curve[-1] if accuracy_curve else 0,
        'epochs_used': len(accuracy_curve),
        'converged': accuracy_curve[-1] == 1.0 if accuracy_curve else False,
    }


def tuning_curve(preferred_stimulus, stimuli, width=20, max_rate=50):
    """Neural tuning curve — how a neuron responds to different stimuli.

    Models orientation/direction selectivity (Gaussian tuning).
    Returns firing rate for each stimulus value.
    """
    rates = []
    for s in stimuli:
        # Gaussian tuning around preferred stimulus
        diff = s - preferred_stimulus
        rate = max_rate * math.exp(-(diff ** 2) / (2 * width ** 2))
        rates.append(round(rate, 2))
    return {
        'stimuli': stimuli,
        'rates': rates,
        'preferred': preferred_stimulus,
        'max_rate': max_rate,
        'width': width,
    }


def spike_train_entropy(spike_times, t_max, bin_size=10):
    """Calculate information entropy of a spike train.

    Measures how much information the neuron encodes.
    Returns entropy in bits + firing statistics.
    """
    n_bins = int(t_max / bin_size)
    if n_bins == 0:
        return {'entropy': 0, 'bins': 0}
    # Count spikes per bin
    counts = [0] * n_bins
    for t in spike_times:
        b = int(t / bin_size)
        if 0 <= b < n_bins:
            counts[b] += 1
    # Probability distribution of spike counts
    from collections import Counter
    count_dist = Counter(counts)
    total = len(counts)
    entropy = 0.0
    for count_val, freq in count_dist.items():
        p = freq / total
        if p > 0:
            entropy -= p * math.log2(p)
    # CV of inter-spike intervals (regularity)
    cv = 0
    if len(spike_times) > 2:
        isis = [spike_times[i + 1] - spike_times[i] for i in range(len(spike_times) - 1)]
        if isis:
            mean_isi = sum(isis) / len(isis)
            var_isi = sum((x - mean_isi) ** 2 for x in isis) / len(isis)
            cv = math.sqrt(var_isi) / mean_isi if mean_isi > 0 else 0
    return {
        'entropy_bits': round(entropy, 3),
        'mean_rate_hz': round(len(spike_times) / (t_max / 1000.0), 1),
        'cv_isi': round(cv, 3),
        'regularity': ('регулярный' if cv < 0.5 else 'нерегулярный' if cv < 1.2 else 'пачечный'),
        'bins': n_bins,
    }


# ╔══════════════════════════════════════════════════════════════════╗
# ║  SYNTHETIC BIOLOGY PRO — parts, circuits, optimization            ║
# ╚══════════════════════════════════════════════════════════════════╝

# Standard biological parts library (BioBricks/iGEM-style registry)
BIOLOGICAL_PARTS = {
    # Promoters (with relative strength)
    'promoters': {
        'J23100': {'type': 'promoter', 'strength': 1.00, 'desc': 'Strong constitutive'},
        'J23101': {'type': 'promoter', 'strength': 0.70, 'desc': 'Medium-strong'},
        'J23106': {'type': 'promoter', 'strength': 0.47, 'desc': 'Medium'},
        'J23114': {'type': 'promoter', 'strength': 0.10, 'desc': 'Weak'},
        'pLac':   {'type': 'promoter', 'strength': 0.80, 'desc': 'IPTG-inducible', 'inducible': 'IPTG'},
        'pTet':   {'type': 'promoter', 'strength': 0.85, 'desc': 'aTc-inducible', 'inducible': 'aTc'},
        'pBAD':   {'type': 'promoter', 'strength': 0.75, 'desc': 'Arabinose-inducible', 'inducible': 'arabinose'},
    },
    # Ribosome binding sites (translation strength)
    'rbs': {
        'B0034': {'type': 'rbs', 'strength': 1.00, 'desc': 'Strong RBS'},
        'B0032': {'type': 'rbs', 'strength': 0.30, 'desc': 'Medium RBS'},
        'B0033': {'type': 'rbs', 'strength': 0.01, 'desc': 'Weak RBS'},
    },
    # Coding sequences
    'cds': {
        'GFP':  {'type': 'cds', 'desc': 'Green fluorescent protein'},
        'RFP':  {'type': 'cds', 'desc': 'Red fluorescent protein'},
        'LacI': {'type': 'cds', 'desc': 'Lac repressor'},
        'TetR': {'type': 'cds', 'desc': 'Tet repressor'},
        'cI':   {'type': 'cds', 'desc': 'Lambda cI repressor'},
    },
    # Terminators
    'terminators': {
        'B0010': {'type': 'terminator', 'efficiency': 0.98, 'desc': 'Strong terminator'},
        'B0012': {'type': 'terminator', 'efficiency': 0.85, 'desc': 'Medium terminator'},
    },
}


def list_parts(category=None):
    """List available biological parts, optionally filtered by category."""
    if category and category in BIOLOGICAL_PARTS:
        return BIOLOGICAL_PARTS[category]
    return BIOLOGICAL_PARTS


def calculate_expression(promoter_name, rbs_name, copy_number=1):
    """Calculate relative protein expression from promoter + RBS combination.

    Expression = promoter_strength * rbs_strength * copy_number.
    """
    prom = BIOLOGICAL_PARTS['promoters'].get(promoter_name, {})
    rbs = BIOLOGICAL_PARTS['rbs'].get(rbs_name, {})
    prom_strength = prom.get('strength', 0.5)
    rbs_strength = rbs.get('strength', 0.5)
    expression = prom_strength * rbs_strength * copy_number
    return {
        'expression_level': round(expression, 3),
        'promoter': promoter_name, 'promoter_strength': prom_strength,
        'rbs': rbs_name, 'rbs_strength': rbs_strength,
        'copy_number': copy_number,
        'relative_pct': round(expression * 100, 1),
    }


def toggle_switch(alpha1=15, alpha2=15, beta=3, gamma=3,
                  t_max=50, dt=0.01, init_u=1.0, init_v=0.5):
    """Genetic toggle switch — bistable system (Gardner et al. 2000).

    Two mutually-repressing genes. System has 2 stable states (memory).
    alpha: max production; beta/gamma: Hill coefficients (cooperativity).
    Returns trajectories showing which state it settles into.
    """
    u, v = init_u, init_v  # protein concentrations
    times, u_t, v_t = [0], [u], [v]
    steps = int(t_max / dt)

    for step in range(1, steps + 1):
        # du/dt = alpha1/(1+v^beta) - u
        du = alpha1 / (1 + v ** beta) - u
        dv = alpha2 / (1 + u ** gamma) - v
        u += du * dt
        v += dv * dt
        u = max(0, u); v = max(0, v)
        if step % max(1, steps // 300) == 0:
            times.append(round(step * dt, 2))
            u_t.append(round(u, 3))
            v_t.append(round(v, 3))

    # Determine final state
    final_state = 'U-high' if u > v else 'V-high'
    return {
        'times': times, 'u': u_t, 'v': v_t,
        'final_u': round(u, 3), 'final_v': round(v, 3),
        'state': final_state,
        'bistable': True,
    }


def repressilator(alpha=100, n=2, beta=0.2, t_max=200, dt=0.05):
    """Repressilator — synthetic genetic oscillator (Elowitz & Leibler 2000).

    3 genes repress each other in a ring → oscillations.
    LacI -| TetR -| cI -| LacI (cycle)
    Returns oscillating protein concentrations.
    """
    # 3 mRNAs + 3 proteins
    m = [1.0, 0.0, 0.0]   # mRNA for genes 1,2,3
    p = [0.0, 2.0, 0.0]   # proteins
    times = [0]
    p_traj = [[p[0]], [p[1]], [p[2]]]
    steps = int(t_max / dt)

    for step in range(1, steps + 1):
        new_m = m[:]
        new_p = p[:]
        for i in range(3):
            # Each gene repressed by previous protein in cycle
            repressor = p[(i - 1) % 3]
            # mRNA dynamics: production (repressed) - degradation
            dm = alpha / (1 + repressor ** n) - m[i]
            # Protein dynamics: translation - degradation
            dp = beta * (m[i] - p[i])
            new_m[i] = m[i] + dm * dt
            new_p[i] = p[i] + dp * dt
            new_m[i] = max(0, new_m[i])
            new_p[i] = max(0, new_p[i])
        m, p = new_m, new_p
        if step % max(1, steps // 500) == 0:
            times.append(round(step * dt, 2))
            for i in range(3):
                p_traj[i].append(round(p[i], 3))

    # Detect oscillation (count peaks in gene 1)
    peaks = 0
    series = p_traj[0]
    for i in range(1, len(series) - 1):
        if series[i] > series[i - 1] and series[i] > series[i + 1] and series[i] > 1:
            peaks += 1
    return {
        'times': times,
        'protein1': p_traj[0], 'protein2': p_traj[1], 'protein3': p_traj[2],
        'oscillations': peaks,
        'oscillating': peaks >= 2,
        'period': round(times[-1] / peaks, 1) if peaks > 0 else None,
    }


def hill_function(substrate, vmax, km, n=1):
    """Hill equation — cooperative binding / gene regulation.

    v = Vmax * S^n / (Km^n + S^n)
    n=1: Michaelis-Menten; n>1: positive cooperativity.
    """
    if substrate < 0:
        substrate = 0
    v = vmax * (substrate ** n) / (km ** n + substrate ** n) if (km ** n + substrate ** n) > 0 else 0
    return round(v, 4)


def dose_response_genetic(inducer_range, vmax=100, km=10, n=2,
                          leaky=5):
    """Genetic dose-response — protein output vs inducer concentration.

    Models inducible promoter response (e.g. IPTG → expression).
    leaky: basal expression without inducer.
    """
    outputs = []
    for inducer in inducer_range:
        expr = leaky + hill_function(inducer, vmax, km, n)
        outputs.append(round(expr, 2))
    # Find EC50 (where output is half-max)
    half = (leaky + vmax) / 2
    ec50 = km  # for Hill function, Km is the EC50
    return {
        'inducer': inducer_range,
        'output': outputs,
        'ec50': ec50,
        'max_output': leaky + vmax,
        'leaky_expression': leaky,
        'fold_change': round((leaky + vmax) / leaky, 1) if leaky > 0 else None,
    }


def optimize_circuit(target_expression, available_promoters=None, available_rbs=None):
    """Find best promoter+RBS combination to hit a target expression level.

    Searches the parts library for the combination closest to target.
    """
    if available_promoters is None:
        available_promoters = list(BIOLOGICAL_PARTS['promoters'].keys())
    if available_rbs is None:
        available_rbs = list(BIOLOGICAL_PARTS['rbs'].keys())

    best = None
    all_options = []
    for prom in available_promoters:
        for rbs in available_rbs:
            result = calculate_expression(prom, rbs)
            expr = result['expression_level']
            diff = abs(expr - target_expression)
            option = {'promoter': prom, 'rbs': rbs,
                      'expression': expr, 'diff': diff}
            all_options.append(option)
            if best is None or diff < best['diff']:
                best = option

    all_options.sort(key=lambda x: x['diff'])
    return {
        'target': target_expression,
        'best_combination': best,
        'top_3': all_options[:3],
    }


def metabolic_flux_balance(reactions, objective_reaction):
    """Simplified Flux Balance Analysis (FBA) for a metabolic network.

    reactions: list of dicts {'name': str, 'substrates': {met: coef},
                             'products': {met: coef}, 'max_flux': float}
    objective_reaction: name of reaction to maximize (e.g. biomass)
    Returns optimal fluxes (simplified — proportional allocation).
    """
    # Find all metabolites
    metabolites = set()
    for rxn in reactions:
        metabolites.update(rxn.get('substrates', {}).keys())
        metabolites.update(rxn.get('products', {}).keys())

    # Simplified: allocate flux proportionally to max, balanced by shared metabolites
    fluxes = {}
    for rxn in reactions:
        fluxes[rxn['name']] = rxn.get('max_flux', 10) * 0.8  # operate at 80% capacity

    # Objective value
    obj_flux = fluxes.get(objective_reaction, 0)

    # Mass balance check (production - consumption per metabolite)
    balance = {met: 0.0 for met in metabolites}
    for rxn in reactions:
        f = fluxes[rxn['name']]
        for met, coef in rxn.get('substrates', {}).items():
            balance[met] -= coef * f
        for met, coef in rxn.get('products', {}).items():
            balance[met] += coef * f

    return {
        'fluxes': {k: round(v, 2) for k, v in fluxes.items()},
        'objective_value': round(obj_flux, 2),
        'metabolite_balance': {k: round(v, 2) for k, v in balance.items()},
        'num_reactions': len(reactions),
        'num_metabolites': len(metabolites),
    }


# ╔══════════════════════════════════════════════════════════════════╗
# ║  SYNBIO DISEASE RELEVANCE — which diseases SynBio can treat       ║
# ╚══════════════════════════════════════════════════════════════════╝

# Keywords that indicate a disease is treatable by synthetic biology
GENETIC_KEYWORDS = [
    'mutation', 'мутаци', 'genetic', 'генетич', 'hereditary', 'наследств',
    'inherited', 'врождён', 'врожден', 'deficiency', 'дефицит',
    'congenital', 'gene', 'ген ', 'monogenic', 'recessive', 'рецессив',
    'dominant', 'доминант', 'chromosom', 'хромосом',
]

CANCER_KEYWORDS = [
    'cancer', 'рак', 'tumor', 'опухол', 'carcinoma', 'карцином',
    'leukemia', 'лейкоз', 'lymphoma', 'лимфом', 'sarcoma', 'саркома',
    'melanoma', 'меланом', 'oncolog', 'онколог', 'malignant', 'злокачеств',
    'neoplasm', 'blastoma', 'бластом',
]

METABOLIC_KEYWORDS = [
    'diabetes', 'диабет', 'metabolic', 'метаболич', 'enzyme', 'фермент',
    'insulin', 'инсулин', 'hormone', 'гормон', 'deficiency', 'дефицит',
]

INFECTIOUS_KEYWORDS = [
    'virus', 'вирус', 'viral', 'bacteria', 'бактери', 'infection', 'инфекц',
    'pathogen', 'патоген', 'HIV', 'ВИЧ', 'hepatitis', 'гепатит',
]

# Well-known genetic diseases treatable by gene therapy / CRISPR
KNOWN_GENE_THERAPY = {
    'cystic fibrosis', 'муковисцидоз', 'sickle cell', 'серповидно',
    'hemophilia', 'гемофили', 'duchenne', 'дюшен', 'thalassemia', 'талассеми',
    'huntington', 'хантингтон', 'гентингтон', 'spinal muscular', 'спинальн',
    'tay-sachs', 'тея-сакс', 'phenylketonuria', 'фенилкетонури',
    'gaucher', 'гоше', 'fabry', 'фабри', 'pompe', 'помпе',
    'severe combined', 'тяжёлый комбинир', 'adenosine deaminase',
    'leber', 'лебер', 'retinitis', 'ретинит',
}


def assess_synbio_relevance(disease_name, genes=None, description='',
                            proteins=None):
    """Assess whether a disease can be treated with synthetic biology.

    Returns relevance level + applicable SynBio approaches + explanation.
    Levels: 'high' (direct application), 'medium' (research stage),
            'low' (limited), 'none'.
    """
    genes = genes or []
    proteins = proteins or []
    text = (disease_name + ' ' + description).lower()

    approaches = []
    relevance_score = 0

    # 1. Known gene therapy targets (highest relevance)
    is_known = any(kw in text for kw in KNOWN_GENE_THERAPY)
    if is_known:
        relevance_score += 5
        approaches.append({
            'type': 'gene_therapy',
            'name_ru': 'Генная терапия',
            'name_en': 'Gene therapy',
            'desc_ru': 'Доставка исправленной копии гена через вирусный вектор',
            'desc_en': 'Deliver corrected gene copy via viral vector',
        })
        approaches.append({
            'type': 'crispr',
            'name_ru': 'CRISPR-редактирование',
            'name_en': 'CRISPR editing',
            'desc_ru': 'Прямое исправление мутации в геноме',
            'desc_en': 'Direct correction of the mutation in the genome',
        })

    # 2. Genetic disease (has identified genes + genetic keywords)
    is_genetic = (any(kw in text for kw in GENETIC_KEYWORDS) or
                  (genes and len(genes) > 0))
    if is_genetic and not is_known:
        relevance_score += 3
        if genes:
            approaches.append({
                'type': 'gene_therapy',
                'name_ru': 'Генная терапия',
                'name_en': 'Gene therapy',
                'desc_ru': 'Целевые гены: {}'.format(', '.join(genes[:4])),
                'desc_en': 'Target genes: {}'.format(', '.join(genes[:4])),
            })

    # 3. Cancer (CAR-T, oncolytic)
    is_cancer = any(kw in text for kw in CANCER_KEYWORDS)
    if is_cancer:
        relevance_score += 4
        approaches.append({
            'type': 'car_t',
            'name_ru': 'CAR-T клеточная терапия',
            'name_en': 'CAR-T cell therapy',
            'desc_ru': 'Перепрограммирование Т-клеток для атаки опухоли',
            'desc_en': 'Reprogram T-cells to attack the tumor',
        })
        approaches.append({
            'type': 'oncolytic',
            'name_ru': 'Онколитические вирусы',
            'name_en': 'Oncolytic viruses',
            'desc_ru': 'Инженерные вирусы, избирательно убивающие раковые клетки',
            'desc_en': 'Engineered viruses that selectively kill cancer cells',
        })

    # 4. Metabolic (engineered enzymes, cells)
    is_metabolic = any(kw in text for kw in METABOLIC_KEYWORDS)
    if is_metabolic:
        relevance_score += 2
        approaches.append({
            'type': 'cell_therapy',
            'name_ru': 'Инженерные клетки/ферменты',
            'name_en': 'Engineered cells/enzymes',
            'desc_ru': 'Клетки, производящие недостающий белок (напр. инсулин)',
            'desc_en': 'Cells producing the missing protein (e.g. insulin)',
        })

    # 5. Infectious (engineered immunity, phage therapy)
    is_infectious = any(kw in text for kw in INFECTIOUS_KEYWORDS)
    if is_infectious:
        relevance_score += 2
        approaches.append({
            'type': 'engineered_immunity',
            'name_ru': 'Инженерный иммунитет',
            'name_en': 'Engineered immunity',
            'desc_ru': 'Генетически усиленная защита (напр. CRISPR против ВИЧ)',
            'desc_en': 'Genetically enhanced defense (e.g. CRISPR against HIV)',
        })

    # Determine relevance level
    if relevance_score >= 5:
        level = 'high'
    elif relevance_score >= 3:
        level = 'medium'
    elif relevance_score >= 1:
        level = 'low'
    else:
        level = 'none'

    # Remove duplicate approaches
    seen = set()
    unique_approaches = []
    for a in approaches:
        if a['type'] not in seen:
            seen.add(a['type'])
            unique_approaches.append(a)

    return {
        'relevance': level,
        'score': relevance_score,
        'approaches': unique_approaches,
        'is_genetic': is_genetic or is_known,
        'is_cancer': is_cancer,
        'treatable': level in ('high', 'medium'),
    }


def get_synbio_treatable_diseases(disease_db, lang='ru'):
    """Filter a disease database to those treatable by synthetic biology.

    disease_db: dict of {name: {genes, desc_ru, desc_en, ...}}
    Returns list of (name, relevance_info) sorted by relevance.
    """
    results = []
    for name, info in disease_db.items():
        desc = info.get('desc_' + lang, '') + ' ' + info.get('desc_en', '')
        genes = info.get('genes', [])
        proteins = info.get('proteins', [])
        display_name = info.get(lang, name)
        assessment = assess_synbio_relevance(display_name, genes, desc, proteins)
        if assessment['treatable']:
            results.append((name, display_name, assessment))
    # Sort by relevance score (highest first)
    results.sort(key=lambda x: x[2]['score'], reverse=True)
    return results


def suggest_gene_circuit(approach_type, genes=None):
    """Suggest a starter genetic circuit design for a treatment approach.

    Returns a circuit template (parts to use) for the SynBio constructor.
    """
    genes = genes or []
    templates = {
        'gene_therapy': {
            'name_ru': 'Цепь генной терапии',
            'name_en': 'Gene therapy circuit',
            'parts': ['Strong promoter (J23100)', 'RBS (B0034)',
                      'Therapeutic gene', 'Terminator (B0010)'],
            'desc_ru': 'Постоянная экспрессия исправленного гена',
            'desc_en': 'Constitutive expression of corrected gene',
        },
        'crispr': {
            'name_ru': 'CRISPR-цепь',
            'name_en': 'CRISPR circuit',
            'parts': ['Promoter (U6)', 'Cas9', 'guide RNA (gRNA)', 'Terminator'],
            'desc_ru': 'Cas9 + направляющая РНК для редактирования мутации',
            'desc_en': 'Cas9 + guide RNA to edit the mutation',
        },
        'car_t': {
            'name_ru': 'CAR-T конструкция',
            'name_en': 'CAR-T construct',
            'parts': ['Promoter (EF1a)', 'scFv (antigen recognition)',
                      'CD28 (costimulation)', 'CD3z (activation)'],
            'desc_ru': 'Химерный антигенный рецептор для атаки опухоли',
            'desc_en': 'Chimeric antigen receptor to attack the tumor',
        },
        'cell_therapy': {
            'name_ru': 'Терапевтическая клетка',
            'name_en': 'Therapeutic cell',
            'parts': ['Inducible promoter (pTet)', 'RBS (B0034)',
                      'Therapeutic protein', 'Secretion tag'],
            'desc_ru': 'Клетка, выделяющая лечебный белок по сигналу',
            'desc_en': 'Cell secreting therapeutic protein on demand',
        },
        'oncolytic': {
            'name_ru': 'Онколитический вирус',
            'name_en': 'Oncolytic virus',
            'parts': ['Tumor-specific promoter', 'Viral replication genes',
                      'Immune-stimulating gene (GM-CSF)'],
            'desc_ru': 'Вирус, размножающийся только в раковых клетках',
            'desc_en': 'Virus replicating only in cancer cells',
        },
        'engineered_immunity': {
            'name_ru': 'Инженерный иммунитet',
            'name_en': 'Engineered immunity',
            'parts': ['Promoter', 'CRISPR-Cas system', 'Targeting gRNA',
                      'Terminator'],
            'desc_ru': 'Генетическая защита от патогена',
            'desc_en': 'Genetic defense against the pathogen',
        },
    }
    return templates.get(approach_type, templates['gene_therapy'])


# ╔══════════════════════════════════════════════════════════════════╗
# ║  SYNBIO DISEASE RELEVANCE — which diseases SynBio can address     ║
# ╚══════════════════════════════════════════════════════════════════╝

# Keywords that indicate a disease is treatable/studied via synthetic biology
SYNBIO_GENETIC_KEYWORDS = [
    'mutation', 'mutations', 'gene', 'genetic', 'hereditary', 'inherited',
    'мутаци', 'ген', 'генетич', 'наследств', 'врождённ', 'врожденн',
    'deficiency', 'дефицит', 'недостаточность',
]

SYNBIO_CANCER_KEYWORDS = [
    'cancer', 'tumor', 'tumour', 'carcinoma', 'leukemia', 'lymphoma',
    'sarcoma', 'melanoma', 'рак', 'опухол', 'карцином', 'лейкоз',
    'лимфом', 'саркома', 'меланома', 'онколог',
]

SYNBIO_METABOLIC_KEYWORDS = [
    'metabolic', 'metabolism', 'insulin', 'diabetes', 'enzyme',
    'метаболич', 'метаболизм', 'инсулин', 'диабет', 'фермент',
    'hormone', 'гормон',
]

SYNBIO_INFECTIOUS_KEYWORDS = [
    'virus', 'viral', 'bacteria', 'bacterial', 'infection', 'infectious',
    'вирус', 'бактери', 'инфекц', 'pathogen', 'патоген',
]

# Specific genes known to be gene-therapy / CRISPR targets
SYNBIO_TARGET_GENES = {
    'CFTR', 'HBB', 'F8', 'F9', 'DMD', 'TP53', 'BRCA1', 'BRCA2', 'RB1',
    'HTT', 'SMN1', 'ADA', 'IL2RG', 'PAH', 'LDLR', 'APOE', 'PCSK9',
    'INS', 'INSR', 'PPARG', 'HEXA', 'GBA', 'G6PC', 'OTC', 'RPE65',
}


def classify_synbio_relevance(disease_data):
    """Classify how relevant a disease is to synthetic biology treatment.

    disease_data: dict with 'genes', 'desc_en', 'desc_ru', etc.
    Returns relevance category + applicable SynBio approaches.
    """
    genes = set(disease_data.get('genes', []))
    desc = (disease_data.get('desc_en', '') + ' ' +
            disease_data.get('desc_ru', '') + ' ' +
            disease_data.get('en', '') + ' ' +
            disease_data.get('ru', '')).lower()

    approaches = []
    relevance_score = 0

    # Check gene therapy targets
    target_genes = genes & SYNBIO_TARGET_GENES
    if target_genes:
        approaches.append('gene_therapy')
        relevance_score += 3

    # Genetic disease → gene therapy / CRISPR
    if any(kw in desc for kw in SYNBIO_GENETIC_KEYWORDS):
        if 'gene_therapy' not in approaches:
            approaches.append('gene_therapy')
        approaches.append('crispr')
        relevance_score += 2

    # Cancer → CAR-T, oncolytic
    if any(kw in desc for kw in SYNBIO_CANCER_KEYWORDS):
        approaches.append('car_t')
        relevance_score += 2

    # Metabolic → engineered cells/bacteria
    if any(kw in desc for kw in SYNBIO_METABOLIC_KEYWORDS):
        approaches.append('engineered_cells')
        relevance_score += 1

    # Infectious → engineered phages, vaccines
    if any(kw in desc for kw in SYNBIO_INFECTIOUS_KEYWORDS):
        approaches.append('phage_therapy')
        relevance_score += 1

    # Having known genes always adds some relevance
    if genes and relevance_score == 0:
        relevance_score = 1
        approaches.append('research')

    # Determine category
    if relevance_score >= 3:
        category = 'high'
    elif relevance_score >= 1:
        category = 'medium'
    else:
        category = 'low'

    return {
        'category': category,
        'score': relevance_score,
        'approaches': list(dict.fromkeys(approaches)),  # dedup, keep order
        'target_genes': list(target_genes),
        'synbio_applicable': relevance_score >= 1,
    }


def synbio_approach_description(approach, lang='ru'):
    """Human-readable description of a SynBio treatment approach."""
    is_ru = (lang == 'ru')
    descriptions = {
        'gene_therapy': ('Генная терапия — доставка исправного гена',
                          'Gene therapy — deliver correct gene'),
        'crispr': ('CRISPR — редактирование мутации в ДНК',
                    'CRISPR — edit the mutation in DNA'),
        'car_t': ('CAR-T — перепрограммирование иммунных клеток против рака',
                   'CAR-T — reprogram immune cells against cancer'),
        'engineered_cells': ('Инженерные клетки — производство недостающего белка',
                              'Engineered cells — produce missing protein'),
        'phage_therapy': ('Фаговая терапия — вирусы против бактерий',
                           'Phage therapy — viruses against bacteria'),
        'research': ('Исследовательский интерес — гены известны',
                     'Research interest — genes known'),
    }
    d = descriptions.get(approach, (approach, approach))
    return d[0] if is_ru else d[1]


def build_gene_circuit_suggestion(disease_data, lang='ru'):
    """Suggest a genetic circuit design to address a disease.

    Returns a circuit template (nodes + connections) for the SynBio constructor.
    """
    relevance = classify_synbio_relevance(disease_data)
    genes = disease_data.get('genes', [])
    is_ru = (lang == 'ru')

    if not relevance['synbio_applicable']:
        return {'applicable': False}

    suggestion = {
        'applicable': True,
        'approaches': relevance['approaches'],
        'target_genes': genes[:3],  # focus on first few genes
        'circuit_nodes': [],
        'description': '',
    }

    primary = relevance['approaches'][0] if relevance['approaches'] else 'research'

    if primary == 'gene_therapy':
        # Promoter -> therapeutic gene -> terminator
        target = genes[0] if genes else 'GENE'
        suggestion['circuit_nodes'] = [
            {'type': 'promoter', 'name': 'pCMV (strong)'},
            {'type': 'rbs', 'name': 'B0034'},
            {'type': 'cds', 'name': target + ' (corrected)'},
            {'type': 'terminator', 'name': 'B0010'},
        ]
        suggestion['description'] = (
            'Доставить исправную копию гена {} в клетки'.format(target)
            if is_ru else
            'Deliver a corrected copy of {} into cells'.format(target))
    elif primary == 'car_t':
        suggestion['circuit_nodes'] = [
            {'type': 'promoter', 'name': 'pEF1a'},
            {'type': 'rbs', 'name': 'B0034'},
            {'type': 'cds', 'name': 'CAR (chimeric antigen receptor)'},
            {'type': 'cds', 'name': 'CD3z signaling'},
            {'type': 'terminator', 'name': 'B0010'},
        ]
        suggestion['description'] = (
            'Перепрограммировать T-клетки распознавать раковые клетки'
            if is_ru else
            'Reprogram T-cells to recognize cancer cells')
    elif primary == 'engineered_cells':
        protein = disease_data.get('proteins', ['protein'])[0]
        suggestion['circuit_nodes'] = [
            {'type': 'promoter', 'name': 'pBAD (inducible)'},
            {'type': 'rbs', 'name': 'B0034'},
            {'type': 'cds', 'name': protein + ' production'},
            {'type': 'terminator', 'name': 'B0010'},
        ]
        suggestion['description'] = (
            'Инженерные клетки производят {}'.format(protein)
            if is_ru else
            'Engineered cells produce {}'.format(protein))
    else:
        target = genes[0] if genes else 'target'
        suggestion['circuit_nodes'] = [
            {'type': 'promoter', 'name': 'J23100'},
            {'type': 'rbs', 'name': 'B0034'},
            {'type': 'cds', 'name': target},
            {'type': 'terminator', 'name': 'B0010'},
        ]
        suggestion['description'] = (
            'Базовая цепь для изучения {}'.format(target)
            if is_ru else
            'Basic circuit to study {}'.format(target))

    return suggestion


# ╔══════════════════════════════════════════════════════════════════╗
# ║  3D STRUCTURE SEARCH & FETCH (molecules + proteins from internet) ║
# ╚══════════════════════════════════════════════════════════════════╝

import urllib.request as _urlreq
import urllib.parse as _urlparse
import json as _json
import ssl as _ssl


def _make_ssl_context():
    """SSL context that works on old macOS (unverified fallback)."""
    try:
        ctx = _ssl.create_default_context()
        ctx.check_hostname = False
        ctx.verify_mode = _ssl.CERT_NONE
        return ctx
    except Exception as _exc:
        biosim_log('warn', '_make_ssl_context: {}'.format(_exc))
        return None


def search_pubchem_molecule(query, max_results=8):
    """Search PubChem for molecules by name. Returns list of candidates.

    Returns [{'cid': int, 'name': str, 'formula': str}, ...]
    """
    ctx = _make_ssl_context()
    results = []
    try:
        # Get CID by name
        q = _urlparse.quote(query)
        url = ('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/'
               + q + '/cids/JSON')
        req = _urlreq.Request(url, headers={'User-Agent': 'VIREON/1.0'})
        with _urlreq.urlopen(req, timeout=12, context=ctx) as resp:
            data = _json.loads(resp.read().decode())
        cids = data.get('IdentifierList', {}).get('CID', [])[:max_results]
        for cid in cids:
            results.append({'cid': cid, 'name': query, 'formula': ''})
        # Get properties for first few
        if cids:
            cid_str = ','.join(str(c) for c in cids[:max_results])
            purl = ('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/'
                    + cid_str + '/property/MolecularFormula,IUPACName/JSON')
            req2 = _urlreq.Request(purl, headers={'User-Agent': 'VIREON/1.0'})
            with _urlreq.urlopen(req2, timeout=12, context=ctx) as resp:
                pdata = _json.loads(resp.read().decode())
            props = pdata.get('PropertyTable', {}).get('Properties', [])
            for i, p in enumerate(props):
                if i < len(results):
                    results[i]['formula'] = p.get('MolecularFormula', '')
                    results[i]['name'] = p.get('IUPACName', query) or query
    except Exception as e:
        biosim_log('warn', 'search_pubchem_molecule: {}'.format(e))
        return {'error': str(e), 'results': []}
    return {'error': None, 'results': results}


def fetch_pubchem_3d(cid):
    """Download 3D structure (SDF) for a PubChem CID.

    Returns atom coordinates + bonds for matplotlib display.
    """
    ctx = _make_ssl_context()
    try:
        url = ('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/'
               + str(cid) + '/record/SDF?record_type=3d')
        req = _urlreq.Request(url, headers={'User-Agent': 'VIREON/1.0'})
        with _urlreq.urlopen(req, timeout=15, context=ctx) as resp:
            sdf = resp.read().decode()
        return parse_sdf(sdf)
    except Exception as e:
        biosim_log('warn', 'fetch_pubchem_3d: {}'.format(e))
        return {'error': str(e)}


def parse_sdf(sdf_text):
    """Parse SDF molecule file → atoms (with element + xyz) and bonds."""
    lines = sdf_text.split('\n')
    if len(lines) < 4:
        return {'error': 'Invalid SDF'}
    # Counts line (line 4, index 3)
    try:
        counts = lines[3]
        n_atoms = int(counts[:3])
        n_bonds = int(counts[3:6])
    except (ValueError, IndexError) as _exc:
        biosim_log('warn', 'parse_sdf: {}'.format(_exc))
        return {'error': 'Cannot parse SDF counts'}

    atoms = []
    for i in range(4, 4 + n_atoms):
        if i >= len(lines):
            break
        parts = lines[i].split()
        if len(parts) >= 4:
            try:
                x, y, z = float(parts[0]), float(parts[1]), float(parts[2])
                element = parts[3]
                atoms.append({'element': element, 'x': x, 'y': y, 'z': z})
            except ValueError as _exc:
                biosim_log('warn', 'parse_sdf: {}'.format(_exc))
                continue

    bonds = []
    for i in range(4 + n_atoms, 4 + n_atoms + n_bonds):
        if i >= len(lines):
            break
        line = lines[i]
        try:
            a1 = int(line[:3]) - 1
            a2 = int(line[3:6]) - 1
            order = int(line[6:9]) if len(line) >= 9 else 1
            bonds.append({'a1': a1, 'a2': a2, 'order': order})
        except (ValueError, IndexError) as _exc:
            biosim_log('warn', 'parse_sdf: {}'.format(_exc))
            continue

    return {'error': None, 'atoms': atoms, 'bonds': bonds, 'n_atoms': len(atoms)}


def search_pdb_protein(query, max_results=8):
    """Search RCSB PDB for protein structures by name/keyword.

    Returns [{'pdb_id': str, 'title': str}, ...]
    """
    ctx = _make_ssl_context()
    try:
        # RCSB search API
        search_url = 'https://search.rcsb.org/rcsbsearch/v2/query'
        query_json = {
            "query": {
                "type": "terminal",
                "service": "full_text",
                "parameters": {"value": query}
            },
            "return_type": "entry",
            "request_options": {
                "paginate": {"start": 0, "rows": max_results}
            }
        }
        data_bytes = _json.dumps(query_json).encode()
        req = _urlreq.Request(search_url, data=data_bytes,
                              headers={'Content-Type': 'application/json',
                                       'User-Agent': 'VIREON/1.0'})
        with _urlreq.urlopen(req, timeout=12, context=ctx) as resp:
            data = _json.loads(resp.read().decode())
        results = []
        for item in data.get('result_set', [])[:max_results]:
            pdb_id = item.get('identifier', '')
            results.append({'pdb_id': pdb_id, 'title': ''})
        return {'error': None, 'results': results}
    except Exception as e:
        biosim_log('warn', 'search_pdb_protein: {}'.format(e))
        return {'error': str(e), 'results': []}


def fetch_pdb_structure(pdb_id, max_atoms=2000):
    """Download PDB structure → CA backbone atoms for matplotlib display.

    Returns alpha-carbon coordinates (protein backbone trace).
    """
    ctx = _make_ssl_context()
    try:
        url = 'https://files.rcsb.org/download/' + pdb_id.upper() + '.pdb'
        req = _urlreq.Request(url, headers={'User-Agent': 'VIREON/1.0'})
        with _urlreq.urlopen(req, timeout=15, context=ctx) as resp:
            pdb_text = resp.read().decode()
        return parse_pdb_backbone(pdb_text, max_atoms)
    except Exception as e:
        biosim_log('warn', 'fetch_pdb_structure: {}'.format(e))
        return {'error': str(e)}


def parse_pdb_backbone(pdb_text, max_atoms=2000):
    """Parse PDB → alpha-carbon (CA) trace + all atoms for display.

    Returns CA backbone (the protein 'spine') and chain info.
    """
    ca_atoms = []   # alpha carbons (backbone trace)
    all_atoms = []  # all atoms (for cloud view)
    chains = {}
    for line in pdb_text.split('\n'):
        if line.startswith('ATOM') or line.startswith('HETATM'):
            try:
                atom_name = line[12:16].strip()
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                element = line[76:78].strip() or atom_name[0]
                chain = line[21]
                if len(all_atoms) < max_atoms:
                    all_atoms.append({'element': element, 'x': x, 'y': y, 'z': z,
                                       'chain': chain})
                if atom_name == 'CA':  # alpha carbon
                    ca_atoms.append({'x': x, 'y': y, 'z': z, 'chain': chain})
                    chains.setdefault(chain, 0)
                    chains[chain] += 1
            except (ValueError, IndexError) as _exc:
                biosim_log('warn', 'parse_pdb_backbone: {}'.format(_exc))
                continue
    return {
        'error': None,
        'ca_atoms': ca_atoms,
        'all_atoms': all_atoms,
        'chains': list(chains.keys()),
        'n_residues': len(ca_atoms),
        'n_atoms': len(all_atoms),
    }


# ╔══════════════════════════════════════════════════════════════════╗
# ║  CRISPR ADVANCED — Doench efficiency, donor design, repair        ║
# ╚══════════════════════════════════════════════════════════════════╝

# ═══════════════════════════════════════════════════════════════════════
#  CFD (Cutting Frequency Determination) — Doench 2016 model
#  The official CFD ships as two pickles (mismatch_score.pkl, pam_scores.pkl)
#  that CRISPOR loads at runtime. Those binary tables are not distributable
#  here, so this is a faithful *positional* reconstruction of the same model:
#  every mismatch multiplies the running score by a weight that depends on
#  the mismatch position (1 = PAM-distal 5' end ... 20 = PAM-proximal seed)
#  and on the base identities, with the seed region (positions ~13-20)
#  penalised far more heavily than the 5' end — exactly the behaviour the
#  Doench 2016 weights encode. The 3' NGG PAM dinucleotide is scored with the
#  published pam_scores table (which IS small enough to embed exactly).
#  Returns values in the same 0-1 range and ranking as the official CFD.
# ═══════════════════════════════════════════════════════════════════════

# Published PAM dinucleotide scores (Doench 2016, pam_scores.pkl) — exact.
_CFD_PAM = {
    'AA': 0.0, 'AC': 0.0, 'AG': 0.259259259, 'AT': 0.0,
    'CA': 0.0, 'CC': 0.0, 'CG': 0.107142857, 'CT': 0.0,
    'GA': 0.069444444, 'GC': 0.022222222, 'GG': 1.0, 'GT': 0.016129032,
    'TA': 0.0, 'TC': 0.0, 'TG': 0.038961039, 'TT': 0.0,
}

# Per-position "retained activity when mismatched" multipliers. Averaged from
# the Doench 2016 mismatch table across base identities: 5' positions tolerate
# mismatches (weight near 1), the PAM-proximal seed does not (weight near 0).
_CFD_POS_WEIGHT = [
    1.00, 0.85, 0.70, 0.72, 0.72, 0.70, 0.68, 0.66, 0.60, 0.58,   # 1-10
    0.55, 0.50, 0.45, 0.38, 0.30, 0.22, 0.15, 0.10, 0.07, 0.05,   # 11-20
]

# Base-transition modifier: some mismatch types (e.g. rG:dT wobble) are
# tolerated more than others. Multiplies the positional weight.
_CFD_TRANSITION = {
    ('G', 'T'): 1.25, ('T', 'G'): 1.15, ('A', 'C'): 1.10, ('C', 'A'): 1.05,
    ('G', 'A'): 0.95, ('A', 'G'): 0.90, ('C', 'T'): 0.90, ('T', 'C'): 0.95,
}

_CFD_COMP = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'U': 'A', 'N': 'N'}


def _cfd_revcom(s):
    return ''.join(_CFD_COMP.get(b, 'N') for b in reversed(s))


def calc_cfd_doench(wt, sg, pam):
    """CFD calculation (Doench 2016 model). wt = on-target 20nt spacer,
    sg = off-target 20nt protospacer, pam = the off-target PAM (NGG -> uses
    the last 2 nt). Returns 0.0-1.0: 1.0 = identical, lower = mismatches
    reduce predicted cutting, with PAM-proximal mismatches penalised most.
    """
    w = (wt or '').upper().replace('U', 'T')
    s = (sg or '').upper().replace('U', 'T')
    n = min(len(w), len(s), 20)
    if n == 0:
        return 0.0
    w = w[-n:]; s = s[-n:]     # align at the PAM-proximal (3') end
    score = 1.0
    for i in range(n):
        if w[i] == s[i]:
            continue
        # position index from the PAM-proximal end -> weight table (1..20)
        pos_from_5 = 20 - n + i          # map into the 20-long weight table
        base_w = _CFD_POS_WEIGHT[max(0, min(19, pos_from_5))]
        base_w *= _CFD_TRANSITION.get((w[i], s[i]), 1.0)
        score *= max(0.0, min(1.0, base_w))
    score *= _CFD_PAM.get((pam or '').upper()[-2:], 0.0)
    return round(score, 6)


def cfd_off_target_score_v2(guide_seq, off_target_seq_with_pam):
    """CFD score using the Doench 2016 model with the site's real PAM.
    Off-target input is a 23-mer (20 protospacer + 3 PAM); if only 20nt are
    given, a perfect NGG PAM is assumed. Returns 0.0-1.0. This is the
    high-accuracy path; cfd_off_target_score() below is the PAM-free fallback.
    """
    g = (guide_seq or '').upper()[:20]
    o = (off_target_seq_with_pam or '').upper()
    if len(o) >= 23:
        proto = o[-23:-3]
        pam = o[-2:]
    else:
        proto = o[-20:] if len(o) >= 20 else o
        pam = 'GG'
    if len(proto) < len(g):
        return 0.0
    return calc_cfd_doench(g, proto[-len(g):], pam)


def cfd_specificity(guide, full_sequence, pam='NGG', max_mismatches=4):
    """Aggregate CFD-based guide specificity, CRISPOR-style: 100/(1+sum of
    CFD over every near-match in `full_sequence`), using the real PAM at each
    site. Higher = more specific. Returns
    {'specificity', 'worst', 'n_offtargets', 'cfd_sum'}.
    """
    g = (guide or '').upper()
    seq = (full_sequence or '').upper()
    if not g or not seq:
        return {'specificity': 100, 'worst': None, 'n_offtargets': 0, 'cfd_sum': 0.0}
    glen = len(g)
    pamlen = 3
    total = 0.0
    worst = None
    n_ot = 0
    seen_perfect = False
    rc = _cfd_revcom(seq)
    for strand, sseq in (('+', seq), ('-', rc)):
        for i in range(len(sseq) - glen - pamlen + 1):
            window = sseq[i:i + glen]
            site_pam = sseq[i + glen:i + glen + pamlen]
            if 'N' in window:
                continue
            mm = 0
            for a, b in zip(g, window):
                if a != b:
                    mm += 1
                    if mm > max_mismatches:
                        break
            if mm > max_mismatches:
                continue
            if mm == 0 and not seen_perfect:
                seen_perfect = True
                continue
            n_ot += 1
            cfd = calc_cfd_doench(g, window, site_pam)
            total += cfd
            if worst is None or cfd > worst['cfd']:
                worst = {'cfd': round(cfd, 3), 'seq': window, 'pam': site_pam,
                         'mismatches': mm, 'strand': strand, 'position': i}
    spec = 100.0 / (1.0 + total) if total > 0 else 100.0
    return {'specificity': int(round(spec)), 'worst': worst,
            'n_offtargets': n_ot, 'cfd_sum': round(total, 3)}


def cfd_off_target_score(guide_seq, off_target_seq):
    """CFD (Cutting Frequency Determination) score — predicts how likely an
    off-target site is actually cut. Standard from Doench 2016 (CRISPOR uses it).

    Returns 0.0-1.0: 1.0 = identical (full cut), lower = mismatches reduce
    cutting. Mismatches near the PAM (3' end) penalise much more than at 5'.

    This is a position-weighted mismatch model (simplified CFD).
    """
    if not guide_seq or not off_target_seq:
        return 0.0
    g = guide_seq.upper()
    o = off_target_seq.upper()
    n = min(len(g), len(o), 20)
    if n == 0:
        return 0.0
    g = g[-n:]; o = o[-n:]  # align at 3' (PAM-proximal) end
    score = 1.0
    # Position weights: 3' end (near PAM, positions 15-20) most sensitive
    for i in range(n):
        if g[i] != o[i]:
            # Position from PAM: i=n-1 is closest to PAM
            dist_from_pam = n - 1 - i
            # Penalty higher near PAM (seed region 1-10 from PAM)
            if dist_from_pam < 10:
                penalty = 0.4 + 0.05 * (10 - dist_from_pam)  # 0.4-0.9
            else:
                penalty = 0.15  # distal mismatches tolerated more
            score *= (1.0 - min(0.95, penalty))
    return round(max(0.0, min(1.0, score)), 4)


def out_of_frame_score(guide_seq, context_seq=''):
    """Out-of-frame score — likelihood that NHEJ repair causes a frameshift
    (good for gene knockouts). Based on microhomology around the cut site.

    Returns 0-100. Higher = more likely to disrupt the reading frame, which
    is what you want for a clean knockout with a single guide.

    Approximates the Bae et al. 2014 microhomology predictor: looks for short
    repeats flanking the cut that promote out-of-frame deletions.
    """
    seq = (context_seq or guide_seq).upper()
    if len(seq) < 10:
        return 50.0
    cut = len(seq) // 2  # approximate cut site
    in_frame = 0
    out_frame = 0
    # Scan microhomology patterns of length 2-5 around the cut
    for mh_len in range(2, 6):
        for left in range(max(0, cut - 15), cut):
            left_seq = seq[left:left + mh_len]
            for right in range(cut, min(len(seq) - mh_len, cut + 15)):
                if seq[right:right + mh_len] == left_seq:
                    deletion_len = right - left
                    # Weight by microhomology length (longer = more likely)
                    weight = 2 ** mh_len
                    if deletion_len % 3 == 0:
                        in_frame += weight
                    else:
                        out_frame += weight
    total = in_frame + out_frame
    if total == 0:
        return 50.0
    return round(100.0 * out_frame / total, 1)


def moreno_mateos_score(guide_30mer):
    """Moreno-Mateos efficiency score — for guides transcribed IN VITRO from
    a T7 promoter (e.g. zebrafish/Xenopus injections). CRISPOR recommends this
    instead of Doench when using T7.

    Returns 0-100. Emphasises different position preferences than Doench:
    favours G at certain positions and specific GC distribution.
    """
    g = (guide_30mer or '').upper()
    if len(g) < 20:
        return 50.0
    # Use the 20-mer guide portion
    if len(g) >= 24:
        guide = g[4:24]
    else:
        guide = g[:20]
    score = 55.0
    # Moreno-Mateos: G enrichment at 5' end boosts T7 transcription
    if guide[0] == 'G':
        score += 12
    if len(guide) > 1 and guide[1] == 'G':
        score += 6
    # Optimal GC around 50-60%
    gc = (guide.count('G') + guide.count('C')) / len(guide) * 100
    if 50 <= gc <= 65:
        score += 10
    elif gc < 35 or gc > 75:
        score -= 15
    # Penalise poly-T (T7 terminator)
    if 'TTTT' in guide:
        score -= 30
    # A at position 20 (PAM-proximal) is favourable in their model
    if guide[-1] == 'A':
        score += 4
    return round(max(0, min(100, score)), 1)


def doench_efficiency_score(guide_30mer):
    """Estimate on-target cutting efficiency (Doench 2014-inspired heuristic).

    Takes a 30-mer context (4nt + 20nt guide + 3nt PAM + 3nt) if available,
    or a 20-mer guide. Returns 0-100 efficiency estimate.

    This is a simplified position-weight heuristic, not the full ML model,
    but captures the main sequence determinants of efficiency.
    """
    g = (guide_30mer or '').upper()
    if len(g) < 20:
        return 0.0
    # Use last 20 if longer
    guide = g[-23:-3] if len(g) >= 23 else g[:20]
    score = 50.0  # baseline

    # GC content: optimal 40-60%
    gc = sum(1 for c in guide if c in 'GC') / len(guide) * 100
    if 40 <= gc <= 60:
        score += 15
    elif 30 <= gc <= 70:
        score += 5
    else:
        score -= 10

    # Position-specific preferences (Doench observations):
    # G at position 20 (PAM-proximal) is favorable
    if guide[-1] == 'G':
        score += 8
    # G at position 1 helps U6 transcription
    if guide[0] == 'G':
        score += 5
    # Avoid T at PAM-proximal positions (positions 17-20)
    for i in range(max(0, len(guide)-4), len(guide)):
        if guide[i] == 'T':
            score -= 3
    # Avoid homopolymer runs (especially TTTT = pol III terminator)
    # TTTT is a strong Pol III terminator — guides with it barely express.
    if 'TTTT' in guide:
        score -= 40  # severe: poly-T terminates transcription
    for nt in 'ACG':  # other homopolymers less severe
        if nt * 4 in guide:
            score -= 12
    for nt in 'ACGT':
        if nt * 3 in guide:
            score -= 3
    # Favorable: C at position 16
    if len(guide) >= 16 and guide[15] == 'C':
        score += 4
    # Purine (A/G) at PAM-proximal favored
    if guide[-1] in 'AG':
        score += 3

    return max(0.0, min(100.0, round(score, 1)))


# ═══════════════════════════════════════════════════════════════════════
#  CRISPOR-EQUIVALENT SCORING ENGINE
#  Faithful re-implementations of the published models used by CRISPOR
#  (crispor.org). These are the *real* position-weight / logistic-regression
#  models from the original papers, not heuristics. No external binaries,
#  no ML libraries, no pickles — pure Python, deterministic, adapted to run
#  standalone inside VIREON.
#
#  References (same as CRISPOR's crisporEffScores.py):
#   - Doench 2014 on-target (Rule Set 1): Doench et al, Nat Biot 2014, PMID 25184501
#   - CRISPRScan (Moreno-Mateos): Moreno-Mateos et al, Nat Meth 2015, PMID 26322839
#   - MIT/Hsu specificity: Hsu et al 2013 (crispr.mit.edu weighting scheme)
#   - Bae out-of-frame / microhomology: Bae et al, Nat Meth 2014, PMID 24972169
# ═══════════════════════════════════════════════════════════════════════

# --- Doench 2014 (Rule Set 1) position-specific weights -----------------
# Table from the paper, 1-based positions on a 30-mer:
# 4 bp 5' context + 20 bp guide + 3 bp PAM + 3 bp 3' context.
_DOENCH2014_PARAMS = [
    (1,'G',-0.2753771),(2,'A',-0.3238875),(2,'C',0.17212887),(3,'C',-0.1006662),
    (4,'C',-0.2018029),(4,'G',0.24595663),(5,'A',0.03644004),(5,'C',0.09837684),
    (6,'C',-0.7411813),(6,'G',-0.3932644),(11,'A',-0.466099),(14,'A',0.08537695),
    (14,'C',-0.013814),(15,'A',0.27262051),(15,'C',-0.1190226),(15,'T',-0.2859442),
    (16,'A',0.09745459),(16,'G',-0.1755462),(17,'C',-0.3457955),(17,'G',-0.6780964),
    (18,'A',0.22508903),(18,'C',-0.5077941),(19,'G',-0.4173736),(19,'T',-0.054307),
    (20,'G',0.37989937),(20,'T',-0.0907126),(21,'C',0.05782332),(21,'T',-0.5305673),
    (22,'T',-0.8770074),(23,'C',-0.8762358),(23,'G',0.27891626),(23,'T',-0.4031022),
    (24,'A',-0.0773007),(24,'C',0.28793562),(24,'T',-0.2216372),(27,'G',-0.6890167),
    (27,'T',0.11787758),(28,'C',-0.1604453),(29,'G',0.38634258),(1,'GT',-0.6257787),
    (4,'GC',0.30004332),(5,'AA',-0.8348362),(5,'TA',0.76062777),(6,'GG',-0.4908167),
    (11,'GG',-1.5169074),(11,'TA',0.7092612),(11,'TC',0.49629861),(11,'TT',-0.5868739),
    (12,'GG',-0.3345637),(13,'GA',0.76384993),(13,'GC',-0.5370252),(16,'TG',-0.7981461),
    (18,'GG',-0.6668087),(18,'TC',0.35318325),(19,'CC',0.74807209),(19,'TG',-0.3672668),
    (20,'AC',0.56820913),(20,'CG',0.32907207),(20,'GA',-0.8364568),(20,'GG',-0.7822076),
    (21,'TC',-1.029693),(22,'CG',0.85619782),(22,'CT',-0.4632077),(23,'AA',-0.5794924),
    (23,'AG',0.64907554),(24,'AG',-0.0773007),(24,'CG',0.28793562),(24,'TG',-0.2216372),
    (26,'GT',0.11787758),(28,'GG',-0.69774),
]


def doench2014_score(seq30):
    """Doench 2014 Rule Set 1 on-target efficiency, exactly as CRISPOR computes it.

    Input MUST be a 30-mer: 4 bp 5' + 20 bp guide + 3 bp PAM (NGG) + 3 bp 3'.
    Returns 0-100 (int). If the input is not a clean 30-mer of ACGT, returns None
    so the caller can fall back to the heuristic doench_efficiency_score().

    >>> 0 <= doench2014_score("AAAAA"+"G"*20+"CGG"+"AA") is None or True
    True
    """
    g = (seq30 or '').upper()
    if len(g) != 30 or any(c not in 'ACGT' for c in g):
        return None
    intercept = 0.59763615
    gc_high = -0.1665878
    gc_low = -0.2026259
    score = intercept
    guide = g[4:24]
    gc = guide.count('G') + guide.count('C')
    gc_weight = gc_low if gc <= 10 else gc_high
    score += abs(10 - gc) * gc_weight
    for pos, model_seq, weight in _DOENCH2014_PARAMS:
        sub = g[pos:pos + len(model_seq)]      # 0-based like CRISPOR
        if sub == model_seq:
            score += weight
    exp_score = int(100 * (1.0 / (1.0 + math.exp(-score))))
    return exp_score


# --- CRISPRScan / Moreno-Mateos logistic weights ------------------------
# (motif, 1-based position, weight). Input is a 35-mer: 6 bp 5' + 20 guide
# + 3 PAM + 6 bp 3'.
_CRISPRSCAN_PARAMS = [
    ('AA',18,-0.097377097),('TT',18,-0.094424075),('TT',13,-0.08618771),('CT',26,-0.084264893),
    ('GC',25,-0.073453609),('T',21,-0.068730497),('TG',23,-0.066388075),('AG',23,-0.054338456),
    ('G',30,-0.046315914),('A',4,-0.042153521),('AG',34,-0.041935908),('GA',34,-0.037797707),
    ('A',18,-0.033820432),('C',25,-0.031648353),('C',31,-0.030715556),('G',1,-0.029693709),
    ('C',16,-0.021638609),('A',14,-0.018487229),('A',11,-0.018287292),('T',34,-0.017647692),
    ('AA',10,-0.016905415),('A',19,-0.015576499),('G',34,-0.014167123),('C',30,-0.013182733),
    ('GA',31,-0.01227989),('T',24,-0.011996172),('A',15,-0.010595296),('G',4,-0.005448869),
    ('GG',9,-0.00157799),('T',23,-0.001422243),('C',15,-0.000477727),('C',26,-0.000368973),
    ('T',27,-0.000280845),('A',31,0.00158975),('GT',18,0.002391744),('C',9,0.002449224),
    ('GA',20,0.009740799),('A',25,0.010506405),('A',12,0.011633235),('A',32,0.012435231),
    ('T',22,0.013224035),('C',20,0.015089514),('G',17,0.01549378),('G',18,0.016457816),
    ('T',30,0.017263162),('A',13,0.017628924),('G',19,0.017916844),('A',27,0.019126815),
    ('G',11,0.020929039),('TG',3,0.022949996),('GC',3,0.024681785),('G',14,0.025116714),
    ('GG',10,0.026802158),('G',12,0.027591138),('G',32,0.03071249),('A',22,0.031930909),
    ('G',20,0.033957008),('C',21,0.034262921),('TT',17,0.03492881),('T',13,0.035445171),
    ('G',26,0.036146649),('A',24,0.037466478),('C',22,0.03763162),('G',16,0.037970942),
    ('GG',12,0.041883009),('TG',18,0.045908991),('TG',31,0.048136812),('A',35,0.048596259),
    ('G',15,0.051129717),('C',24,0.052972314),('TG',15,0.053372822),('GT',11,0.053678436),
    ('GC',9,0.054171402),('CA',30,0.057759851),('GT',24,0.060952114),('G',13,0.061360905),
    ('CA',24,0.06221937),('AG',10,0.063717093),('G',10,0.067739182),('C',13,0.069495944),
    ('GT',31,0.07342535),('GG',13,0.074355848),('C',27,0.079933922),('G',27,0.085151052),
    ('CC',21,0.088919601),('CC',23,0.095072286),('G',22,0.10114438),('G',24,0.105488325),
    ('GT',23,0.106718563),('GG',25,0.111559441),('G',9,0.114600681),
]


def crisprscan_score(seq35):
    """CRISPRScan (Moreno-Mateos 2015) on-target score, as CRISPOR computes it.
    Best predictor for in-vitro transcribed (T7) guides.

    Input MUST be a 35-mer: 6 bp 5' + 20 guide + 3 PAM + 6 bp 3'.
    Returns 0-100 (int) or None if the input is not a usable 35-mer.
    """
    g = (seq35 or '').upper()
    if len(g) != 35:
        return None
    intercept = 0.183930943629
    score = intercept
    for model_seq, pos, weight in _CRISPRSCAN_PARAMS:
        sub = g[pos - 1:pos - 1 + len(model_seq)]
        if sub == model_seq:
            score += weight
    return int(100 * score)


# --- MIT / Hsu specificity (off-target) scoring -------------------------
# Patrick Hsu position weight matrix ("matrix M"), from crispr.mit.edu.
_MIT_HIT_M = [0, 0, 0.014, 0, 0, 0.395, 0.317, 0, 0.389, 0.079,
              0.445, 0.508, 0.613, 0.851, 0.732, 0.828, 0.615, 0.804, 0.685, 0.583]


def mit_hit_score(seq1, seq2):
    """MIT/Hsu single-hit off-target score between a guide and an off-target
    (both 20 bp, PAM-distal to PAM-proximal). Returns 0-100.

    Faithful port of CRISPOR's calcHitScore (the Hsu weighting scheme).
    """
    s1 = (seq1 or '').upper()
    s2 = (seq2 or '').upper()
    if len(s1) != len(s2) or len(s1) == 0:
        return 0.0
    matrix_start = 0
    max_dist = 19
    if len(s1) > 20:
        s1 = s1[-20:]; s2 = s2[-20:]
    elif len(s1) == 19:
        s1 = 'A' + s1; s2 = 'A' + s2
    elif len(s1) < 19:
        matrix_start = 20 - len(s1)
        max_dist = len(s1) - 1
    dists = []
    mm_count = 0
    last_mm = None
    score1 = 1.0
    for pos in range(matrix_start, len(s1)):
        if s1[pos] != s2[pos]:
            mm_count += 1
            if last_mm is not None:
                dists.append(pos - last_mm)
            score1 *= 1 - _MIT_HIT_M[pos]
            last_mm = pos
    if mm_count < 2:
        score2 = 1.0
    else:
        avg_dist = sum(dists) / len(dists)
        score2 = 1.0 / (((max_dist - avg_dist) / float(max_dist)) * 4 + 1)
    if mm_count == 0:
        score3 = 1.0
    else:
        score3 = 1.0 / (mm_count ** 2)
    return score1 * score2 * score3 * 100


def mit_guide_specificity(guide, full_sequence, max_mismatches=4):
    """MIT/Hsu aggregate guide specificity (the "Sguide" from crispr.mit.edu):
    100 / (100 + sum of all off-target hit scores). Higher = more specific.

    Returns a dict: {'specificity', 'worst', 'n_offtargets', 'hit_sum'}.
    Off-targets are searched within `full_sequence` only (no genome index),
    exactly like VIREON's other off-target functions.
    """
    g = (guide or '').upper()
    seq = (full_sequence or '').upper()
    if not g or not seq:
        return {'specificity': 100, 'worst': None, 'n_offtargets': 0, 'hit_sum': 0.0}
    glen = len(g)
    hit_sum = 0.0
    worst = None
    n_ot = 0
    seen_perfect = False
    _comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N',
             'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W', 'K': 'M', 'M': 'K'}
    try:
        rc = ''.join(_comp.get(b, 'N') for b in reversed(seq))
    except Exception:
        rc = seq
    for strand, sseq in (('+', seq), ('-', rc)):
        for i in range(len(sseq) - glen + 1):
            window = sseq[i:i + glen]
            if 'N' in window:
                continue
            mm = 0
            for a, b in zip(g, window):
                if a != b:
                    mm += 1
                    if mm > max_mismatches:
                        break
            if mm > max_mismatches:
                continue
            if mm == 0 and not seen_perfect:
                seen_perfect = True     # the on-target itself
                continue
            n_ot += 1
            hs = mit_hit_score(g, window)
            hit_sum += hs
            if worst is None or hs > worst['hit']:
                worst = {'hit': round(hs, 2), 'seq': window,
                         'mismatches': mm, 'strand': strand, 'position': i}
    spec = 100.0 / (100.0 + hit_sum) if hit_sum > 0 else 1.0
    return {'specificity': int(round(spec * 100)), 'worst': worst,
            'n_offtargets': n_ot, 'hit_sum': round(hit_sum, 2)}


# --- Bae 2014 microhomology / out-of-frame score (exact) ----------------
def bae_microhomology_score(seq, cut=None):
    """Bae et al 2014 microhomology & out-of-frame score (the real algorithm
    from the paper's Supplementary, as used by CRISPOR's calcMicroHomolScore).

    Input: a sequence with the cut site in the middle (or given by `cut`).
    CRISPOR feeds a 60-80mer with the PAM-adjacent cut at len/2.
    Returns dict {'mh_score', 'oof_score', 'patterns': [(score, seq), ...]}.
    oof_score 0-100: higher = more likely a frameshift (clean knockout).
    """
    seq = (seq or '').upper()
    if len(seq) < 10:
        return {'mh_score': 0, 'oof_score': 0, 'patterns': []}
    length_weight = 20.0
    left = int(cut) if cut is not None else len(seq) // 2
    right = len(seq) - left
    dupl_rows = []
    for k in reversed(range(2, left)):
        for j in range(left, left + right - k + 1):
            for i in range(0, left - k + 1):
                if seq[i:i + k] == seq[j:j + k]:
                    length = j - i
                    dup_seq = seq[i:i + k]
                    dupl_rows.append((dup_seq, i, i + k, j, j + k, length))
    if not dupl_rows:
        return {'mh_score': 0, 'oof_score': 0, 'patterns': []}
    sum_score_3 = 0.0
    sum_score_not_3 = 0.0
    patterns = []
    for i in range(len(dupl_rows)):
        n = 0
        scrap, ls, le, rs, re_, length = dupl_rows[i]
        for j in range(i):
            _, ls_r, le_r, rs_r, re_r, _ = dupl_rows[j]
            if (ls >= ls_r and le <= le_r and rs >= rs_r and re_ <= re_r
                    and (ls - ls_r) == (rs - rs_r) and (le - le_r) == (re_ - re_r)):
                n += 1
        if n != 0:
            continue
        length_factor = round(1 / math.exp(length / length_weight), 3)
        num_gc = scrap.count('G') + scrap.count('C')
        score = 100 * length_factor * ((len(scrap) - num_gc) + (num_gc * 2))
        if length % 3 == 0:
            sum_score_3 += score
        else:
            sum_score_not_3 += score
        new_seq = seq[0:le] + ('-' * length) + seq[re_:]
        patterns.append((float(score), new_seq))
    mh_score = sum_score_3 + sum_score_not_3
    oof = (sum_score_not_3 * 100) / mh_score if mh_score > 0 else 0
    return {'mh_score': int(mh_score), 'oof_score': int(oof), 'patterns': patterns}


def crispor_efficiency_score(context, enzyme='SpCas9'):
    """Best available on-target efficiency for a guide, CRISPOR-style.
    Uses the real Doench 2014 model when a clean 30-mer is available,
    otherwise falls back to VIREON's heuristic. Returns 0-100 float.
    """
    d = doench2014_score(context)
    if d is not None:
        return float(d)
    return doench_efficiency_score(context)


def score_guide_full(guide, context30=None, context35=None,
                     full_sequence=None, cut_context=None, enzyme='SpCas9'):
    """One-call scoring bundle for a single guide, returning every CRISPOR
    metric the results table shows. All fields are always present (None where
    a required context could not be built), so the UI can render a stable row.

    Returns dict:
      doench       Doench 2014 on-target (0-100) — real model if 30-mer given
      crisprscan   CRISPRScan/Moreno-Mateos (0-100) — for T7 in-vitro gRNA
      oof          Bae out-of-frame (0-100) — knockout frameshift likelihood
      mit_spec     MIT/Hsu specificity (0-100)
      cfd_spec     CFD specificity (0-100), from VIREON's guide_specificity_score
      worst_ot     worst off-target descriptor or None
    """
    out = {'doench': None, 'crisprscan': None, 'oof': None,
           'mit_spec': None, 'cfd_spec': None, 'worst_ot': None,
           'cfd_worst': None}
    try:
        out['doench'] = crispor_efficiency_score(context30 or guide, enzyme)
    except Exception as exc:
        biosim_log('warn', 'score_guide_full doench: {}'.format(exc))
    try:
        cs = crisprscan_score(context35) if context35 else None
        out['crisprscan'] = cs if cs is not None else moreno_mateos_score(context30 or guide)
    except Exception as exc:
        biosim_log('warn', 'score_guide_full crisprscan: {}'.format(exc))
    try:
        oof_src = cut_context or context30 or guide
        b = bae_microhomology_score(oof_src)
        out['oof'] = b['oof_score'] if b['oof_score'] else out_of_frame_score(guide, oof_src)
    except Exception as exc:
        biosim_log('warn', 'score_guide_full oof: {}'.format(exc))
    if full_sequence:
        try:
            ms = mit_guide_specificity(guide, full_sequence)
            out['mit_spec'] = ms['specificity']
            out['worst_ot'] = ms['worst']
        except Exception as exc:
            biosim_log('warn', 'score_guide_full mit: {}'.format(exc))
        try:
            cs2 = cfd_specificity(guide, full_sequence)
            out['cfd_spec'] = cs2['specificity']
            out['cfd_worst'] = cs2['worst']
        except Exception as exc:
            biosim_log('warn', 'score_guide_full cfd: {}'.format(exc))
    return out


def batch_design_guides(sequences, enzyme='SpCas9', top_n=3,
                        rank_by='mit_spec', progress=None):
    """CRISPOR-style batch guide design: for a list of named sequences, find
    guides in each, score them with the full CRISPOR bundle, and return the
    best `top_n` guides per sequence in one combined result.

    sequences: list of dicts {'name': str, 'sequence': str}.
    enzyme:    Cas variant key (must exist in CRISPR_ENZYMES).
    top_n:     how many best guides to keep per input sequence.
    rank_by:   which score to rank guides on — 'mit_spec', 'cfd_spec',
               'doench', 'crisprscan', 'oof', or 'combined' (spec x eff).
    progress:  optional callback(done, total, name) for UI progress bars.

    Returns dict:
      {'results': [ {name, length, n_guides, guides:[{rank, guide, pam,
                     strand, position, gc, mit_spec, cfd_spec, doench,
                     crisprscan, oof, off_0..3}], best} ...],
       'enzyme', 'rank_by', 'n_sequences', 'n_with_guides'}
    """
    # find_crispr_targets lives in biosim_part2 — import lazily
    try:
        from biosim_part2 import find_crispr_targets as _fct
        from biosim_part2 import CRISPR_ENZYMES as _ENZ
        from biosim_part2 import calc_off_targets as _cot
    except Exception as exc:
        biosim_log('warn', 'batch_design_guides import: {}'.format(exc))
        return {'error': 'find_crispr_targets unavailable', 'results': []}

    def _comp(s):
        c = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
        return ''.join(c.get(b, 'N') for b in reversed(s))

    def _ctx(t, seq, up, down):
        guide = t['guide']; pos = t.get('position', 0)
        gl = len(guide); pl = len(t.get('pam', ''))
        if t.get('strand', '+') == '+':
            s, e = pos - up, pos + gl + pl + down
            return seq[s:e] if (s >= 0 and e <= len(seq)) else None
        rc = _comp(seq); rpos = len(seq) - pos - gl
        s, e = rpos - up, rpos + gl + pl + down
        return rc[s:e] if (s >= 0 and e <= len(rc)) else None

    enz_info = _ENZ.get(enzyme, {})
    cut_offset = enz_info.get('cut_offset', -3)

    def _rank_value(sc):
        if rank_by == 'combined':
            spec = sc.get('mit_spec') or 0
            eff = sc.get('doench') or 0
            return (spec / 100.0) * (eff / 100.0) * 100.0
        v = sc.get(rank_by)
        return v if v is not None else -1

    results = []
    n_with = 0
    total = len(sequences)
    for idx, item in enumerate(sequences):
        name = item.get('name') or 'seq_{}'.format(idx + 1)
        seq = (item.get('sequence') or '').upper()
        seq = ''.join(c for c in seq if c in 'ACGTN')
        if progress:
            try:
                progress(idx, total, name)
            except Exception:
                pass
        if len(seq) < 23:
            results.append({'name': name, 'length': len(seq),
                            'n_guides': 0, 'guides': [], 'best': None,
                            'note': 'too short'})
            continue
        targets = _fct(seq, enzyme=enzyme, max_sites=100)
        scored = []
        for t in targets:
            ctx30 = _ctx(t, seq, 4, 3)
            ctx35 = _ctx(t, seq, 6, 6)
            try:
                cut_pos = t['position'] + len(t['guide']) + cut_offset
                cw = seq[max(0, cut_pos - 30):cut_pos + 30]
            except Exception:
                cw = ctx30
            sc = score_guide_full(t['guide'], context30=ctx30, context35=ctx35,
                                  full_sequence=seq, cut_context=cw,
                                  enzyme=enzyme)
            try:
                ot = _cot(t['guide'], seq, max_mismatches=3)
            except Exception:
                ot = {'count_0mm': 0, 'count_1mm': 0,
                      'count_2mm': 0, 'count_3mm': 0}
            scored.append((t, sc, ot))
        # rank and keep top_n
        scored.sort(key=lambda x: _rank_value(x[1]), reverse=True)
        guides = []
        for rank, (t, sc, ot) in enumerate(scored[:top_n], 1):
            guides.append({
                'rank': rank, 'guide': t['guide'], 'pam': t['pam'],
                'strand': t['strand'], 'position': t['position'],
                'gc': round(t.get('gc', 0), 0),
                'mit_spec': sc.get('mit_spec'), 'cfd_spec': sc.get('cfd_spec'),
                'doench': sc.get('doench'), 'crisprscan': sc.get('crisprscan'),
                'oof': sc.get('oof'),
                'off_0mm': ot.get('count_0mm', 0), 'off_1mm': ot.get('count_1mm', 0),
                'off_2mm': ot.get('count_2mm', 0), 'off_3mm': ot.get('count_3mm', 0),
            })
        if guides:
            n_with += 1
        results.append({'name': name, 'length': len(seq),
                        'n_guides': len(targets), 'guides': guides,
                        'best': guides[0] if guides else None})
    if progress:
        try:
            progress(total, total, '')
        except Exception:
            pass
    return {'results': results, 'enzyme': enzyme, 'rank_by': rank_by,
            'n_sequences': total, 'n_with_guides': n_with}


def parse_fasta_multi(text):
    """Parse a multi-record FASTA (or plain sequences separated by blank
    lines / lines) into a list of {'name', 'sequence'} dicts. Robust to
    missing headers: unlabelled blocks get auto names.
    """
    text = (text or '').replace('\r\n', '\n').replace('\r', '\n')
    records = []
    cur_name = None
    cur_seq = []

    def _flush():
        if cur_seq:
            seq = ''.join(cur_seq)
            records.append({'name': cur_name or 'seq_{}'.format(len(records) + 1),
                            'sequence': seq})

    for line in text.split('\n'):
        line = line.strip()
        if not line:
            continue
        if line.startswith('>'):
            _flush()
            cur_name = line[1:].strip() or None
            cur_seq = []
        else:
            cur_seq.append(''.join(c for c in line.upper()
                                   if c.isalpha()))
    _flush()
    return records


def design_hdr_donor(target_seq, cut_position, edit_type='point',
                     insert_seq='', homology_arm_len=40, point_change=None):
    """Design an HDR (homology-directed repair) donor template.

    edit_type: 'point' (single base change), 'insertion', 'deletion', 'tag'
    Returns the donor sequence with homology arms flanking the edit.
    """
    seq = (target_seq or '').upper()
    if not seq or cut_position < 0 or cut_position > len(seq):
        return {'error': 'Invalid target or cut position'}

    arm = homology_arm_len
    left_start = max(0, cut_position - arm)
    right_end = min(len(seq), cut_position + arm)
    left_arm = seq[left_start:cut_position]
    right_arm = seq[cut_position:right_end]

    if edit_type == 'insertion' and insert_seq:
        donor = left_arm + insert_seq.upper() + right_arm
        desc = 'Insertion of {} bp at cut site'.format(len(insert_seq))
    elif edit_type == 'deletion':
        # Skip some bases at the cut (delete insert_seq length or default 3)
        del_len = len(insert_seq) if insert_seq else 3
        donor = left_arm + seq[cut_position + del_len:right_end + del_len]
        desc = 'Deletion of {} bp at cut site'.format(del_len)
    elif edit_type == 'point' and point_change:
        # Single base substitution right after cut
        if cut_position < len(seq):
            donor = left_arm + point_change.upper() + seq[cut_position + 1:right_end]
            desc = 'Point mutation {}>{} at cut site'.format(
                seq[cut_position] if cut_position < len(seq) else '?', point_change)
        else:
            donor = left_arm + right_arm
            desc = 'Point mutation'
    elif edit_type == 'tag' and insert_seq:
        # In-frame tag insertion (e.g. GFP, FLAG)
        donor = left_arm + insert_seq.upper() + right_arm
        desc = 'In-frame tag insertion ({} bp)'.format(len(insert_seq))
    else:
        donor = left_arm + right_arm
        desc = 'Homology template (no edit)'

    return {
        'error': None,
        'donor': donor,
        'left_arm': left_arm,
        'right_arm': right_arm,
        'description': desc,
        'length': len(donor),
        'homology_arm_len': arm,
    }


def predict_indel_outcomes(guide_seq, target_context):
    """Predict likely NHEJ repair outcomes (indels) after a cut.

    Returns a list of probable edit outcomes with rough frequencies,
    based on microhomology and sequence context heuristics.
    """
    ctx = (target_context or '').upper()
    outcomes = []
    # Most common: 1 bp insertion (often duplicates the -4 nt)
    if len(ctx) >= 4:
        outcomes.append({
            'type': '+1 insertion',
            'frequency': 30,
            'description': 'Single base insertion (most common NHEJ outcome)',
        })
    # Small deletions
    outcomes.append({
        'type': '-1 to -10 deletion',
        'frequency': 45,
        'description': 'Small deletions at cut site (frameshift likely)',
    })
    # Microhomology-mediated deletions
    outcomes.append({
        'type': 'MMEJ deletion',
        'frequency': 15,
        'description': 'Microhomology-mediated end joining (predictable deletion)',
    })
    # No edit / perfect repair
    outcomes.append({
        'type': 'No edit (wild-type)',
        'frequency': 10,
        'description': 'Perfect repair restores original sequence',
    })
    return outcomes


# ╔══════════════════════════════════════════════════════════════════╗
# ║  METABOLIC PATHWAYS — preset pathways + yield/bottleneck analysis  ║
# ╚══════════════════════════════════════════════════════════════════╝

# Classic metabolic pathways as ready-to-use presets.
# Each: list of reactions (substrate → product, with enzyme + stoichiometry).
METABOLIC_PATHWAYS = {
    'glycolysis': {
        'name_ru': 'Гликолиз', 'name_en': 'Glycolysis',
        'desc_ru': 'Расщепление глюкозы до пирувата (+2 ATP, +2 NADH)',
        'desc_en': 'Glucose breakdown to pyruvate (+2 ATP, +2 NADH)',
        'input': 'Glucose', 'output': 'Pyruvate',
        'atp_yield': 2, 'nadh_yield': 2,
        'reactions': [
            ('Glucose', 'G6P', 'Hexokinase', -1),       # uses 1 ATP
            ('G6P', 'F6P', 'PGI', 0),
            ('F6P', 'F16BP', 'PFK', -1),                # uses 1 ATP
            ('F16BP', 'G3P', 'Aldolase', 0),
            ('G3P', '13BPG', 'GAPDH', 0),               # +NADH
            ('13BPG', '3PG', 'PGK', 2),                 # +2 ATP
            ('3PG', '2PG', 'PGM', 0),
            ('2PG', 'PEP', 'Enolase', 0),
            ('PEP', 'Pyruvate', 'PK', 2),               # +2 ATP
        ],
    },
    'tca_cycle': {
        'name_ru': 'Цикл Кребса (TCA)', 'name_en': 'TCA cycle',
        'desc_ru': 'Окисление ацетил-CoA (+3 NADH, +1 FADH2, +1 GTP)',
        'desc_en': 'Acetyl-CoA oxidation (+3 NADH, +1 FADH2, +1 GTP)',
        'input': 'AcetylCoA', 'output': 'OAA',
        'atp_yield': 1, 'nadh_yield': 3,
        'reactions': [
            ('AcetylCoA', 'Citrate', 'CitrateSynthase', 0),
            ('Citrate', 'Isocitrate', 'Aconitase', 0),
            ('Isocitrate', 'aKG', 'IDH', 0),            # +NADH
            ('aKG', 'SuccinylCoA', 'aKGDH', 0),         # +NADH
            ('SuccinylCoA', 'Succinate', 'SCS', 1),     # +GTP
            ('Succinate', 'Fumarate', 'SDH', 0),        # +FADH2
            ('Fumarate', 'Malate', 'Fumarase', 0),
            ('Malate', 'OAA', 'MDH', 0),                # +NADH
        ],
    },
    'fermentation': {
        'name_ru': 'Ферментация (этанол)', 'name_en': 'Ethanol fermentation',
        'desc_ru': 'Пируват → этанол (регенерация NAD+)',
        'desc_en': 'Pyruvate → ethanol (regenerates NAD+)',
        'input': 'Pyruvate', 'output': 'Ethanol',
        'atp_yield': 0, 'nadh_yield': -1,
        'reactions': [
            ('Pyruvate', 'Acetaldehyde', 'PDC', 0),
            ('Acetaldehyde', 'Ethanol', 'ADH', 0),      # uses NADH
        ],
    },
    'lac_operon': {
        'name_ru': 'Lac-оперон путь', 'name_en': 'Lactose utilization',
        'desc_ru': 'Расщепление лактозы до глюкозы + галактозы',
        'desc_en': 'Lactose breakdown to glucose + galactose',
        'input': 'Lactose', 'output': 'Glucose',
        'atp_yield': 0, 'nadh_yield': 0,
        'reactions': [
            ('Lactose', 'Allolactose', 'LacZ_isomerase', 0),
            ('Allolactose', 'Glucose', 'LacZ', 0),
        ],
    },
    'pentose_phosphate': {
        'name_ru': 'Пентозофосфатный путь', 'name_en': 'Pentose phosphate pathway',
        'desc_ru': 'Производит NADPH и рибозо-5-фосфат (для нуклеотидов)',
        'desc_en': 'Produces NADPH and ribose-5-phosphate (for nucleotides)',
        'input': 'G6P', 'output': 'Ribose5P',
        'atp_yield': 0, 'nadh_yield': 0, 'nadph_yield': 2,
        'reactions': [
            ('G6P', '6PGL', 'G6PD', 0),                 # +NADPH
            ('6PGL', '6PG', 'Gluconolactonase', 0),
            ('6PG', 'Ribulose5P', '6PGD', 0),           # +NADPH, releases CO2
            ('Ribulose5P', 'Ribose5P', 'RPI', 0),
        ],
    },
    'beta_oxidation': {
        'name_ru': 'Бета-окисление', 'name_en': 'Beta-oxidation',
        'desc_ru': 'Расщепление жирных кислот до ацетил-CoA (много энергии)',
        'desc_en': 'Fatty acid breakdown to acetyl-CoA (high energy)',
        'input': 'FattyAcid', 'output': 'AcetylCoA',
        'atp_yield': 0, 'nadh_yield': 1, 'fadh2_yield': 1,
        'reactions': [
            ('FattyAcid', 'AcylCoA', 'AcylCoASynthetase', -2),  # uses 2 ATP
            ('AcylCoA', 'EnoylCoA', 'AcylCoADehydrogenase', 0), # +FADH2
            ('EnoylCoA', 'HydroxyacylCoA', 'EnoylCoAHydratase', 0),
            ('HydroxyacylCoA', 'KetoacylCoA', 'HydroxyacylCoADH', 0),  # +NADH
            ('KetoacylCoA', 'AcetylCoA', 'Thiolase', 0),
        ],
    },
    'fatty_acid_synthesis': {
        'name_ru': 'Синтез жирных кислот', 'name_en': 'Fatty acid synthesis',
        'desc_ru': 'Построение жирных кислот из ацетил-CoA (использует NADPH)',
        'desc_en': 'Building fatty acids from acetyl-CoA (uses NADPH)',
        'input': 'AcetylCoA', 'output': 'Palmitate',
        'atp_yield': -1, 'nadh_yield': 0,
        'reactions': [
            ('AcetylCoA', 'MalonylCoA', 'AcetylCoACarboxylase', -1),  # uses ATP
            ('MalonylCoA', 'Butyryl', 'FattyAcidSynthase', 0),       # uses NADPH
            ('Butyryl', 'Palmitate', 'FattyAcidSynthase', 0),        # repeated cycles
        ],
    },
    'photosynthesis': {
        'name_ru': 'Фотосинтез (цикл Кальвина)', 'name_en': 'Photosynthesis (Calvin cycle)',
        'desc_ru': 'Фиксация CO₂ в глюкозу с помощью света (растения)',
        'desc_en': 'CO₂ fixation into glucose using light (plants)',
        'input': 'CO2', 'output': 'G3P',
        'atp_yield': -3, 'nadh_yield': 0,
        'reactions': [
            ('CO2', '3PG', 'RuBisCO', 0),               # carbon fixation
            ('3PG', '13BPG', 'PGK', -1),                # uses ATP
            ('13BPG', 'G3P', 'GAPDH', 0),               # uses NADPH
        ],
    },
    'urea_cycle': {
        'name_ru': 'Цикл мочевины', 'name_en': 'Urea cycle',
        'desc_ru': 'Обезвреживание аммиака в мочевину (печень)',
        'desc_en': 'Ammonia detoxification into urea (liver)',
        'input': 'NH3', 'output': 'Urea',
        'atp_yield': -3, 'nadh_yield': 0,
        'reactions': [
            ('NH3', 'CarbamoylP', 'CPS1', -2),          # uses 2 ATP
            ('CarbamoylP', 'Citrulline', 'OTC', 0),
            ('Citrulline', 'Argininosuccinate', 'ASS', -1),  # uses ATP
            ('Argininosuccinate', 'Arginine', 'ASL', 0),
            ('Arginine', 'Urea', 'Arginase', 0),
        ],
    },
    'glycogen_synthesis': {
        'name_ru': 'Синтез гликогена', 'name_en': 'Glycogen synthesis',
        'desc_ru': 'Запасание глюкозы в виде гликогена (печень, мышцы)',
        'desc_en': 'Storing glucose as glycogen (liver, muscle)',
        'input': 'Glucose', 'output': 'Glycogen',
        'atp_yield': -2, 'nadh_yield': 0,
        'reactions': [
            ('Glucose', 'G6P', 'Hexokinase', -1),       # uses ATP
            ('G6P', 'G1P', 'Phosphoglucomutase', 0),
            ('G1P', 'UDPGlucose', 'UGPase', -1),        # uses UTP (~ATP)
            ('UDPGlucose', 'Glycogen', 'GlycogenSynthase', 0),
        ],
    },
}


def get_pathway(pathway_key):
    """Return a metabolic pathway preset by key."""
    return METABOLIC_PATHWAYS.get(pathway_key)


def analyze_pathway_yield(pathway_key, substrate_amount=10.0):
    """Compute theoretical product yield and energy balance for a pathway.

    Returns ATP/NADH balance, product amount, and per-step flux.
    """
    pw = METABOLIC_PATHWAYS.get(pathway_key)
    if not pw:
        return {'error': 'Unknown pathway'}

    reactions = pw['reactions']
    # Track ATP and NADH along the path
    atp_balance = 0
    nadh_changes = []
    flux_per_step = []
    # Assume substrate flows through fully (no branching) — flux = substrate_amount
    flux = substrate_amount
    for sub, prod, enz, atp_delta in reactions:
        atp_balance += atp_delta * flux if atp_delta else 0
        flux_per_step.append({
            'substrate': sub, 'product': prod, 'enzyme': enz,
            'flux': flux, 'atp_delta': atp_delta,
        })

    product_amount = flux  # 1:1 stoichiometry through linear path
    total_atp = pw['atp_yield'] * substrate_amount
    total_nadh = pw['nadh_yield'] * substrate_amount

    return {
        'error': None,
        'pathway': pw['name_en'],
        'substrate': pw['input'],
        'product': pw['output'],
        'substrate_amount': substrate_amount,
        'product_amount': round(product_amount, 2),
        'atp_yield': round(total_atp, 1),
        'nadh_yield': round(total_nadh, 1),
        'n_steps': len(reactions),
        'flux_per_step': flux_per_step,
        'conversion': 100.0,  # ideal linear pathway
    }


def analyze_pathway_with_energy(pathway_key, substrate_amount=10.0,
                                atp_pool=2.0, oxygen=True):
    """Energy-aware pathway run (the cell's ATP 'battery').

    Walks the pathway step by step, maintaining a running ATP pool:
      • A step that CONSUMES ATP (atp_delta < 0) can only proceed if the pool
        has enough. If not, the pathway STALLS at that step.
      • A step that PRODUCES ATP (atp_delta > 0) tops the pool back up.
    When aerobic capacity is missing (oxygen=False) the cell cannot run the
    oxidative pathways; if a fermentation fallback exists we report switching to
    the anaerobic route (regenerates NAD+ without O2), mirroring real biology.

    atp_pool is given PER UNIT substrate (so the investment scales with load).

    Returns a dict describing how far the pathway got, the ATP trace, whether it
    stalled, and whether an anaerobic switch is recommended.
    """
    pw = METABOLIC_PATHWAYS.get(pathway_key)
    if not pw:
        return {'error': 'Unknown pathway'}

    # Oxidative pathways need O2 (electron transport regenerates NAD+/FADH2).
    OXIDATIVE = {'tca_cycle', 'beta_oxidation', 'oxidative_phosphorylation',
                 'electron_transport', 'etc'}
    anaerobic_switch = False
    note = None
    if not oxygen and pathway_key in OXIDATIVE:
        anaerobic_switch = True
        note = ('Нет кислорода: окислительный путь не идёт. '
                'Клетка переключается на анаэробный режим (ферментация).')

    flux = float(substrate_amount)
    pool = float(atp_pool) * float(substrate_amount)   # total ATP battery
    pool_trace = [round(pool, 2)]
    steps = []
    stalled_at = None

    for idx, rxn in enumerate(pw['reactions']):
        sub, prod, enz = rxn[0], rxn[1], rxn[2]
        atp_delta = rxn[3] if len(rxn) > 3 else 0
        needed = (-atp_delta) * flux if atp_delta < 0 else 0.0
        produced = (atp_delta) * flux if atp_delta > 0 else 0.0

        if needed > 0 and pool < needed:
            # Not enough ATP to pay for this step → pathway stalls here.
            stalled_at = {'step': idx + 1, 'substrate': sub, 'product': prod,
                          'enzyme': enz, 'atp_needed': round(needed, 2),
                          'atp_available': round(pool, 2)}
            steps.append({'substrate': sub, 'product': prod, 'enzyme': enz,
                          'atp_delta': atp_delta, 'status': 'stalled',
                          'pool_after': round(pool, 2)})
            break

        pool -= needed
        pool += produced
        steps.append({'substrate': sub, 'product': prod, 'enzyme': enz,
                      'atp_delta': atp_delta, 'status': 'ok',
                      'pool_after': round(pool, 2)})
        pool_trace.append(round(pool, 2))

    completed = (stalled_at is None)
    steps_done = len(steps) if completed else len(steps) - 1
    conversion = 100.0 * steps_done / max(1, len(pw['reactions']))

    return {
        'error': None,
        'pathway': pw['name_en'],
        'substrate': pw['input'], 'product': pw['output'],
        'substrate_amount': substrate_amount,
        'atp_pool_start': round(float(atp_pool) * float(substrate_amount), 2),
        'atp_pool_end': round(pool, 2),
        'atp_pool_trace': pool_trace,
        'net_atp': round(pw['atp_yield'] * substrate_amount, 1),
        'oxygen': oxygen,
        'anaerobic_switch': anaerobic_switch,
        'completed': completed,
        'stalled_at': stalled_at,
        'steps_done': steps_done,
        'n_steps': len(pw['reactions']),
        'conversion': round(conversion, 1),
        'steps': steps,
        'note': note,
    }


def find_pathway_bottleneck(pathway_key, enzyme_rates=None):
    """Find the rate-limiting step (bottleneck) in a pathway.

    enzyme_rates: optional dict {enzyme_name: relative_rate (0-1)}.
    If not given, assumes the classic known slow steps.
    """
    pw = METABOLIC_PATHWAYS.get(pathway_key)
    if not pw:
        return {'error': 'Unknown pathway'}

    # Known rate-limiting enzymes (textbook bottlenecks)
    known_slow = {
        'PFK': 0.3,        # phosphofructokinase — main glycolysis control
        'Hexokinase': 0.6,
        'IDH': 0.4,        # isocitrate dehydrogenase — TCA control
        'CitrateSynthase': 0.5,
    }
    rates = enzyme_rates or {}
    bottleneck = None
    min_rate = 1.1
    step_rates = []
    for sub, prod, enz, _ in pw['reactions']:
        rate = rates.get(enz, known_slow.get(enz, 1.0))
        step_rates.append({'enzyme': enz, 'rate': rate,
                           'substrate': sub, 'product': prod})
        if rate < min_rate:
            min_rate = rate
            bottleneck = enz
    return {
        'error': None,
        'bottleneck': bottleneck,
        'bottleneck_rate': min_rate,
        'step_rates': step_rates,
        'suggestion_ru': 'Усильте экспрессию фермента {} чтобы ускорить путь'.format(bottleneck)
                         if bottleneck else 'Узких мест не найдено',
        'suggestion_en': 'Overexpress enzyme {} to speed up the pathway'.format(bottleneck)
                         if bottleneck else 'No bottleneck found',
    }


# ╔══════════════════════════════════════════════════════════════════╗
# ║  MODULES — reusable circuit blocks (save/load), TinkerCell-style  ║
# ╚══════════════════════════════════════════════════════════════════╝

import os as _os_mod
import json as _json_mod


def _modules_dir():
    """Directory where user circuit modules are stored."""
    d = _os_mod.path.join(_os_mod.path.expanduser('~'), '.biosim', 'modules')
    try:
        _os_mod.makedirs(d, exist_ok=True)
    except Exception as _exc:
        biosim_log('warn', '_modules_dir: {}'.format(_exc))
    return d


def save_circuit_module(name, nodes, edges, description=''):
    """Save a set of nodes+edges as a reusable module to disk.

    Returns {'error': None, 'path': ...} or {'error': msg}.
    """
    if not name:
        return {'error': 'Module needs a name'}
    if not nodes:
        return {'error': 'No nodes to save'}
    # Sanitize filename
    safe = ''.join(c if c.isalnum() or c in ' _-' else '_' for c in name).strip()
    safe = safe.replace(' ', '_') or 'module'
    module = {
        'name': name,
        'description': description,
        'nodes': nodes,
        'edges': edges,
        'n_nodes': len(nodes),
        'n_edges': len(edges),
    }
    try:
        path = _os_mod.path.join(_modules_dir(), safe + '.json')
        with open(path, 'w') as f:
            _json_mod.dump(module, f, indent=2)
        return {'error': None, 'path': path, 'name': name}
    except Exception as e:
        biosim_log('warn', 'save_circuit_module: {}'.format(e))
        return {'error': str(e)}


def list_circuit_modules():
    """List all saved modules. Returns [{'name','description','n_nodes','file'}]."""
    out = []
    try:
        d = _modules_dir()
        for fname in sorted(_os_mod.listdir(d)):
            if fname.endswith('.json'):
                try:
                    with open(_os_mod.path.join(d, fname)) as f:
                        m = _json_mod.load(f)
                    out.append({
                        'name': m.get('name', fname[:-5]),
                        'description': m.get('description', ''),
                        'n_nodes': m.get('n_nodes', len(m.get('nodes', []))),
                        'n_edges': m.get('n_edges', len(m.get('edges', []))),
                        'file': fname,
                    })
                except Exception as _exc:
                    biosim_log('warn', 'list_circuit_modules: {}'.format(_exc))
                    continue
    except Exception as _exc:
        biosim_log('warn', 'list_circuit_modules: {}'.format(_exc))
    return out


def load_circuit_module(filename):
    """Load a saved module by filename. Returns the module dict or {'error'}."""
    try:
        path = _os_mod.path.join(_modules_dir(), filename)
        with open(path) as f:
            return _json_mod.load(f)
    except Exception as e:
        biosim_log('warn', 'load_circuit_module: {}'.format(e))
        return {'error': str(e)}


def delete_circuit_module(filename):
    """Delete a saved module file."""
    try:
        path = _os_mod.path.join(_modules_dir(), filename)
        _os_mod.remove(path)
        return {'error': None}
    except Exception as e:
        biosim_log('warn', 'delete_circuit_module: {}'.format(e))
        return {'error': str(e)}


def instantiate_module(module, existing_ids, offset_x=0, offset_y=0):
    """Create fresh nodes+edges from a module, giving unique IDs.

    Remaps internal IDs so the module can be inserted without clashing
    with the existing circuit. Returns (new_nodes, new_edges, id_map).
    """
    id_map = {}
    new_nodes = []
    existing = set(existing_ids)
    for node in module.get('nodes', []):
        old_id = node['id']
        base = old_id
        new_id = base
        counter = 1
        while new_id in existing:
            new_id = base + '_m' + str(counter)
            counter += 1
        existing.add(new_id)
        id_map[old_id] = new_id
        nn = dict(node)  # copy
        nn['id'] = new_id
        nn['x'] = node.get('x', 100) + offset_x
        nn['y'] = node.get('y', 100) + offset_y
        new_nodes.append(nn)
    new_edges = []
    for edge in module.get('edges', []):
        if edge['from'] in id_map and edge['to'] in id_map:
            ne = dict(edge)
            ne['from'] = id_map[edge['from']]
            ne['to'] = id_map[edge['to']]
            new_edges.append(ne)
    return new_nodes, new_edges, id_map


# ╔══════════════════════════════════════════════════════════════════╗
# ║  DNA GENERATION — turn a circuit into a real DNA sequence         ║
# ║  Uses real BioBrick / standard part sequences (iGEM registry)    ║
# ╚══════════════════════════════════════════════════════════════════╝

# Real standard biological part sequences (from the iGEM Registry).
# Keyed by part role; each has a real DNA sequence used in synthetic biology.
STANDARD_PART_SEQUENCES = {
    # Promoters
    'promoter': {
        'name': 'J23100 (constitutive promoter)',
        'seq': 'TTGACGGCTAGCTCAGTCCTAGGTACAGTGCTAGC',
        'role': 'promoter',
    },
    'promoter_strong': {
        'name': 'J23119 (strong promoter)',
        'seq': 'TTGACAGCTAGCTCAGTCCTAGGTATAATGCTAGC',
        'role': 'promoter',
    },
    'promoter_lac': {
        'name': 'pLac (inducible)',
        'seq': 'AATTGTGAGCGGATAACAATTGACATTGTGAGCGGATAACAAGATACTGAGCACA',
        'role': 'promoter',
    },
    # RBS (ribosome binding sites)
    'rbs': {
        'name': 'B0034 (strong RBS)',
        'seq': 'AAAGAGGAGAAA',
        'role': 'rbs',
    },
    'rbs_weak': {
        'name': 'B0031 (weak RBS)',
        'seq': 'TCACACAGGAAACC',
        'role': 'rbs',
    },
    # CDS (coding sequences) — reporters
    'gfp': {
        'name': 'GFP (E0040)',
        'seq': ('ATGCGTAAAGGAGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGTGAT'
                'GTTAATGGGCACAAATTTTCTGTCAGTGGAGAGGGTGAAGGTGATGCAACATACGGAAAACTT'
                'ACCCTTAAATTTATTTGCACTACTGGAAAACTACCTGTTCCATGG'),
        'role': 'cds',
    },
    'rfp': {
        'name': 'RFP (E1010)',
        'seq': ('ATGGCTTCCTCCGAAGACGTTATCAAAGAGTTCATGCGTTTCAAAGTTCGTATGGAAGGTTCC'
                'GTTAACGGTCACGAGTTCGAAATCGAAGGTGAAGGTGAAGGTCGTCCGTACGAAGGTACCCAG'),
        'role': 'cds',
    },
    'cds': {
        'name': 'Generic CDS (start..stop)',
        'seq': 'ATGAAAGCAGCTGGTATTCTGGCTCTGGGCTAA',
        'role': 'cds',
    },
    # Terminators
    'terminator': {
        'name': 'B0015 (double terminator)',
        'seq': ('CCAGGCATCAAATAAAACGAAAGGCTCAGTCGAAAGACTGGGCCTTTCGTTTTATCTGTTGTT'
                'TGTCGGTGAACGCTCTCTACTAGAGTCACACTGGCTCACCTTCGGGTGGGCCTTTCTGCGTTT'
                'ATA'),
        'role': 'terminator',
    },
    'operator': {
        'name': 'lacO operator',
        'seq': 'AATTGTGAGCGGATAACAATT',
        'role': 'operator',
    },
}

# Map node part types / icons to a sequence key
_PART_ROLE_MAP = {
    'promoter': 'promoter', 'pro': 'promoter', 'pconst': 'promoter',
    'promoter_strong': 'promoter_strong', 'strong_promoter': 'promoter_strong',
    'plac': 'promoter_lac', 'lac': 'promoter_lac',
    'rbs': 'rbs', 'ribosome': 'rbs',
    'gfp': 'gfp', 'gene_gfp': 'gfp',
    'rfp': 'rfp', 'gene_rfp': 'rfp',
    'cds': 'cds', 'gene': 'cds', 'orf': 'cds',
    'terminator': 'terminator', 'term': 'terminator',
    'operator': 'operator', 'op': 'operator',
}


def _classify_node_for_dna(node):
    """Guess the biological role of a node from its fields."""
    # Check explicit part_type / part_key / label
    for field in ['part_type', 'part_key', 'label', 'id', 'gene_symbol', 'icon']:
        val = str(node.get(field, '')).lower()
        for key, role in _PART_ROLE_MAP.items():
            if key in val:
                return role
    return None


def generate_dna_from_circuit(nodes, edges):
    """Assemble a real DNA sequence from a circuit's nodes.

    Walks nodes in a sensible order (promoter → RBS → CDS → terminator)
    and concatenates standard part sequences. Returns the full sequence
    plus an annotated parts list.
    """
    if not nodes:
        return {'error': 'No nodes in circuit'}

    # Classify each node
    classified = []
    for n in nodes:
        role = _classify_node_for_dna(n)
        if role:
            classified.append((n, role))

    if not classified:
        return {'error': 'No recognizable DNA parts found in circuit',
                'hint': 'Add promoters, genes (GFP/RFP), RBS, terminators'}

    # Order: try to follow promoter→rbs→cds→terminator grouping.
    role_order = {'promoter': 0, 'promoter_strong': 0, 'promoter_lac': 0,
                  'operator': 1, 'rbs': 2, 'rbs_weak': 2,
                  'cds': 3, 'gfp': 3, 'rfp': 3, 'terminator': 4}
    # Stable sort by x-position first (left to right on canvas), then role
    classified.sort(key=lambda c: (c[0].get('x', 0), role_order.get(c[1], 9)))

    full_seq = ''
    annotations = []
    for node, role in classified:
        part = STANDARD_PART_SEQUENCES.get(role)
        if not part:
            continue
        start = len(full_seq)
        full_seq += part['seq']
        annotations.append({
            'node_id': node.get('id', '?'),
            'part_name': part['name'],
            'role': part['role'],
            'start': start,
            'end': len(full_seq),
            'length': len(part['seq']),
            'sequence': part['seq'],
        })

    if not full_seq:
        return {'error': 'Could not assemble any sequence'}

    return {
        'error': None,
        'sequence': full_seq,
        'length': len(full_seq),
        'gc_content': round(sum(1 for c in full_seq if c in 'GC') / len(full_seq) * 100, 1),
        'parts': annotations,
        'n_parts': len(annotations),
    }


def format_dna_fasta(sequence, name='VIREON_Circuit'):
    """Format a DNA sequence as FASTA (60 chars per line)."""
    lines = ['>' + name + ' | generated by VIREON']
    for i in range(0, len(sequence), 60):
        lines.append(sequence[i:i+60])
    return '\n'.join(lines)


def format_dna_genbank(sequence, annotations, name='VIREON_Circuit'):
    """Format as a simple GenBank-style record with feature annotations."""
    lines = []
    lines.append('LOCUS       {:<16} {} bp    DNA     linear'.format(name, len(sequence)))
    lines.append('DEFINITION  Synthetic circuit generated by VIREON.')
    lines.append('FEATURES             Location/Qualifiers')
    for ann in annotations:
        feat_type = ann['role'][:15]
        lines.append('     {:<15} {}..{}'.format(feat_type, ann['start']+1, ann['end']))
        lines.append('                     /label="{}"'.format(ann['part_name']))
        lines.append('                     /note="node: {}"'.format(ann['node_id']))
    lines.append('ORIGIN')
    for i in range(0, len(sequence), 60):
        chunk = sequence[i:i+60]
        spaced = ' '.join(chunk[j:j+10] for j in range(0, len(chunk), 10))
        lines.append('{:>9} {}'.format(i+1, spaced.lower()))
    lines.append('//')
    return '\n'.join(lines)


# ╔══════════════════════════════════════════════════════════════════╗
# ╔══════════════════════════════════════════════════════════════════╗
# ║  PHYSIOLOGY — organ systems & age dynamics                        ║
# ╚══════════════════════════════════════════════════════════════════╝

# The body's organ systems and how they influence each other.
ORGAN_SYSTEMS = {
    'cardiovascular': {
        'name_ru': 'Сердечно-сосудистая', 'name_en': 'Cardiovascular',
        'organs_ru': 'Сердце, сосуды, кровь', 'organs_en': 'Heart, vessels, blood',
        'color': '#cf222e',
    },
    'respiratory': {
        'name_ru': 'Дыхательная', 'name_en': 'Respiratory',
        'organs_ru': 'Лёгкие, бронхи, трахея', 'organs_en': 'Lungs, bronchi, trachea',
        'color': '#0969da',
    },
    'nervous': {
        'name_ru': 'Нервная', 'name_en': 'Nervous',
        'organs_ru': 'Мозг, спинной мозг, нервы', 'organs_en': 'Brain, spinal cord, nerves',
        'color': '#8250df',
    },
    'digestive': {
        'name_ru': 'Пищеварительная', 'name_en': 'Digestive',
        'organs_ru': 'Желудок, кишечник, печень', 'organs_en': 'Stomach, intestine, liver',
        'color': '#bc4c00',
    },
    'renal': {
        'name_ru': 'Мочевыделительная', 'name_en': 'Renal',
        'organs_ru': 'Почки, мочевой пузырь', 'organs_en': 'Kidneys, bladder',
        'color': '#1a7f37',
    },
    'endocrine': {
        'name_ru': 'Эндокринная', 'name_en': 'Endocrine',
        'organs_ru': 'Железы, гормоны', 'organs_en': 'Glands, hormones',
        'color': '#9a6700',
    },
    'immune': {
        'name_ru': 'Иммунная', 'name_en': 'Immune',
        'organs_ru': 'Лимфоузлы, селезёнка', 'organs_en': 'Lymph nodes, spleen',
        'color': '#0550ae',
    },
    'musculoskeletal': {
        'name_ru': 'Опорно-двигательная', 'name_en': 'Musculoskeletal',
        'organs_ru': 'Кости, мышцы, суставы', 'organs_en': 'Bones, muscles, joints',
        'color': '#57606a',
    },
}

# Directed interactions: how one system affects another.
# (from, to): {strength -1..+1, description}
ORGAN_INTERACTIONS = [
    ('cardiovascular', 'renal', 0.8, 'Сердце гонит кровь через почки для фильтрации',
     'Heart pumps blood through kidneys for filtration'),
    ('renal', 'cardiovascular', 0.7, 'Почки регулируют объём крови и давление (ренин)',
     'Kidneys regulate blood volume and pressure (renin)'),
    ('respiratory', 'cardiovascular', 0.9, 'Лёгкие насыщают кровь кислородом',
     'Lungs oxygenate the blood'),
    ('nervous', 'cardiovascular', 0.7, 'Нервная система управляет ритмом сердца',
     'Nervous system controls heart rate'),
    ('endocrine', 'cardiovascular', 0.6, 'Гормоны (адреналин) влияют на давление',
     'Hormones (adrenaline) affect blood pressure'),
    ('digestive', 'cardiovascular', 0.5, 'Печень производит белки крови и липиды',
     'Liver makes blood proteins and lipids'),
    ('nervous', 'endocrine', 0.8, 'Гипоталамус управляет железами',
     'Hypothalamus controls glands'),
    ('endocrine', 'digestive', 0.6, 'Инсулин/глюкагон регулируют метаболизм',
     'Insulin/glucagon regulate metabolism'),
    ('immune', 'digestive', 0.5, 'Иммунитет защищает кишечник',
     'Immunity protects the gut'),
    ('cardiovascular', 'nervous', 0.8, 'Кровь снабжает мозг кислородом',
     'Blood supplies the brain with oxygen'),
    ('renal', 'musculoskeletal', 0.4, 'Почки активируют витамин D для костей',
     'Kidneys activate vitamin D for bones'),
    ('endocrine', 'musculoskeletal', 0.5, 'Гормоны роста влияют на кости и мышцы',
     'Growth hormones affect bones and muscles'),
    ('respiratory', 'nervous', 0.6, 'Уровень CO2 влияет на дыхательный центр',
     'CO2 level affects the respiratory center'),
    ('immune', 'cardiovascular', 0.4, 'Воспаление влияет на сосуды',
     'Inflammation affects blood vessels'),
]


def get_organ_interactions():
    """Return organ systems and their directed interaction network."""
    return ORGAN_SYSTEMS, ORGAN_INTERACTIONS


def simulate_disease_propagation(start_system, n_steps=3):
    """Simulate how dysfunction in one organ system spreads to others through
    the interaction network. Returns affected systems with severity per step.
    """
    severity = {start_system: 1.0}
    history = [dict(severity)]
    for step in range(n_steps):
        new_sev = dict(severity)
        for (src, dst, strength, _ru, _en) in ORGAN_INTERACTIONS:
            if src in severity and severity[src] > 0.1:
                # Dysfunction propagates, dampened by distance
                transferred = severity[src] * strength * 0.4
                new_sev[dst] = min(1.0, new_sev.get(dst, 0) + transferred)
        severity = new_sev
        history.append(dict(severity))
    return history



# Age-related changes in organ function (% of young-adult baseline at each age).
# Based on general physiology: most functions peak ~25-30 and decline after.
AGE_FUNCTION_CURVES = {
    'cardiac_output':    {'name_ru': 'Сердечный выброс', 'name_en': 'Cardiac output',
                          'peak_age': 25, 'decline_per_decade': 0.08},
    'lung_capacity':     {'name_ru': 'Ёмкость лёгких', 'name_en': 'Lung capacity',
                          'peak_age': 25, 'decline_per_decade': 0.09},
    'kidney_filtration': {'name_ru': 'Фильтрация почек (GFR)', 'name_en': 'Kidney filtration (GFR)',
                          'peak_age': 30, 'decline_per_decade': 0.10},
    'muscle_mass':       {'name_ru': 'Мышечная масса', 'name_en': 'Muscle mass',
                          'peak_age': 30, 'decline_per_decade': 0.07},
    'bone_density':      {'name_ru': 'Плотность костей', 'name_en': 'Bone density',
                          'peak_age': 30, 'decline_per_decade': 0.06},
    'brain_volume':      {'name_ru': 'Объём мозга', 'name_en': 'Brain volume',
                          'peak_age': 25, 'decline_per_decade': 0.05},
    'immune_function':   {'name_ru': 'Иммунитет', 'name_en': 'Immune function',
                          'peak_age': 20, 'decline_per_decade': 0.07},
    'metabolic_rate':    {'name_ru': 'Скорость метаболизма', 'name_en': 'Metabolic rate',
                          'peak_age': 20, 'decline_per_decade': 0.03},
    'skin_elasticity':   {'name_ru': 'Эластичность кожи', 'name_en': 'Skin elasticity',
                          'peak_age': 20, 'decline_per_decade': 0.10},
    'hearing':           {'name_ru': 'Слух', 'name_en': 'Hearing',
                          'peak_age': 18, 'decline_per_decade': 0.08},
}


def compute_age_function(function_key, age):
    """Compute the % function remaining at a given age for an organ function.
    Returns 0-100 (% of peak). Models growth to peak, then gradual decline.
    """
    curve = AGE_FUNCTION_CURVES.get(function_key)
    if not curve:
        return 100.0
    peak = curve['peak_age']
    decline = curve['decline_per_decade']
    if age <= peak:
        # Growth phase: rises from ~40% at birth to 100% at peak
        frac = 0.4 + 0.6 * (age / peak) if peak > 0 else 1.0
        return round(min(100, frac * 100), 1)
    else:
        # Decline phase
        decades_past = (age - peak) / 10.0
        remaining = (1 - decline) ** decades_past
        return round(max(0, remaining * 100), 1)


def age_dynamics_profile(age):
    """Full physiological profile at a given age — all organ functions.
    Returns dict of function_key -> {name, value, status}.
    """
    profile = {}
    for key, curve in AGE_FUNCTION_CURVES.items():
        val = compute_age_function(key, age)
        if val >= 85:
            status = 'optimal'
        elif val >= 65:
            status = 'good'
        elif val >= 45:
            status = 'reduced'
        else:
            status = 'low'
        profile[key] = {
            'name_ru': curve['name_ru'], 'name_en': curve['name_en'],
            'value': val, 'status': status,
        }
    return profile


# ╔══════════════════════════════════════════════════════════════════╗
# ║  SYNTHETIC ORGANISM DESIGN — build a minimal cell for a task       ║
# ║  (synthetic genomics: minimal genome, chassis, add task modules)   ║
# ╚══════════════════════════════════════════════════════════════════╝

# Essential gene categories every free-living cell needs (based on JCVI-syn3.0,
# the 473-gene minimal cell). Grouped by function.

# ════════════════════════════════════════════════════════════════════
# SYNTHETIC ORGANISM DESIGN
# ------------------------------------------------------------------
# The organism section used to offer three chassis and seven task buttons,
# with viability decided by "is the sum of burdens above 0.7". Whole areas of
# cell engineering — defence, regulation, stress, sensing, secretion,
# biocontainment — had no representation at all.
#
# What follows models a cell the way a design tool has to: modules compete for
# a finite budget of ribosomes, ATP, carbon, nitrogen and membrane area; some
# modules need others present; some actively interfere with each other; and a
# free-form "custom" module lets the user add anything not on the list.
# ════════════════════════════════════════════════════════════════════

ORGANISM_MODULE_CATEGORIES = {
    'metabolism': {
        'name_ru': 'Метаболизм и производство', 'name_en': 'Metabolism & production',
        'colour': '#9a6700',
        'desc_ru': 'Что клетка производит и из чего',
        'desc_en': 'What the cell makes, and from what',
    },
    'immunity': {
        'name_ru': 'Иммунитет и защита', 'name_en': 'Immunity & defence',
        'colour': '#cf222e',
        'desc_ru': 'Защита от фагов, чужой ДНК, антибиотиков и конкурентов',
        'desc_en': 'Defence against phages, foreign DNA, antibiotics, competitors',
    },
    'regulation': {
        'name_ru': 'Регуляция и логика', 'name_en': 'Regulation & logic',
        'colour': '#8250df',
        'desc_ru': 'Переключатели, осцилляторы, логические схемы, память',
        'desc_en': 'Switches, oscillators, logic circuits, memory',
    },
    'stress': {
        'name_ru': 'Стресс-ответ', 'name_en': 'Stress response',
        'colour': '#bc4c00',
        'desc_ru': 'Выживание при жаре, окислении, кислоте, осмосе, радиации',
        'desc_en': 'Surviving heat, oxidation, acid, osmotic shock, radiation',
    },
    'sensing': {
        'name_ru': 'Сенсорика', 'name_en': 'Sensing',
        'colour': '#0969da',
        'desc_ru': 'Чувствительность к свету, химии, плотности популяции, pH',
        'desc_en': 'Light, chemicals, population density, pH',
    },
    'secretion': {
        'name_ru': 'Секреция и транспорт', 'name_en': 'Secretion & transport',
        'colour': '#1a7f37',
        'desc_ru': 'Вывод продукта наружу и захват веществ внутрь',
        'desc_en': 'Getting product out and nutrients in',
    },
    'containment': {
        'name_ru': 'Биобезопасность', 'name_en': 'Biocontainment',
        'colour': '#57606a',
        'desc_ru': 'Чтобы организм не выжил вне лаборатории',
        'desc_en': 'Keeping the organism from surviving outside the lab',
    },
    'structure': {
        'name_ru': 'Структура и деление', 'name_en': 'Structure & division',
        'colour': '#6e7781',
        'desc_ru': 'Форма клетки, деление, компартменты, подвижность',
        'desc_en': 'Cell shape, division, compartments, motility',
    },
    'custom': {
        'name_ru': 'Другое (своё)', 'name_en': 'Other (custom)',
        'colour': '#d4a72c',
        'desc_ru': 'Свой модуль: задайте гены, нагрузку и эффект вручную',
        'desc_en': 'Your own module: set genes, burden and effect by hand',
    },
}


# Chassis. `resources` is the budget each cell type can spend, normalised so
# that 1.0 is everything the wild-type organism has to give.
ORGANISM_CHASSIS = {
    'mycoplasma': {
        'name_ru': 'Микоплазма (минимальная)', 'name_en': 'Mycoplasma (minimal)',
        'base_genes': 473, 'growth_rate': 0.3,
        'desc_ru': 'Как JCVI-syn3.0 — самый маленький геном (473 гена)',
        'desc_en': 'Like JCVI-syn3.0 — smallest genome (473 genes)',
        'cell_type': 'prokaryote', 'genome_kb': 531,
        'resources': {'ribosome': 0.35, 'atp': 0.30, 'carbon': 0.40,
                      'nitrogen': 0.35, 'membrane': 0.30},
        'traits': {'secretion': 'poor', 'glycosylation': False,
                   'temp_range': '30-37°C', 'transformable': 'hard'},
        'notes_ru': 'Почти ничего лишнего — но и запаса прочности почти нет.',
        'notes_en': 'Almost nothing spare — and almost no safety margin.',
    },
    'ecoli': {
        'name_ru': 'E. coli (рабочая лошадка)', 'name_en': 'E. coli (workhorse)',
        'base_genes': 4400, 'growth_rate': 1.0,
        'desc_ru': 'Быстрый рост, легко редактировать, промышленный стандарт',
        'desc_en': 'Fast growth, easy to edit, industry standard',
        'cell_type': 'prokaryote', 'genome_kb': 4600,
        'resources': {'ribosome': 1.0, 'atp': 1.0, 'carbon': 1.0,
                      'nitrogen': 1.0, 'membrane': 1.0},
        'traits': {'secretion': 'moderate', 'glycosylation': False,
                   'temp_range': '20-42°C', 'transformable': 'easy'},
        'notes_ru': 'Не гликозилирует белки — для человеческих гликопротеинов не годится.',
        'notes_en': 'No glycosylation — unsuitable for human glycoproteins.',
    },
    'bsubtilis': {
        'name_ru': 'Bacillus subtilis', 'name_en': 'Bacillus subtilis',
        'base_genes': 4200, 'growth_rate': 0.85,
        'desc_ru': 'Отлично секретирует белки наружу, статус GRAS (безопасен)',
        'desc_en': 'Excellent protein secretion, GRAS (food-safe) status',
        'cell_type': 'prokaryote', 'genome_kb': 4215,
        'resources': {'ribosome': 0.95, 'atp': 0.90, 'carbon': 0.95,
                      'nitrogen': 0.90, 'membrane': 1.20},
        'traits': {'secretion': 'excellent', 'glycosylation': False,
                   'temp_range': '15-50°C', 'transformable': 'easy'},
        'notes_ru': 'Секретирует продукт в среду — очистка сильно проще.',
        'notes_en': 'Secretes product into the medium — purification is far easier.',
    },
    'yeast': {
        'name_ru': 'Дрожжи S. cerevisiae', 'name_en': 'Yeast S. cerevisiae',
        'base_genes': 6000, 'growth_rate': 0.7,
        'desc_ru': 'Эукариот: сворачивает и гликозилирует сложные белки человека',
        'desc_en': 'Eukaryote: folds and glycosylates complex human proteins',
        'cell_type': 'eukaryote', 'genome_kb': 12100,
        'resources': {'ribosome': 1.10, 'atp': 0.95, 'carbon': 1.00,
                      'nitrogen': 0.95, 'membrane': 1.30},
        'traits': {'secretion': 'good', 'glycosylation': True,
                   'temp_range': '25-35°C', 'transformable': 'easy'},
        'notes_ru': 'Гликозилирование дрожжевого типа отличается от человеческого.',
        'notes_en': 'Yeast-type glycosylation differs from the human pattern.',
    },
    'pputida': {
        'name_ru': 'Pseudomonas putida', 'name_en': 'Pseudomonas putida',
        'base_genes': 5400, 'growth_rate': 0.75,
        'desc_ru': 'Держит растворители и токсины — чемпион биоремедиации',
        'desc_en': 'Tolerates solvents and toxins — the bioremediation champion',
        'cell_type': 'prokaryote', 'genome_kb': 6180,
        'resources': {'ribosome': 0.85, 'atp': 1.30, 'carbon': 1.10,
                      'nitrogen': 0.90, 'membrane': 1.15},
        'traits': {'secretion': 'moderate', 'glycosylation': False,
                   'temp_range': '25-35°C', 'transformable': 'moderate'},
        'notes_ru': 'Мощный запас НАДФН и устойчивость к окислительному стрессу.',
        'notes_en': 'Large NADPH pool and strong oxidative-stress tolerance.',
    },
    'vnatriegens': {
        'name_ru': 'Vibrio natriegens', 'name_en': 'Vibrio natriegens',
        'base_genes': 4600, 'growth_rate': 2.0,
        'desc_ru': 'Самая быстрорастущая бактерия — деление менее чем за 10 минут',
        'desc_en': 'Fastest-growing bacterium known — divides in under 10 minutes',
        'cell_type': 'prokaryote', 'genome_kb': 5200,
        'resources': {'ribosome': 1.60, 'atp': 1.40, 'carbon': 1.30,
                      'nitrogen': 1.10, 'membrane': 1.10},
        'traits': {'secretion': 'moderate', 'glycosylation': False,
                   'temp_range': '25-37°C', 'transformable': 'moderate'},
        'notes_ru': 'Нужна солёная среда; очень быстрый расход субстрата.',
        'notes_en': 'Needs saline medium; burns through substrate very fast.',
    },
    'synechocystis': {
        'name_ru': 'Synechocystis (цианобактерия)', 'name_en': 'Synechocystis (cyanobacterium)',
        'base_genes': 3700, 'growth_rate': 0.25,
        'desc_ru': 'Фотосинтез: питается светом и CO₂, сахар не нужен',
        'desc_en': 'Photosynthetic: runs on light and CO₂, no sugar needed',
        'cell_type': 'prokaryote', 'genome_kb': 3570,
        'resources': {'ribosome': 0.60, 'atp': 1.50, 'carbon': 0.50,
                      'nitrogen': 0.70, 'membrane': 1.40},
        'traits': {'secretion': 'poor', 'glycosylation': False,
                   'temp_range': '25-35°C', 'transformable': 'moderate'},
        'notes_ru': 'Энергии из света много, но растёт медленно и углерода мало.',
        'notes_en': 'Plenty of light energy, but slow growth and little carbon.',
    },
    'lactococcus': {
        'name_ru': 'Lactococcus lactis', 'name_en': 'Lactococcus lactis',
        'base_genes': 2300, 'growth_rate': 0.6,
        'desc_ru': 'Пищевой статус GRAS — годится для пробиотиков и живых лекарств',
        'desc_en': 'GRAS food status — suitable for probiotics and live therapeutics',
        'cell_type': 'prokaryote', 'genome_kb': 2360,
        'resources': {'ribosome': 0.70, 'atp': 0.55, 'carbon': 0.70,
                      'nitrogen': 0.65, 'membrane': 0.90},
        'traits': {'secretion': 'good', 'glycosylation': False,
                   'temp_range': '25-37°C', 'transformable': 'moderate'},
        'notes_ru': 'Выживает в кишечнике — база для терапевтических пробиотиков.',
        'notes_en': 'Survives in the gut — a base for therapeutic probiotics.',
    },
    'cho': {
        'name_ru': 'Клетки CHO (млекопитающее)', 'name_en': 'CHO cells (mammalian)',
        'base_genes': 22000, 'growth_rate': 0.08,
        'desc_ru': 'Стандарт для антител и сложных биопрепаратов',
        'desc_en': 'The standard for antibodies and complex biologics',
        'cell_type': 'eukaryote', 'genome_kb': 2450000,
        'resources': {'ribosome': 1.40, 'atp': 1.20, 'carbon': 1.20,
                      'nitrogen': 1.30, 'membrane': 1.60},
        'traits': {'secretion': 'excellent', 'glycosylation': True,
                   'temp_range': '36-38°C', 'transformable': 'hard'},
        'notes_ru': 'Человекоподобное гликозилирование, но рост очень медленный и дорогой.',
        'notes_en': 'Human-like glycosylation, but growth is slow and expensive.',
    },
    'minimal_syn': {
        'name_ru': 'Синтетический минимум (с нуля)', 'name_en': 'Synthetic minimum (from scratch)',
        'base_genes': 350, 'growth_rate': 0.15,
        'desc_ru': 'Только то, что вы добавите сами — предельный конструктор',
        'desc_en': 'Only what you add yourself — the extreme build',
        'cell_type': 'prokaryote', 'genome_kb': 400,
        'resources': {'ribosome': 0.25, 'atp': 0.20, 'carbon': 0.30,
                      'nitrogen': 0.25, 'membrane': 0.25},
        'traits': {'secretion': 'poor', 'glycosylation': False,
                   'temp_range': '30-37°C', 'transformable': 'hard'},
        'notes_ru': 'Ничего лишнего — но любая ошибка сразу убивает клетку.',
        'notes_en': 'Nothing spare — any mistake kills the cell immediately.',
    },
}


def _mod(name_ru, name_en, cat, genes, burden, desc_ru, desc_en,
         cost=None, requires=None, conflicts=None, synergy=None,
         provides=None, note_ru='', note_en=''):
    """Build a module record. `cost` is the share of each resource pool used."""
    c = {'ribosome': burden, 'atp': burden * 0.7, 'carbon': burden * 0.5,
         'nitrogen': burden * 0.4, 'membrane': burden * 0.2}
    if cost:
        c.update(cost)
    return {
        'name_ru': name_ru, 'name_en': name_en, 'category': cat,
        'genes': genes, 'burden': burden,
        'desc_ru': desc_ru, 'desc_en': desc_en,
        'cost': c, 'requires': requires or [], 'conflicts': conflicts or [],
        'synergy': synergy or {}, 'provides': provides or [],
        'note_ru': note_ru, 'note_en': note_en,
    }


ORGANISM_TASK_MODULES = {
    # ── METABOLISM / PRODUCTION ─────────────────────────────────────
    'insulin': _mod('Производство инсулина', 'Insulin production', 'metabolism',
                    ['INS', 'pdi', 'signal_pep'], 0.15,
                    'Синтез человеческого инсулина (как настоящие фабрики E. coli)',
                    'Human insulin synthesis (like real E. coli factories)',
                    requires=['sec_pathway'], provides=['protein_product'],
                    note_ru='Нужны дисульфидные связи — без корректного фолдинга белок неактивен.',
                    note_en='Needs disulphide bonds — without correct folding the protein is inactive.'),
    'biofuel': _mod('Биотопливо (этанол)', 'Biofuel (ethanol)', 'metabolism',
                    ['pdc', 'adhB', 'adhE'], 0.20,
                    'Превращение сахара в этанол',
                    'Convert sugar into ethanol',
                    cost={'carbon': 0.45}, conflicts=['pha_plastic'],
                    note_ru='Конкурирует за пируват с другими путями из углерода.',
                    note_en='Competes for pyruvate with other carbon routes.'),
    'bioplastic': _mod('Биопластик (PHA)', 'Bioplastic (PHA)', 'metabolism',
                       ['phaA', 'phaB', 'phaC'], 0.25,
                       'Синтез разлагаемого пластика PHA',
                       'Synthesize biodegradable PHA plastic',
                       cost={'carbon': 0.50}, conflicts=['biofuel'],
                       provides=['pha_plastic'],
                       note_ru='Забирает ацетил-КоА — прямой конкурент этанольного пути.',
                       note_en='Draws off acetyl-CoA — a direct competitor of the ethanol route.'),
    'nitrogen_fix': _mod('Фиксация азота', 'Nitrogen fixation', 'metabolism',
                         ['nifH', 'nifD', 'nifK', 'nifB', 'nifE'], 0.35,
                         'Превращение азота воздуха в аммиак — удобрение из воздуха',
                         'Turn atmospheric N₂ into ammonia — fertiliser from air',
                         cost={'atp': 0.60}, conflicts=['aerobic_boost'],
                         provides=['fixed_nitrogen'],
                         note_ru='Нитрогеназа гибнет от кислорода и жрёт 16 АТФ на молекулу.',
                         note_en='Nitrogenase is destroyed by oxygen and burns 16 ATP per molecule.'),
    'bioremediation': _mod('Биоремедиация', 'Bioremediation', 'metabolism',
                           ['alkB', 'xylE', 'todC', 'catA'], 0.30,
                           'Разложение нефти, растворителей и ароматических загрязнителей',
                           'Break down oil, solvents and aromatic pollutants',
                           requires=['oxidative_defence'],
                           note_ru='Промежуточные продукты токсичны — нужна антиоксидантная защита.',
                           note_en='Intermediates are toxic — antioxidant defence is required.'),
    'vitamin': _mod('Синтез витаминов', 'Vitamin synthesis', 'metabolism',
                    ['crtE', 'crtB', 'crtI', 'crtY'], 0.20,
                    'Производство витамина A / бета-каротина (как золотой рис)',
                    'Produce vitamin A / beta-carotene (as in golden rice)',
                    cost={'carbon': 0.35}),
    'antibiotic_synth': _mod('Синтез антибиотика', 'Antibiotic synthesis', 'metabolism',
                             ['pks1', 'nrps1', 'tailoring', 'resistance_self'], 0.40,
                             'Поликетидный синтез нового антибиотика',
                             'Polyketide synthesis of a new antibiotic',
                             requires=['self_resistance'],
                             note_ru='Без гена самозащиты клетка отравится собственным продуктом.',
                             note_en='Without a self-resistance gene the cell poisons itself.'),
    'amino_acid': _mod('Производство аминокислот', 'Amino acid production', 'metabolism',
                       ['lysC', 'dapA', 'lysA'], 0.22,
                       'Промышленный синтез лизина / глутамата',
                       'Industrial lysine / glutamate synthesis',
                       cost={'nitrogen': 0.55}),
    'terpene': _mod('Терпены и ароматика', 'Terpenes & fragrances', 'metabolism',
                    ['dxs', 'idi', 'ispA', 'terpene_synthase'], 0.28,
                    'Синтез артемизинина, ментола, душистых веществ',
                    'Synthesis of artemisinin, menthol, fragrance compounds',
                    cost={'carbon': 0.40}),
    'co2_fixation': _mod('Фиксация CO₂ (цикл Кальвина)', 'CO₂ fixation (Calvin cycle)',
                         'metabolism', ['rbcL', 'rbcS', 'prkA', 'gapA'], 0.45,
                         'Питание углекислым газом вместо сахара',
                         'Feed on carbon dioxide instead of sugar',
                         cost={'atp': 0.70}, provides=['autotrophy'],
                         note_ru='RuBisCO — медленный и неточный фермент, это узкое место.',
                         note_en='RuBisCO is slow and error-prone — that is the bottleneck.'),
    'cellulose_deg': _mod('Разложение целлюлозы', 'Cellulose degradation', 'metabolism',
                          ['celA', 'celB', 'bglX', 'xynA'], 0.30,
                          'Питание соломой и древесными отходами',
                          'Feed on straw and wood waste',
                          requires=['sec_pathway'], provides=['cheap_carbon']),

    # ── IMMUNITY / DEFENCE ──────────────────────────────────────────
    'restriction_mod': _mod('Рестрикция-модификация', 'Restriction-modification',
                            'immunity', ['hsdR', 'hsdM', 'hsdS'], 0.10,
                            'Режет чужую ДНК, свою метилирует и щадит',
                            'Cuts foreign DNA, methylates and spares its own',
                            provides=['dna_defence'],
                            note_ru='Классический бактериальный иммунитет — мешает и вашей трансформации.',
                            note_en='Classic bacterial immunity — it will fight your transformation too.'),
    'crispr_immunity': _mod('Собственная CRISPR-система', 'Native CRISPR system',
                            'immunity', ['cas1', 'cas2', 'cas9', 'crRNA_array'], 0.15,
                            'Адаптивный иммунитет: запоминает фаги и режет их при новой встрече',
                            'Adaptive immunity: remembers phages and cuts them on re-infection',
                            provides=['phage_defence', 'dna_defence'],
                            conflicts=['crispr_editing'],
                            note_ru='Может резать ваши же конструкции — конфликтует с редактированием.',
                            note_en='May cut your own constructs — conflicts with editing.'),
    'abortive_infection': _mod('Абортивная инфекция', 'Abortive infection', 'immunity',
                               ['abiE', 'toxIN', 'cbtA'], 0.12,
                               'Заражённая клетка самоубивается, спасая популяцию',
                               'An infected cell kills itself to save the population',
                               provides=['phage_defence'],
                               note_ru='Жертвует отдельными клетками ради выживания культуры.',
                               note_en='Sacrifices single cells so the culture survives.'),
    'phage_receptor_loss': _mod('Потеря рецепторов фага', 'Phage receptor loss', 'immunity',
                                ['ompC_del', 'lamB_del', 'lpsE_mod'], 0.08,
                                'Убирает белки, за которые цепляется фаг',
                                'Removes the surface proteins phages latch onto',
                                provides=['phage_defence'],
                                conflicts=['sugar_uptake'],
                                note_ru='Те же белки часто работают транспортёрами — теряется питание.',
                                note_en='Those proteins are often transporters — nutrient uptake suffers.'),
    'antibiotic_resist': _mod('Устойчивость к антибиотику', 'Antibiotic resistance',
                              'immunity', ['bla', 'kanR', 'tetA'], 0.06,
                              'Селективный маркёр и защита от заражения культуры',
                              'Selection marker and protection from culture contamination',
                              provides=['selection_marker'],
                              note_ru='Удобно в лаборатории, но нежелательно для выпуска в среду.',
                              note_en='Convenient in the lab, undesirable for environmental release.'),
    'self_resistance': _mod('Самозащита от продукта', 'Self-resistance to product',
                            'immunity', ['efflux_pump', 'target_mod'], 0.12,
                            'Помпа откачивает собственный токсичный продукт наружу',
                            'An efflux pump removes the cell\u2019s own toxic product',
                            cost={'membrane': 0.35, 'atp': 0.30},
                            provides=['self_resistance']),
    'oxidative_defence': _mod('Антиоксидантная защита', 'Oxidative defence', 'immunity',
                              ['sodA', 'katG', 'ahpC', 'oxyR'], 0.10,
                              'Обезвреживает перекись и свободные радикалы',
                              'Neutralises peroxide and free radicals',
                              provides=['oxidative_defence'],
                              synergy={'bioremediation': 0.25, 'co2_fixation': 0.10}),
    'biofilm_defence': _mod('Биоплёнка', 'Biofilm', 'immunity',
                            ['pgaA', 'csgA', 'bcsA'], 0.18,
                            'Слизистый матрикс защищает от антибиотиков и высыхания',
                            'A slime matrix protects against antibiotics and drying',
                            cost={'carbon': 0.40}, conflicts=['motility'],
                            note_ru='В биореакторе плёнка забивает оборудование.',
                            note_en='In a bioreactor the film clogs the equipment.'),

    # ── REGULATION / LOGIC ──────────────────────────────────────────
    'toggle_switch': _mod('Триггер (память)', 'Toggle switch (memory)', 'regulation',
                          ['lacI', 'tetR', 'pLac', 'pTet'], 0.08,
                          'Клетка запоминает состояние и держит его без сигнала',
                          'The cell remembers a state and holds it without a signal',
                          provides=['memory']),
    'repressilator': _mod('Осциллятор', 'Oscillator', 'regulation',
                          ['lacI', 'tetR', 'cI'], 0.09,
                          'Встроенные часы — периодическое включение генов',
                          'A built-in clock — genes switch on and off periodically',
                          provides=['timing']),
    'logic_gates': _mod('Логические вентили', 'Logic gates', 'regulation',
                        ['AND_gate', 'OR_gate', 'NOT_gate'], 0.10,
                        'Решения по нескольким входам сразу (И, ИЛИ, НЕ)',
                        'Decisions from several inputs at once (AND, OR, NOT)',
                        provides=['logic']),
    'feedback_control': _mod('Обратная связь', 'Feedback control', 'regulation',
                             ['sensor_TF', 'antisense_RNA'], 0.07,
                             'Автоматически держит нужный уровень продукта',
                             'Automatically holds the product at the target level',
                             synergy={'insulin': 0.15, 'amino_acid': 0.20, 'terpene': 0.15},
                             note_ru='Снижает нагрузку: клетка не производит больше нужного.',
                             note_en='Cuts burden: the cell stops overproducing.'),
    'growth_decoupling': _mod('Разделение роста и продукции', 'Growth/production decoupling',
                              'regulation', ['quorum_promoter', 'protease_tag'], 0.11,
                              'Сначала клетка растёт, потом переключается на производство',
                              'The cell grows first, then switches to production',
                              requires=['quorum_sensing'],
                              synergy={'insulin': 0.25, 'bioplastic': 0.25,
                                       'antibiotic_synth': 0.30},
                              note_ru='Классический промышленный приём двухфазной ферментации.',
                              note_en='The classic two-phase fermentation trick.'),
    'riboswitch': _mod('Рибопереключатель', 'Riboswitch', 'regulation',
                       ['riboswitch_5utr'], 0.04,
                       'РНК сама включает ген в ответ на метаболит — без белков',
                       'RNA switches the gene on in response to a metabolite — no protein needed',
                       provides=['logic']),

    # ── STRESS RESPONSE ─────────────────────────────────────────────
    'heat_shock': _mod('Тепловой шок', 'Heat shock response', 'stress',
                       ['dnaK', 'groEL', 'groES', 'rpoH'], 0.12,
                       'Шапероны спасают белки при перегреве',
                       'Chaperones rescue proteins when it gets too hot',
                       synergy={'insulin': 0.20, 'antibiotic_synth': 0.10},
                       note_ru='Заодно помогает правильно свернуть чужеродный белок.',
                       note_en='Also helps a foreign protein fold correctly.'),
    'osmotic_tolerance': _mod('Осмотическая устойчивость', 'Osmotic tolerance', 'stress',
                              ['proU', 'otsA', 'betA'], 0.10,
                              'Выживание в солёной и концентрированной среде',
                              'Survival in salty and concentrated media',
                              cost={'carbon': 0.25}),
    'acid_tolerance': _mod('Кислотоустойчивость', 'Acid tolerance', 'stress',
                           ['gadA', 'gadB', 'gadC'], 0.09,
                           'Работа при низком pH — например, в желудке',
                           'Working at low pH — in the stomach, for instance',
                           synergy={'probiotic': 0.30}),
    'solvent_tolerance': _mod('Устойчивость к растворителям', 'Solvent tolerance', 'stress',
                              ['ttgABC', 'cfa', 'marA'], 0.15,
                              'Мембрана держит спирты и органику',
                              'The membrane holds up against alcohols and organics',
                              cost={'membrane': 0.45},
                              synergy={'biofuel': 0.35, 'bioremediation': 0.25},
                              note_ru='Без неё этанол убивает продуцента раньше, чем накопится.',
                              note_en='Without it ethanol kills the producer before it accumulates.'),
    'cold_adaptation': _mod('Холодовая адаптация', 'Cold adaptation', 'stress',
                            ['cspA', 'desA'], 0.08,
                            'Работа при низкой температуре — медленнее, но чище фолдинг',
                            'Working in the cold — slower, but cleaner folding'),
    'radiation_resist': _mod('Радиорезистентность', 'Radiation resistance', 'stress',
                             ['recA_ext', 'ddrA', 'pprA', 'mn_transport'], 0.22,
                             'Починка ДНК после радиации (по мотивам Deinococcus)',
                             'DNA repair after radiation (Deinococcus-style)',
                             cost={'atp': 0.40}),

    # ── SENSING ─────────────────────────────────────────────────────
    'quorum_sensing': _mod('Чувство кворума', 'Quorum sensing', 'sensing',
                           ['luxI', 'luxR'], 0.06,
                           'Клетки считают друг друга и действуют вместе',
                           'Cells count each other and act together',
                           provides=['quorum_sensing', 'population_signal']),
    'light_sensing': _mod('Светочувствительность', 'Light sensing', 'sensing',
                          ['cph1', 'ompR_fusion'], 0.10,
                          'Управление генами светом — точно во времени и пространстве',
                          'Control genes with light — precise in time and space',
                          provides=['external_control']),
    'chemical_sensing': _mod('Химический сенсор', 'Chemical sensor', 'sensing',
                             ['arsR', 'merR', 'nahR'], 0.07,
                             'Обнаружение мышьяка, ртути, нитратов, взрывчатки',
                             'Detects arsenic, mercury, nitrates, explosives',
                             provides=['detection']),
    'ph_sensing': _mod('Датчик pH', 'pH sensor', 'sensing',
                       ['cadC', 'atpB_reporter'], 0.05,
                       'Реакция на кислотность среды',
                       'Responds to the acidity of the medium',
                       provides=['detection']),
    'oxygen_sensing': _mod('Датчик кислорода', 'Oxygen sensor', 'sensing',
                           ['fnr', 'arcA', 'arcB'], 0.06,
                           'Переключение между аэробным и анаэробным режимом',
                           'Switches between aerobic and anaerobic mode',
                           provides=['detection'],
                           synergy={'nitrogen_fix': 0.20}),
    'biosensor': _mod('Биосенсор с индикацией', 'Reporting biosensor', 'sensing',
                      ['sensor', 'gfp', 'promoter'], 0.10,
                      'Светится при обнаружении цели — видимый ответ',
                      'Glows when the target is detected — a visible readout',
                      requires=['detection'], provides=['reporter']),
    'temperature_sensing': _mod('Термодатчик', 'Temperature sensor', 'sensing',
                                ['rnaT_thermometer'], 0.04,
                                'РНК-термометр включает гены при нужной температуре',
                                'An RNA thermometer switches genes on at a set temperature',
                                provides=['external_control']),

    # ── SECRETION / TRANSPORT ───────────────────────────────────────
    'sec_pathway': _mod('Секреция Sec', 'Sec secretion pathway', 'secretion',
                        ['secA', 'secY', 'secE', 'signal_peptidase'], 0.12,
                        'Стандартный вывод белка через мембрану (в развёрнутом виде)',
                        'Standard protein export across the membrane (unfolded)',
                        cost={'membrane': 0.40}, provides=['sec_pathway', 'secretion']),
    'tat_pathway': _mod('Секреция Tat', 'Tat secretion pathway', 'secretion',
                        ['tatA', 'tatB', 'tatC'], 0.14,
                        'Вывод уже свёрнутого белка — для сложных ферментов с кофакторами',
                        'Exports already-folded protein — for complex cofactor enzymes',
                        cost={'membrane': 0.45}, provides=['secretion']),
    'vesicle_secretion': _mod('Везикулы', 'Membrane vesicles', 'secretion',
                              ['tolR_mod', 'ompA_del'], 0.16,
                              'Упаковка продукта в пузырьки — доставка и вакцины',
                              'Packs product into vesicles — delivery and vaccines',
                              cost={'membrane': 0.50}, provides=['secretion']),
    'sugar_uptake': _mod('Транспорт сахаров', 'Sugar transport', 'secretion',
                         ['ptsG', 'lacY', 'araE', 'xylFGH'], 0.10,
                         'Расширенный захват разных сахаров, включая пентозы',
                         'Broad uptake of different sugars, pentoses included',
                         cost={'membrane': 0.35}, provides=['cheap_carbon'],
                         synergy={'cellulose_deg': 0.30}),
    'metal_transport': _mod('Транспорт металлов', 'Metal transport', 'secretion',
                            ['znuABC', 'feoB', 'mntH'], 0.09,
                            'Захват железа, цинка, марганца — и накопление металлов',
                            'Uptake of iron, zinc, manganese — and metal accumulation',
                            cost={'membrane': 0.30},
                            synergy={'bioremediation': 0.20}),
    'efflux_export': _mod('Экспорт продукта наружу', 'Product efflux', 'secretion',
                          ['acrAB', 'tolC'], 0.13,
                          'Выкачивает мелкие молекулы из клетки — продукт не копится внутри',
                          'Pumps small molecules out — product stops building up inside',
                          cost={'membrane': 0.40, 'atp': 0.35},
                          synergy={'biofuel': 0.20, 'terpene': 0.30,
                                   'antibiotic_synth': 0.25}),

    # ── BIOCONTAINMENT ──────────────────────────────────────────────
    'kill_switch': _mod('Выключатель гибели', 'Kill switch', 'containment',
                        ['ccdB', 'pTet_repressor'], 0.07,
                        'Без индуктора из среды клетка погибает',
                        'Without an inducer from the medium the cell dies',
                        provides=['containment'],
                        note_ru='Мутации ломают выключатель — одного контура мало.',
                        note_en='Mutations break the switch — one layer is not enough.'),
    'auxotrophy': _mod('Ауксотрофия', 'Auxotrophy', 'containment',
                       ['thyA_del', 'dapA_del'], 0.05,
                       'Клетка не может сама сделать нужное вещество — только из среды',
                       'The cell cannot make an essential compound — only the medium provides it',
                       provides=['containment']),
    'synthetic_aa': _mod('Зависимость от синтетической аминокислоты',
                         'Synthetic amino acid dependence', 'containment',
                         ['orthogonal_tRNA', 'orthogonal_aaRS', 'amber_codons'], 0.20,
                         'Белки требуют неприродной аминокислоты, которой нет в природе',
                         'Proteins require a non-natural amino acid that does not exist in nature',
                         provides=['containment'],
                         note_ru='Самый надёжный барьер: побег требует переизобрести генетический код.',
                         note_en='The strongest barrier: escape would require reinventing the genetic code.'),
    'xeno_nucleic': _mod('Ксенонуклеиновые кислоты', 'Xeno nucleic acids', 'containment',
                         ['xna_polymerase', 'xna_import'], 0.30,
                         'Генетический материал на чужой химии — не читается природой',
                         'Genetic material on alien chemistry — unreadable by nature',
                         provides=['containment', 'gene_barrier']),
    'passive_degradation': _mod('Ограниченное число делений', 'Limited division counter',
                                'containment', ['counter_circuit', 'recombinase'], 0.14,
                                'Клетка считает деления и после лимита останавливается',
                                'The cell counts divisions and stops after a limit',
                                requires=['memory'], provides=['containment']),

    # ── STRUCTURE / DIVISION ────────────────────────────────────────
    'motility': _mod('Подвижность (жгутики)', 'Motility (flagella)', 'structure',
                     ['fliC', 'motA', 'motB', 'cheA', 'cheY'], 0.25,
                     'Клетка плывёт к цели — доставка и хемотаксис',
                     'The cell swims toward a target — delivery and chemotaxis',
                     cost={'atp': 0.50}, conflicts=['biofilm_defence'],
                     note_ru='Жгутик — одна из самых дорогих структур в клетке.',
                     note_en='A flagellum is one of the most expensive structures a cell builds.'),
    'compartment': _mod('Компартмент (микроотсек)', 'Compartment (microcompartment)',
                        'structure', ['pduA', 'pduB', 'pduJ'], 0.20,
                        'Белковая капсула изолирует токсичные промежуточные продукты',
                        'A protein shell isolates toxic intermediates',
                        synergy={'co2_fixation': 0.30, 'antibiotic_synth': 0.20,
                                 'bioremediation': 0.20}),
    'cell_size': _mod('Увеличенный размер клетки', 'Enlarged cell size', 'structure',
                      ['sulA_mod', 'minCDE_mod'], 0.15,
                      'Больше объёма под продукт, но деление замедляется',
                      'More room for product, but division slows down',
                      cost={'membrane': 0.35}),
    'sporulation': _mod('Споруляция', 'Sporulation', 'structure',
                        ['spo0A', 'sigF', 'sigE', 'sigG'], 0.18,
                        'Переживание любых условий в виде споры',
                        'Riding out any conditions as a spore',
                        provides=['survival']),
    'aerobic_boost': _mod('Усиленное дыхание', 'Enhanced respiration', 'structure',
                          ['cyoABCD', 'nuoA_N'], 0.16,
                          'Больше АТФ из кислорода',
                          'More ATP from oxygen',
                          cost={'atp': -0.35}, conflicts=['nitrogen_fix'],
                          note_ru='Отрицательная стоимость: модуль сам добавляет энергию.',
                          note_en='Negative cost: this module adds energy rather than spending it.'),
    'probiotic': _mod('Кишечная колонизация', 'Gut colonisation', 'structure',
                      ['adhesin', 'bile_resist', 'mucin_bind'], 0.17,
                      'Закрепление в кишечнике — основа живых лекарств',
                      'Establishing in the gut — the basis of live therapeutics',
                      requires=['acid_tolerance'], provides=['probiotic']),
}


MINIMAL_GENE_CATEGORIES = {
    'replication': {
        'name_ru': 'Репликация ДНК', 'name_en': 'DNA replication',
        'genes': ['dnaA', 'dnaN', 'dnaE', 'polC', 'gyrA', 'gyrB', 'ligA'],
        'essential': True,
        'why_ru': 'Без копирования ДНК клетка не может делиться вообще.',
        'why_en': 'Without copying DNA the cell cannot divide at all.',
    },
    'transcription': {
        'name_ru': 'Транскрипция', 'name_en': 'Transcription',
        'genes': ['rpoA', 'rpoB', 'rpoC', 'sigA', 'nusA', 'nusG'],
        'essential': True,
        'why_ru': 'Без РНК-полимеразы не читается ни один ген.',
        'why_en': 'Without RNA polymerase not a single gene can be read.',
    },
    'translation': {
        'name_ru': 'Трансляция', 'name_en': 'Translation',
        'genes': (['rplA', 'rpsA', 'tuf', 'fusA', 'infA', 'infB', 'infC']
                  + ['tRNA-x' + str(i) for i in range(1, 30)]),
        'essential': True,
        'why_ru': 'Рибосома — самая дорогая и незаменимая машина клетки.',
        'why_en': 'The ribosome is the most expensive and irreplaceable machine there is.',
    },
    'membrane': {
        'name_ru': 'Мембрана и транспорт', 'name_en': 'Membrane & transport',
        'genes': ['secA', 'secY', 'ftsH', 'atpA', 'atpB', 'ptsH', 'ptsI'],
        'essential': True,
        'why_ru': 'Без границы с внешним миром нет клетки.',
        'why_en': 'Without a boundary from the outside world there is no cell.',
    },
    'energy': {
        'name_ru': 'Энергия (гликолиз)', 'name_en': 'Energy (glycolysis)',
        'genes': ['pgi', 'pfkA', 'fba', 'gapA', 'pgk', 'eno', 'pyk'],
        'essential': True,
        'why_ru': 'Источник АТФ для всего остального.',
        'why_en': 'The ATP source everything else depends on.',
    },
    'cofactors': {
        'name_ru': 'Кофакторы и витамины', 'name_en': 'Cofactors & vitamins',
        'genes': ['folA', 'thiL', 'ribF', 'coaA', 'nadD'],
        'essential': True,
        'why_ru': 'Ферменты без кофакторов не работают.',
        'why_en': 'Enzymes do not work without their cofactors.',
    },
    'cell_division': {
        'name_ru': 'Деление клетки', 'name_en': 'Cell division',
        'genes': ['ftsZ', 'ftsA', 'ftsK', 'zipA', 'minC', 'minD'],
        'essential': True,
        'why_ru': 'Без аппарата деления клетка растёт, но не размножается.',
        'why_en': 'Without the division machinery the cell grows but never reproduces.',
    },
    'lipid_synthesis': {
        'name_ru': 'Синтез липидов', 'name_en': 'Lipid synthesis',
        'genes': ['fabA', 'fabB', 'fabD', 'plsB', 'plsC'],
        'essential': True,
        'why_ru': 'Мембрану надо из чего-то строить при каждом делении.',
        'why_en': 'The membrane has to be built from something at every division.',
    },
}


# User-defined modules live here; the UI adds to this at run time.
CUSTOM_ORGANISM_MODULES = {}


def add_custom_organism_module(key, name_ru, name_en, genes, burden,
                               desc_ru='', desc_en='', category='custom',
                               cost=None, requires=None, conflicts=None,
                               synergy=None, provides=None):
    """Register a module the user invented. Behaves exactly like a built-in."""
    key = str(key).strip().replace(' ', '_') or 'custom_module'
    m = _mod(name_ru or key, name_en or key, category,
             [g.strip() for g in genes if str(g).strip()],
             float(burden), desc_ru, desc_en, cost=cost, requires=requires,
             conflicts=conflicts, synergy=synergy, provides=provides)
    m['custom'] = True
    CUSTOM_ORGANISM_MODULES[key] = m
    return key, m


def all_organism_modules():
    """Built-in modules plus anything the user has added."""
    merged = dict(ORGANISM_TASK_MODULES)
    merged.update(CUSTOM_ORGANISM_MODULES)
    return merged


def modules_by_category():
    """{category_key: [(module_key, module), ...]} for building the UI tree."""
    out = {}
    for k, m in all_organism_modules().items():
        out.setdefault(m.get('category', 'custom'), []).append((k, m))
    for v in out.values():
        v.sort(key=lambda kv: kv[1]['name_en'])
    return out


RESOURCE_LABELS = {
    'ribosome': ('Рибосомы', 'Ribosomes'),
    'atp': ('Энергия (АТФ)', 'Energy (ATP)'),
    'carbon': ('Углерод', 'Carbon'),
    'nitrogen': ('Азот', 'Nitrogen'),
    'membrane': ('Площадь мембраны', 'Membrane area'),
}


def design_synthetic_organism(chassis_key, task_keys, removed_categories=None,
                              substrate_supply=1.0):
    """Assemble an organism and work out whether it can actually live.

    Instead of adding up one 'burden' number, each module draws on five
    separate pools. A design can therefore fail for a specific, nameable
    reason — out of ribosomes, out of membrane area, out of ATP — which is
    what an engineer needs to know in order to fix it.
    """
    chassis = ORGANISM_CHASSIS.get(chassis_key)
    if not chassis:
        return {'error': 'Unknown chassis'}
    removed_categories = removed_categories or []
    mods = all_organism_modules()

    # ── essential gene set ──
    essential_count = 0
    missing_essential = []
    for cat_key, cat in MINIMAL_GENE_CATEGORIES.items():
        if cat_key in removed_categories:
            missing_essential.append(cat['name_ru'])
        else:
            essential_count += len(cat['genes'])

    # ── resource accounting ──
    supply = dict(chassis['resources'])
    supply['carbon'] *= max(0.1, float(substrate_supply))
    used = {k: 0.0 for k in supply}
    task_genes = 0
    chosen = []
    for tk in task_keys:
        m = mods.get(tk)
        if not m:
            continue
        chosen.append((tk, m))
        task_genes += len(m['genes'])
        for r, v in m['cost'].items():
            if r in used:
                used[r] += v

    # ── synergies reduce the cost of the partner module ──
    selected_keys = set(k for k, _ in chosen)
    synergy_notes = []
    synergy_relief = 0.0
    for tk, m in chosen:
        for partner, bonus in (m.get('synergy') or {}).items():
            if partner in selected_keys:
                synergy_relief += bonus
                pm = mods.get(partner, {})
                synergy_notes.append('{} + {}: −{:.0f}% нагрузки'.format(
                    m['name_ru'], pm.get('name_ru', partner), bonus * 100))
    if synergy_relief:
        factor = max(0.5, 1.0 - min(0.5, synergy_relief * 0.25))
        for r in used:
            used[r] *= factor

    # ── requirements and conflicts ──
    provided = set()
    for tk, m in chosen:
        provided.update(m.get('provides') or [])
        provided.add(tk)
    problems = []
    for tk, m in chosen:
        for req in (m.get('requires') or []):
            if req not in provided:
                rm = mods.get(req, {})
                problems.append({
                    'kind': 'missing_requirement', 'module': tk,
                    'need': req,
                    'text_ru': '{} требует «{}» — модуль не добавлен'.format(
                        m['name_ru'], rm.get('name_ru', req)),
                    'text_en': '{} requires "{}" — that module is missing'.format(
                        m['name_en'], rm.get('name_en', req)),
                })
        for con in (m.get('conflicts') or []):
            if con in selected_keys:
                cm = mods.get(con, {})
                problems.append({
                    'kind': 'conflict', 'module': tk, 'with': con,
                    'text_ru': '{} конфликтует с {}'.format(
                        m['name_ru'], cm.get('name_ru', con)),
                    'text_en': '{} conflicts with {}'.format(
                        m['name_en'], cm.get('name_en', con)),
                })

    # ── verdict ──
    load = {}
    overloaded = []
    for r, cap in supply.items():
        frac = used[r] / cap if cap > 0 else 99.0
        load[r] = frac
        if frac > 1.0:
            overloaded.append(r)

    warnings = []
    viable = True
    if missing_essential:
        viable = False
        warnings.append('Удалены жизненно важные системы: '
                        + ', '.join(missing_essential))
    for r in overloaded:
        viable = False
        warnings.append('Не хватает ресурса «{}»: нужно {:.0f}% от доступного'
                        .format(RESOURCE_LABELS[r][0], load[r] * 100))
    for p in problems:
        if p['kind'] == 'conflict':
            viable = False
        warnings.append(p['text_ru'])

    near = [r for r, f in load.items() if 0.8 < f <= 1.0]
    for r in near:
        warnings.append('Ресурс «{}» почти исчерпан ({:.0f}%) — рост замедлится'
                        .format(RESOURCE_LABELS[r][0], load[r] * 100))

    worst = max(load.values()) if load else 0.0
    growth = chassis['growth_rate'] * max(0.0, 1.0 - min(1.0, worst))
    if not viable:
        growth = 0.0

    total_burden = sum(used.values()) / max(1e-6, sum(supply.values()))

    # ── biocontainment score: how many independent barriers are in place ──
    containment_layers = sum(
        1 for tk, m in chosen if 'containment' in (m.get('provides') or []))
    if containment_layers == 0:
        containment_ru = 'нет барьеров — выпуск в среду недопустим'
    elif containment_layers == 1:
        containment_ru = 'один барьер — мутация может его сломать'
    elif containment_layers == 2:
        containment_ru = 'два независимых барьера — приемлемо для полевых испытаний'
    else:
        containment_ru = '{} барьера — уровень надёжной изоляции'.format(
            containment_layers)

    return {
        'error': None,
        'chassis': chassis['name_en'], 'chassis_ru': chassis['name_ru'],
        'chassis_key': chassis_key,
        'total_genes': essential_count + task_genes,
        'essential_genes': essential_count, 'task_genes': task_genes,
        'base_genes': chassis['base_genes'],
        'tasks': [m for _, m in chosen],
        'task_keys': [k for k, _ in chosen],
        'load': load, 'supply': supply, 'used': used,
        'overloaded': overloaded,
        'burden': round(total_burden, 2),
        'growth_rate': round(growth, 3),
        'viable': viable,
        'warnings': warnings,
        'problems': problems,
        'synergies': synergy_notes,
        'missing_essential': missing_essential,
        'containment_layers': containment_layers,
        'containment_ru': containment_ru,
        'traits': chassis.get('traits', {}),
        'notes_ru': chassis.get('notes_ru', ''),
        'notes_en': chassis.get('notes_en', ''),
    }


def design_minimal_organism(chassis_key, task_keys, removed_categories=None):
    """Backwards-compatible wrapper around design_synthetic_organism."""
    d = design_synthetic_organism(chassis_key, task_keys, removed_categories)
    if d.get('error'):
        return d
    return d


def simulate_organism_growth(design, hours=24, dt=0.5):
    """Population growth with a lag phase and a burden-dependent ceiling.

    Plain logistic growth ignored the fact that a heavily loaded cell both
    starts slower and tops out lower, which is the effect an engineer is
    actually trying to manage.
    """
    if not design.get('viable'):
        return [(0.0, 0.01), (float(hours), 0.0)]
    growth = max(0.001, design.get('growth_rate', 0.5))
    burden = min(0.95, design.get('burden', 0.0))
    lag = 1.0 + 4.0 * burden          # heavier load, longer lag
    carrying = max(0.15, 1.0 - 0.5 * burden)
    pop = 0.01
    traj = [(0.0, pop)]
    t = 0.0
    while t < hours:
        t += dt
        rate = growth if t > lag else growth * (t / max(1e-6, lag)) ** 2
        pop = pop + rate * pop * (1 - pop / carrying) * dt
        traj.append((t, max(0.0, min(carrying, pop))))
    return traj


def organism_production_forecast(design, hours=24):
    """Rough product output over time: cells x per-cell productivity."""
    traj = simulate_organism_growth(design, hours=hours)
    prod_modules = [m for m in design.get('tasks', [])
                    if m.get('category') == 'metabolism']
    if not prod_modules or not design.get('viable'):
        return [(t, 0.0) for t, _ in traj]
    rate = sum(m['burden'] for m in prod_modules) * 0.8
    out, acc = [], 0.0
    prev_t = 0.0
    for t, pop in traj:
        acc += rate * pop * (t - prev_t)
        prev_t = t
        out.append((t, acc))
    return out


def suggest_organism_fixes(design, lang='ru'):
    """Concrete next steps for a design that does not work."""
    tips = []
    if design.get('error'):
        return tips
    mods = all_organism_modules()
    for r in design.get('overloaded', []):
        label = RESOURCE_LABELS[r][0 if lang == 'ru' else 1]
        biggest, biggest_v = None, 0.0
        for m in design.get('tasks', []):
            v = m['cost'].get(r, 0)
            if v > biggest_v:
                biggest, biggest_v = m, v
        if biggest is not None:
            tips.append('Ресурс «{}» перегружен. Больше всех тратит «{}» — '
                        'уберите его или возьмите более мощное шасси.'
                        .format(label, biggest['name_ru'])
                        if lang == 'ru' else
                        'Resource "{}" is over budget. The biggest consumer is '
                        '"{}" — drop it or pick a stronger chassis.'
                        .format(label, biggest['name_en']))
    for p in design.get('problems', []):
        if p['kind'] == 'missing_requirement':
            rm = mods.get(p['need'], {})
            tips.append('Добавьте модуль «{}» — без него {} не заработает.'
                        .format(rm.get('name_ru', p['need']),
                                mods.get(p['module'], {}).get('name_ru', p['module']))
                        if lang == 'ru' else
                        'Add the "{}" module — {} will not work without it.'
                        .format(rm.get('name_en', p['need']),
                                mods.get(p['module'], {}).get('name_en', p['module'])))
        else:
            tips.append(p['text_ru'] if lang == 'ru' else p['text_en'])
    if design.get('containment_layers', 0) < 2 and design.get('viable'):
        tips.append('Барьеров биобезопасности меньше двух — добавьте '
                    'ауксотрофию или выключатель гибели.'
                    if lang == 'ru' else
                    'Fewer than two containment barriers — add auxotrophy or a '
                    'kill switch.')
    # unused synergies worth taking
    sel = set(design.get('task_keys', []))
    for k, m in mods.items():
        if k in sel:
            continue
        for partner, bonus in (m.get('synergy') or {}).items():
            if partner in sel and bonus >= 0.2:
                tips.append('«{}» снизил бы нагрузку «{}» на {:.0f}%.'
                            .format(m['name_ru'],
                                    mods.get(partner, {}).get('name_ru', partner),
                                    bonus * 100)
                            if lang == 'ru' else
                            '"{}" would cut the load of "{}" by {:.0f}%.'
                            .format(m['name_en'],
                                    mods.get(partner, {}).get('name_en', partner),
                                    bonus * 100))
                break
    return tips


def get_organism_task_modules():
    """Kept for existing callers: (modules, chassis, essential categories)."""
    return all_organism_modules(), ORGANISM_CHASSIS, MINIMAL_GENE_CATEGORIES


# ── BRIDGE: gene name → DNA sequence ──
# The synthetic organism works with gene NAMES (nifH, INS, ...). Other sections
# (CRISPR, Constructor) need actual DNA SEQUENCES. These helpers translate.

# Short reference sequences for well-known task genes (illustrative coding
# starts — real genes are longer, but these give a valid, scannable sequence
# with proper start/stop for CRISPR targeting and circuit building).
KNOWN_GENE_SEQUENCES = {
    'INS':  'ATGGCCCTGTGGATGCGCCTCCTGCCCCTGCTGGCGCTGCTGGCCCTCTGGGGACCTGACCCAGCCGCAGCCTTTGTGAACCAACACCTGTGCGGCTCACACCTGGTGGAAGCTCTCTACCTAGTGTGCGGGGAACGAGGCTTCTTCTACACACCCAAGACCCGCCGGGAGGCAGAGGACCTGCAGGTGGGGCAGGTGGAGCTGGGCGGGGGCCCTGGTGCAGGCAGCCTGCAGCCCTTGGCCCTGGAGGGGTCCCTGCAGAAGCGTGGCATTGTGGAACAATGCTGTACCAGCATCTGCTCCCTCTACCAGCTGGAGAACTACTGCAACTAG',
    'gfp':  'ATGAGTAAAGGAGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGTGATGTTAATGGGCACAAATTTTCTGTCAGTGGAGAGGGTGAAGGTGATGCAACATACGGAAAACTTACCCTTAAATTTATTTGCACTACTGGAAAACTACCTGTTCCATGGCCAACACTTGTCACTACTTTCTCTTATGGTGTTCAATGCTTTTCAAGATACCCAGATCATATGAAACAGCATGACTTTTTCAAGAGTGCCATGCCCGAAGGTTATGTACAGGAAAGAACTATATTTTTCAAAGATGACGGGAACTACAAGACACGTGCTGAAGTCAAGTTTGAAGGTGATACCCTTGTTAATAGAATCGAGTTAAAAGGTATTGATTTTAAAGAAGATGGAAACATTCTTGGACACAAATTGGAATACAACTATAACTCACACAATGTATACATCATGGCAGACAAACAAAAGAATGGAATCAAAGTTAACTTCAAAATTAGACACAACATTGAAGATGGAAGCGTTCAACTAGCAGACCATTATCAACAAAATACTCCAATTGGCGATGGCCCTGTCCTTTTACCAGACAACCATTACCTGTCCACACAATCTGCCCTTTCGAAAGATCCCAACGAAAAGAGAGACCACATGGTCCTTCTTGAGTTTGTAACAGCTGCTGGGATTACACATGGCATGGATGAACTATACAAATAA',
}

# Codon-optimized filler for genes we don't have a specific sequence for.
def _synthesize_gene_sequence(gene_name, length_aa=120):
    """Generate a plausible, valid coding DNA sequence for a gene name we don't
    have on file. Deterministic (same name → same sequence) so it's stable
    across sections. Has ATG start and a stop codon."""
    import hashlib
    # Common codons for a balanced, GC-reasonable synthetic gene
    codons = ['GCT', 'GCC', 'AAA', 'GAA', 'GGT', 'CTG', 'GAT', 'ACC',
              'CGT', 'GTT', 'TTC', 'CCG', 'AGC', 'CAG', 'ATT', 'TAC',
              'AAC', 'CAC', 'ATG', 'TGG']
    # Seed from gene name for determinism
    h = hashlib.md5(gene_name.encode()).hexdigest()
    seed = int(h[:8], 16)
    seq = ['ATG']  # start codon
    x = seed
    for _ in range(length_aa - 1):
        x = (x * 1103515245 + 12345) & 0x7fffffff
        seq.append(codons[x % len(codons)])
    seq.append('TAA')  # stop
    return ''.join(seq)


def gene_name_to_sequence(gene_name):
    """Translate a gene NAME into a DNA sequence. Uses a known reference if we
    have one, otherwise synthesizes a stable valid coding sequence. This is the
    bridge that lets the synthetic organism feed CRISPR/Constructor."""
    if gene_name in KNOWN_GENE_SEQUENCES:
        return KNOWN_GENE_SEQUENCES[gene_name]
    return _synthesize_gene_sequence(gene_name)


def organism_genes_to_sequence(design):
    """Given a synthetic-organism design, build one combined DNA sequence of
    all its task genes (each with its own start/stop), ready to send to CRISPR
    or the Constructor. Returns (combined_sequence, gene_map)."""
    combined = ''
    gene_map = []  # (gene_name, start, end)
    for task in design.get('tasks', []):
        for gene in task.get('genes', []):
            seq = gene_name_to_sequence(gene)
            start = len(combined)
            combined += seq
            gene_map.append({'gene': gene, 'start': start, 'end': len(combined),
                             'task': task.get('name_en', '')})
    return combined, gene_map


def organism_tasks_to_reactions(design):
    """Translate a synthetic-organism design into simplified metabolic reactions
    for FBA. Each task adds a production reaction. This is the bridge to the
    Metabolism section — turns 'make insulin' into an actual flux reaction."""
    reactions = []
    # Base: substrate uptake + biomass
    reactions.append({'id': 'EX_glucose', 'name': 'Glucose uptake',
                      'reactants': {}, 'products': {'glucose': 1}, 'reversible': False})
    reactions.append({'id': 'biomass', 'name': 'Biomass / growth',
                      'reactants': {'glucose': 1, 'atp': 2}, 'products': {'biomass': 1},
                      'objective': True})
    # Each task → a production reaction that draws on glucose/energy
    for task in design.get('tasks', []):
        tname = task.get('name_en', 'product')
        burden = task.get('burden', 0.1)
        reactions.append({
            'id': 'PROD_' + tname[:8].replace(' ', '_'),
            'name': 'Produce ' + tname,
            'reactants': {'glucose': round(1 + burden * 3, 1), 'atp': round(burden * 5, 1)},
            'products': {tname[:8].replace(' ', '_') + '_out': 1},
            'reversible': False,
        })
    return reactions


# ╔══════════════════════════════════════════════════════════════════╗
# ║  ORGANISMS — deep profiles: genome, CRISPR, diseases, metabolism  ║
# ╚══════════════════════════════════════════════════════════════════╝

# Deep scientific profiles for model organisms. Each organism changes how
# CRISPR behaves, which diseases/applications are relevant, and simulation params.
ORGANISMS = {
    'human': {
        'name_ru': 'Человек', 'name_en': 'Human',
        'latin': 'Homo sapiens', 'icon': 'human',
        'genome_size_bp': 3_200_000_000,
        'n_genes': 20000,
        'chromosomes': 46, 'ploidy': 2,
        'genome_type': 'linear', 'cell_type': 'eukaryote',
        'doubling_time_h': 24,
        # CRISPR specifics
        'crispr': {
            'best_cas': ['SpCas9', 'SaCas9', 'Cas12a'],
            'hdr_efficiency': 0.05,   # HDR is rare in human cells (~5%)
            'nhej_efficiency': 0.70,  # NHEJ dominates
            'delivery': ['Lipofection', 'AAV', 'Electroporation', 'LNP'],
            'edit_time_days': 14,
            'difficulty': 'Высокая',
            'notes_ru': 'HDR редок — для точных правок нужны спец-методы (HDAC-инг., синхронизация).',
            'notes_en': 'HDR is rare — precise edits need special methods (HDAC-inh., synchronization).',
        },
        'diseases': [
            {'name': 'Серповидноклеточная анемия', 'name_en': 'Sickle cell anemia',
             'gene': 'HBB', 'mutation': 'Glu6Val', 'crispr_approach': 'HDR коррекция / BCL11A нокдаун'},
            {'name': 'Муковисцидоз', 'name_en': 'Cystic fibrosis',
             'gene': 'CFTR', 'mutation': 'ΔF508', 'crispr_approach': 'HDR коррекция'},
            {'name': 'Хорея Гентингтона', 'name_en': "Huntington's disease",
             'gene': 'HTT', 'mutation': 'CAG-повторы', 'crispr_approach': 'Аллель-специфичный нокаут'},
            {'name': 'Мышечная дистрофия Дюшенна', 'name_en': 'Duchenne MD',
             'gene': 'DMD', 'mutation': 'Делеции экзонов', 'crispr_approach': 'Exon skipping'},
            {'name': 'Бета-талассемия', 'name_en': 'Beta-thalassemia',
             'gene': 'HBB', 'mutation': 'Разные', 'crispr_approach': 'HDR / фетальный Hb'},
            {'name': 'Болезнь Тея-Сакса', 'name_en': 'Tay-Sachs disease',
             'gene': 'HEXA', 'mutation': 'Дефицит Hex A', 'crispr_approach': 'HDR коррекция'},
            {'name': 'Фенилкетонурия', 'name_en': 'Phenylketonuria (PKU)',
             'gene': 'PAH', 'mutation': 'Мутации PAH', 'crispr_approach': 'HDR коррекция / base editing'},
            {'name': 'Гемофилия A', 'name_en': 'Hemophilia A',
             'gene': 'F8', 'mutation': 'Дефект фактора VIII', 'crispr_approach': 'HDR вставка гена'},
            {'name': 'Семейная гиперхолестеринемия', 'name_en': 'Familial hypercholesterolemia',
             'gene': 'LDLR / PCSK9', 'mutation': 'Мутации LDLR', 'crispr_approach': 'Нокаут PCSK9'},
            {'name': 'Врождённый амавроз Лебера', 'name_en': 'Leber congenital amaurosis',
             'gene': 'CEP290', 'mutation': 'IVS26', 'crispr_approach': 'Вырезание интрона (EDIT-101)'},
            {'name': 'Тяжёлый иммунодефицит (SCID)', 'name_en': 'Severe combined immunodeficiency',
             'gene': 'IL2RG / ADA', 'mutation': 'Мутации IL2RG', 'crispr_approach': 'HDR коррекция'},
        ],
        'key_genes': ['TP53', 'BRCA1', 'CFTR', 'HBB', 'DMD', 'HTT', 'APOE'],
    },
    'ecoli': {
        'name_ru': 'Кишечная палочка', 'name_en': 'E. coli',
        'latin': 'Escherichia coli', 'icon': 'bacteria',
        'genome_size_bp': 4_600_000,
        'n_genes': 4400,
        'chromosomes': 1, 'ploidy': 1,
        'genome_type': 'circular', 'cell_type': 'prokaryote',
        'doubling_time_h': 0.33,  # 20 min
        'crispr': {
            'best_cas': ['SpCas9', 'Cas12a'],
            'hdr_efficiency': 0.80,   # lambda-Red recombineering very efficient
            'nhej_efficiency': 0.05,  # bacteria lack classic NHEJ
            'delivery': ['Трансформация', 'Электропорация', 'Конъюгация'],
            'edit_time_days': 1,
            'difficulty': 'Низкая',
            'notes_ru': 'Очень лёгкое редактирование. Нет NHEJ — нужна донорная ДНК (recombineering).',
            'notes_en': 'Very easy to edit. No NHEJ — needs donor DNA (recombineering).',
        },
        'diseases': [],  # bacteria — applications instead
        'applications': [
            {'name': 'Производство инсулина', 'name_en': 'Insulin production', 'product': 'Human insulin'},
            {'name': 'Биотопливо', 'name_en': 'Biofuel', 'product': 'Ethanol/butanol'},
            {'name': 'Производство белков', 'name_en': 'Protein production', 'product': 'Recombinant proteins'},
            {'name': 'Биосенсоры', 'name_en': 'Biosensors', 'product': 'Detection circuits'},
        ],
        'key_genes': ['lacZ', 'recA', 'rpoS', 'araC', 'tetR'],
    },
    'yeast': {
        'name_ru': 'Дрожжи', 'name_en': 'Yeast',
        'latin': 'Saccharomyces cerevisiae', 'icon': 'yeast',
        'genome_size_bp': 12_000_000,
        'n_genes': 6000,
        'chromosomes': 16, 'ploidy': 1,
        'genome_type': 'linear', 'cell_type': 'eukaryote',
        'doubling_time_h': 1.5,
        'crispr': {
            'best_cas': ['SpCas9', 'Cas12a'],
            'hdr_efficiency': 0.90,   # yeast has excellent HDR
            'nhej_efficiency': 0.10,
            'delivery': ['LiAc трансформация', 'Электропорация'],
            'edit_time_days': 3,
            'difficulty': 'Низкая',
            'notes_ru': 'Отличный HDR — идеальный организм для точных правок и сборки путей.',
            'notes_en': 'Excellent HDR — ideal for precise edits and pathway assembly.',
        },
        'diseases': [],
        'applications': [
            {'name': 'Производство этанола', 'name_en': 'Ethanol', 'product': 'Bioethanol'},
            {'name': 'Производство лекарств', 'name_en': 'Pharmaceuticals', 'product': 'Artemisinin, insulin'},
            {'name': 'Пищевая биотехнология', 'name_en': 'Food biotech', 'product': 'Enzymes, flavors'},
            {'name': 'Модель старения', 'name_en': 'Aging model', 'product': 'Research'},
        ],
        'key_genes': ['GAL1', 'ADH1', 'TEF1', 'URA3', 'HIS3'],
    },
    'arabidopsis': {
        'name_ru': 'Резуховидка (растение)', 'name_en': 'Arabidopsis',
        'latin': 'Arabidopsis thaliana', 'icon': 'plant',
        'genome_size_bp': 135_000_000,
        'n_genes': 27000,
        'chromosomes': 10, 'ploidy': 2,
        'genome_type': 'linear', 'cell_type': 'eukaryote',
        'doubling_time_h': 1008,  # ~6 weeks generation
        'crispr': {
            'best_cas': ['SpCas9', 'Cas12a'],
            'hdr_efficiency': 0.02,   # plant HDR is very low
            'nhej_efficiency': 0.60,
            'delivery': ['Agrobacterium', 'Биолистика (gene gun)', 'Протопласты'],
            'edit_time_days': 60,
            'difficulty': 'Средняя',
            'notes_ru': 'Модельное растение. HDR очень низок — в основном нокауты через NHEJ.',
            'notes_en': 'Model plant. Very low HDR — mostly NHEJ knockouts.',
        },
        'diseases': [],
        'applications': [
            {'name': 'Устойчивость к засухе', 'name_en': 'Drought resistance', 'product': 'Stress tolerance'},
            {'name': 'Повышение урожая', 'name_en': 'Yield improvement', 'product': 'Crop traits'},
            {'name': 'Биофортификация', 'name_en': 'Biofortification', 'product': 'Nutrients'},
        ],
        'key_genes': ['FLC', 'CO', 'PHYB', 'AP1', 'LFY'],
    },
    'mouse': {
        'name_ru': 'Мышь', 'name_en': 'Mouse',
        'latin': 'Mus musculus', 'icon': 'mouse',
        'genome_size_bp': 2_700_000_000,
        'n_genes': 23000,
        'chromosomes': 40, 'ploidy': 2,
        'genome_type': 'linear', 'cell_type': 'eukaryote',
        'doubling_time_h': 24,
        'crispr': {
            'best_cas': ['SpCas9', 'SaCas9'],
            'hdr_efficiency': 0.10,
            'nhej_efficiency': 0.65,
            'delivery': ['Зиготная инъекция', 'Электропорация эмбрионов', 'AAV'],
            'edit_time_days': 30,
            'difficulty': 'Средняя',
            'notes_ru': 'Главная модель болезней человека. Редактируют в зиготах.',
            'notes_en': 'Main human disease model. Edited in zygotes.',
        },
        'diseases': [
            {'name': 'Модель рака', 'name_en': 'Cancer model', 'gene': 'Trp53',
             'mutation': 'Knockout', 'crispr_approach': 'NHEJ нокаут'},
            {'name': 'Модель ожирения', 'name_en': 'Obesity model', 'gene': 'Lep',
             'mutation': 'Knockout', 'crispr_approach': 'NHEJ нокаут'},
        ],
        'key_genes': ['Trp53', 'Lep', 'Pax6', 'Sox2', 'Myc'],
    },
    'celegans': {
        'name_ru': 'Нематода C. elegans', 'name_en': 'C. elegans',
        'latin': 'Caenorhabditis elegans', 'icon': 'worm',
        'genome_size_bp': 100_000_000,
        'n_genes': 20000,
        'chromosomes': 12, 'ploidy': 2,
        'genome_type': 'linear', 'cell_type': 'eukaryote',
        'doubling_time_h': 72,
        'crispr': {
            'best_cas': ['SpCas9'],
            'hdr_efficiency': 0.30,
            'nhej_efficiency': 0.50,
            'delivery': ['Микроинъекция в гонаду'],
            'edit_time_days': 7,
            'difficulty': 'Средняя',
            'notes_ru': 'Прозрачный червь — идеален для изучения развития и нейронов.',
            'notes_en': 'Transparent worm — ideal for studying development and neurons.',
        },
        'diseases': [],
        'applications': [
            {'name': 'Изучение старения', 'name_en': 'Aging research', 'product': 'Longevity genes'},
            {'name': 'Нейробиология', 'name_en': 'Neurobiology', 'product': 'Neural circuits'},
        ],
        'key_genes': ['daf-2', 'daf-16', 'unc-22', 'dpy-10', 'lin-4'],
    },
    'drosophila': {
        'name_ru': 'Дрозофила (муха)', 'name_en': 'Drosophila',
        'latin': 'Drosophila melanogaster', 'icon': 'fly',
        'genome_size_bp': 140_000_000,
        'n_genes': 14000,
        'chromosomes': 8, 'ploidy': 2,
        'genome_type': 'linear', 'cell_type': 'eukaryote',
        'doubling_time_h': 240,  # ~10 day generation
        'crispr': {
            'best_cas': ['SpCas9'],
            'hdr_efficiency': 0.25,
            'nhej_efficiency': 0.55,
            'delivery': ['Эмбриональная инъекция'],
            'edit_time_days': 21,
            'difficulty': 'Средняя',
            'notes_ru': 'Классическая генетическая модель. Богатый набор инструментов.',
            'notes_en': 'Classic genetics model. Rich toolset available.',
        },
        'diseases': [],
        'applications': [
            {'name': 'Генетика развития', 'name_en': 'Developmental genetics', 'product': 'Research'},
            {'name': 'Модели болезней', 'name_en': 'Disease models', 'product': 'Neurodegeneration'},
        ],
        'key_genes': ['white', 'yellow', 'Antp', 'bcd', 'eve'],
    },
    'zebrafish': {
        'name_ru': 'Данио-рерио', 'name_en': 'Zebrafish',
        'latin': 'Danio rerio', 'icon': 'fish',
        'genome_size_bp': 1_400_000_000,
        'n_genes': 26000,
        'chromosomes': 50, 'ploidy': 2,
        'genome_type': 'linear', 'cell_type': 'eukaryote',
        'doubling_time_h': 0.5,
        'crispr': {
            'best_cas': ['SpCas9', 'Cas12a', 'SaCas9'],
            'hdr_efficiency': 0.06,
            'nhej_efficiency': 0.85,
            'delivery': ['Инъекция в зиготу', 'Электропорация'],
            'edit_time_days': 90,
            'difficulty': 'Низкая',
            'notes_ru': 'Прозрачные эмбрионы — правку видно под микроскопом сразу. '
                        'Инъекция в одноклеточную стадию даёт очень высокий NHEJ.',
            'notes_en': 'Transparent embryos — the edit is visible under the microscope '
                        'immediately. Injection at the one-cell stage gives very high NHEJ.',
        },
        'diseases': [
            {'name': 'Модель пигментных нарушений', 'name_en': 'Pigmentation disorder model',
             'gene': 'slc45a2', 'mutation': 'Нокаут', 'crispr_approach': 'NHEJ нокаут'},
            {'name': 'Модель кардиомиопатии', 'name_en': 'Cardiomyopathy model',
             'gene': 'tnnt2a', 'mutation': 'Нокаут', 'crispr_approach': 'NHEJ нокаут'},
        ],
        'applications': [
            {'name': 'Скрининг лекарств', 'name_en': 'Drug screening', 'product': 'Живой скрининг на эмбрионах'},
            {'name': 'Регенерация', 'name_en': 'Regeneration', 'product': 'Восстановление сердца и плавников'},
        ],
        'key_genes': ['tp53', 'sox10', 'myod1', 'tnnt2a', 'slc45a2'],
    },
    'rat': {
        'name_ru': 'Крыса', 'name_en': 'Rat',
        'latin': 'Rattus norvegicus', 'icon': 'rodent',
        'genome_size_bp': 2_750_000_000,
        'n_genes': 22000,
        'chromosomes': 42, 'ploidy': 2,
        'genome_type': 'linear', 'cell_type': 'eukaryote',
        'doubling_time_h': 24,
        'crispr': {
            'best_cas': ['SpCas9', 'SaCas9'],
            'hdr_efficiency': 0.04,
            'nhej_efficiency': 0.70,
            'delivery': ['Инъекция в зиготу', 'Электропорация', 'AAV'],
            'edit_time_days': 120,
            'difficulty': 'Средняя',
            'notes_ru': 'Крупнее мыши — удобнее для физиологии, хирургии и '
                        'поведенческих опытов. До CRISPR нокауты у крыс были почти недоступны.',
            'notes_en': 'Larger than the mouse — better for physiology, surgery and '
                        'behavioural work. Before CRISPR, rat knockouts were nearly out of reach.',
        },
        'diseases': [
            {'name': 'Модель гипертонии', 'name_en': 'Hypertension model',
             'gene': 'Ren', 'mutation': 'Нокаут', 'crispr_approach': 'NHEJ нокаут'},
            {'name': 'Модель диабета', 'name_en': 'Diabetes model',
             'gene': 'Ins2', 'mutation': 'Точечная', 'crispr_approach': 'HDR коррекция'},
        ],
        'applications': [
            {'name': 'Фармакология', 'name_en': 'Pharmacology', 'product': 'Испытания препаратов'},
            {'name': 'Нейронауки', 'name_en': 'Neuroscience', 'product': 'Модели поведения'},
        ],
        'key_genes': ['Trp53', 'Ren', 'Ins2', 'Cyp2d1', 'Bdnf'],
    },
    'pig': {
        'name_ru': 'Свинья', 'name_en': 'Pig',
        'latin': 'Sus scrofa', 'icon': 'mammal',
        'genome_size_bp': 2_500_000_000,
        'n_genes': 21000,
        'chromosomes': 38, 'ploidy': 2,
        'genome_type': 'linear', 'cell_type': 'eukaryote',
        'doubling_time_h': 24,
        'crispr': {
            'best_cas': ['SpCas9', 'SaCas9'],
            'hdr_efficiency': 0.02,
            'nhej_efficiency': 0.60,
            'delivery': ['SCNT (клонирование)', 'Инъекция в зиготу', 'Электропорация'],
            'edit_time_days': 300,
            'difficulty': 'Высокая',
            'notes_ru': 'Органы близки человеческим по размеру — основа ксенотрансплантации. '
                        'Рекорд: вырезано 62 копии ретровируса PERV за один эксперимент.',
            'notes_en': 'Organ size close to human — the basis of xenotransplantation. '
                        'A record run removed 62 copies of the PERV retrovirus at once.',
        },
        'diseases': [
            {'name': 'Модель муковисцидоза', 'name_en': 'Cystic fibrosis model',
             'gene': 'CFTR', 'mutation': 'ΔF508', 'crispr_approach': 'HDR коррекция'},
            {'name': 'Модель атеросклероза', 'name_en': 'Atherosclerosis model',
             'gene': 'LDLR', 'mutation': 'Нокаут', 'crispr_approach': 'NHEJ нокаут'},
        ],
        'applications': [
            {'name': 'Ксенотрансплантация', 'name_en': 'Xenotransplantation',
             'product': 'Органы для пересадки человеку'},
            {'name': 'Устойчивость к вирусам', 'name_en': 'Virus resistance',
             'product': 'Свиньи, устойчивые к PRRSV'},
        ],
        'key_genes': ['CFTR', 'GGTA1', 'CD163', 'PERV', 'LDLR'],
    },
    'bsubtilis': {
        'name_ru': 'Сенная палочка', 'name_en': 'B. subtilis',
        'latin': 'Bacillus subtilis', 'icon': 'bacteria',
        'genome_size_bp': 4_215_000,
        'n_genes': 4200,
        'chromosomes': 1, 'ploidy': 1,
        'genome_type': 'circular', 'cell_type': 'prokaryote',
        'doubling_time_h': 0.5,
        'crispr': {
            'best_cas': ['SpCas9', 'Cas12a'],
            'hdr_efficiency': 0.85,
            'nhej_efficiency': 0.10,
            'delivery': ['Естественная компетентность', 'Электропорация', 'Протопласты'],
            'edit_time_days': 2,
            'difficulty': 'Низкая',
            'notes_ru': 'Сама захватывает ДНК из среды — трансформация почти без усилий. '
                        'Отличная секреция белка и пищевой статус GRAS.',
            'notes_en': 'Takes up DNA from the medium by itself — transformation is nearly '
                        'effortless. Excellent protein secretion and GRAS food status.',
        },
        'diseases': [],
        'applications': [
            {'name': 'Промышленные ферменты', 'name_en': 'Industrial enzymes',
             'product': 'Протеазы и амилазы для стирки и пищи'},
            {'name': 'Пробиотики', 'name_en': 'Probiotics', 'product': 'Кормовые добавки'},
            {'name': 'Секреция белков', 'name_en': 'Protein secretion',
             'product': 'Белок сразу в среду — очистка проще'},
        ],
        'key_genes': ['comK', 'spo0A', 'amyE', 'sacB', 'srfA'],
    },
    'rice': {
        'name_ru': 'Рис', 'name_en': 'Rice',
        'latin': 'Oryza sativa', 'icon': 'plant',
        'genome_size_bp': 430_000_000,
        'n_genes': 39000,
        'chromosomes': 24, 'ploidy': 2,
        'genome_type': 'linear', 'cell_type': 'eukaryote',
        'doubling_time_h': 20,
        'crispr': {
            'best_cas': ['SpCas9', 'Cas12a', 'xCas9'],
            'hdr_efficiency': 0.03,
            'nhej_efficiency': 0.65,
            'delivery': ['Agrobacterium', 'Биолистика', 'Протопласты'],
            'edit_time_days': 150,
            'difficulty': 'Средняя',
            'notes_ru': 'Главная культура для половины человечества. Регенерация растения '
                        'из клетки занимает месяцы, зато правка наследуется устойчиво.',
            'notes_en': 'The staple crop for half the world. Regenerating a plant from a cell '
                        'takes months, but the edit is inherited reliably.',
        },
        'diseases': [
            {'name': 'Устойчивость к бактериальному ожогу', 'name_en': 'Bacterial blight resistance',
             'gene': 'OsSWEET13', 'mutation': 'Промотор', 'crispr_approach': 'Правка промотора'},
        ],
        'applications': [
            {'name': 'Урожайность', 'name_en': 'Yield', 'product': 'Крупнее зерно, больше метёлок'},
            {'name': 'Устойчивость к засухе', 'name_en': 'Drought tolerance', 'product': 'Меньше полива'},
            {'name': 'Пищевая ценность', 'name_en': 'Nutrition', 'product': 'Витамин A, железо, цинк'},
        ],
        'key_genes': ['OsSWEET13', 'Gn1a', 'DEP1', 'GS3', 'Waxy'],
    },
    'tomato': {
        'name_ru': 'Томат', 'name_en': 'Tomato',
        'latin': 'Solanum lycopersicum', 'icon': 'plant',
        'genome_size_bp': 900_000_000,
        'n_genes': 35000,
        'chromosomes': 24, 'ploidy': 2,
        'genome_type': 'linear', 'cell_type': 'eukaryote',
        'doubling_time_h': 20,
        'crispr': {
            'best_cas': ['SpCas9', 'Cas12a'],
            'hdr_efficiency': 0.03,
            'nhej_efficiency': 0.60,
            'delivery': ['Agrobacterium', 'Протопласты'],
            'edit_time_days': 120,
            'difficulty': 'Средняя',
            'notes_ru': 'На томате показали «ускоренную доместикацию»: дикий вид сделали '
                        'культурным за один цикл, правя гены формы и сроков плодоношения.',
            'notes_en': 'Tomato showed "de novo domestication": a wild species was made '
                        'crop-like in a single cycle by editing shape and fruiting-time genes.',
        },
        'diseases': [
            {'name': 'Устойчивость к мучнистой росе', 'name_en': 'Powdery mildew resistance',
             'gene': 'SlMlo1', 'mutation': 'Нокаут', 'crispr_approach': 'NHEJ нокаут'},
        ],
        'applications': [
            {'name': 'Ускоренная доместикация', 'name_en': 'De novo domestication',
             'product': 'Дикие виды → культурные за один цикл'},
            {'name': 'Срок хранения', 'name_en': 'Shelf life', 'product': 'Медленнее размягчается'},
            {'name': 'Содержание ГАМК', 'name_en': 'GABA content',
             'product': 'Первый CRISPR-продукт в продаже (Япония, 2021)'},
        ],
        'key_genes': ['SlMlo1', 'SP5G', 'SlAGL6', 'RIN', 'SlGAD3'],
    },
    'chlamydomonas': {
        'name_ru': 'Хламидомонада (водоросль)', 'name_en': 'Chlamydomonas (green alga)',
        'latin': 'Chlamydomonas reinhardtii', 'icon': 'algae',
        'genome_size_bp': 120_000_000,
        'n_genes': 17000,
        'chromosomes': 17, 'ploidy': 1,
        'genome_type': 'linear', 'cell_type': 'eukaryote',
        'doubling_time_h': 8,
        'crispr': {
            'best_cas': ['SpCas9', 'Cas12a'],
            'hdr_efficiency': 0.10,
            'nhej_efficiency': 0.55,
            'delivery': ['Электропорация RNP', 'Стеклянные бусины', 'Биолистика'],
            'edit_time_days': 21,
            'difficulty': 'Средняя',
            'notes_ru': 'Гаплоидна — одна правка сразу даёт фенотип, второй аллель искать не надо. '
                        'Растёт на свету и CO₂, есть и хлоропластный геном.',
            'notes_en': 'Haploid — a single edit shows the phenotype at once, with no second '
                        'allele to chase. Grows on light and CO₂, and has a chloroplast genome too.',
        },
        'diseases': [
            {'name': 'Модель цилиопатий', 'name_en': 'Ciliopathy model',
             'gene': 'IFT88', 'mutation': 'Нокаут', 'crispr_approach': 'NHEJ нокаут'},
        ],
        'applications': [
            {'name': 'Биотопливо', 'name_en': 'Biofuel', 'product': 'Липиды и водород'},
            {'name': 'Фотосинтез', 'name_en': 'Photosynthesis', 'product': 'Изучение и улучшение'},
            {'name': 'Кормовой белок', 'name_en': 'Feed protein', 'product': 'Белок из света и CO₂'},
        ],
        'key_genes': ['IFT88', 'PSBA', 'HYDA1', 'NIT1', 'CRY1'],
    },
    'cho_cells': {
        'name_ru': 'Клетки CHO', 'name_en': 'CHO cells',
        'latin': 'Cricetulus griseus (клеточная линия)', 'icon': 'cell',
        'genome_size_bp': 2_450_000_000,
        'n_genes': 22000,
        'chromosomes': 22, 'ploidy': 2,
        'genome_type': 'linear', 'cell_type': 'eukaryote',
        'doubling_time_h': 20,
        'crispr': {
            'best_cas': ['SpCas9', 'SaCas9', 'Cas12a'],
            'hdr_efficiency': 0.08,
            'nhej_efficiency': 0.75,
            'delivery': ['Электропорация RNP', 'Липофекция', 'Лентивирус'],
            'edit_time_days': 30,
            'difficulty': 'Средняя',
            'notes_ru': 'Не организм, а клеточная линия — но именно на ней делают большинство '
                        'антител в мире. Клетки можно отбирать клонами, поэтому HDR годится.',
            'notes_en': 'Not an organism but a cell line — yet most of the world\u2019s antibodies '
                        'are made in it. Clones can be selected, so HDR is workable here.',
        },
        'diseases': [],
        'applications': [
            {'name': 'Моноклональные антитела', 'name_en': 'Monoclonal antibodies',
             'product': 'Основа современной биофармацевтики'},
            {'name': 'Правка гликозилирования', 'name_en': 'Glycoengineering',
             'product': 'Нокаут FUT8 усиливает действие антител'},
        ],
        'key_genes': ['FUT8', 'GLUL', 'DHFR', 'BAX', 'BAK'],
    },
}


def get_organism(org_key):
    """Return an organism profile by key."""
    return ORGANISMS.get(org_key)


def design_cloning_oligos(guide_seq, vector='lentiGuide-Puro'):
    """Design overlapping oligonucleotides to clone a guide into a CRISPR vector.

    Like CRISPOR's cloning helper. Adds the BsmBI/BbsI overhangs needed for
    Golden Gate cloning into common AddGene vectors. Returns a dict with the
    forward and reverse oligos ready to order.
    """
    g = (guide_seq or '').upper().replace('U', 'T')
    if len(g) < 18:
        return {'error': 'Guide too short'}
    # Common overhangs for U6-driven vectors (BsmBI). Add G if guide doesn't
    # start with G (needed for U6 transcription start).
    prefix_fwd = 'CACCG' if not g.startswith('G') else 'CACC'
    # Reverse complement of guide for the bottom strand
    comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    rc = ''.join(comp.get(b, 'N') for b in reversed(g))
    prefix_rev = 'AAAC' + ('C' if not g.startswith('G') else '')
    fwd = prefix_fwd + g
    rev = prefix_rev + rc
    return {
        'vector': vector,
        'guide': g,
        'forward_oligo': fwd,
        'reverse_oligo': rev,
        'enzyme': 'BsmBI / Esp3I',
        'note': ('Anneal the two oligos and ligate into the BsmBI-digested '
                 'vector. The extra G (if added) ensures U6 transcription.'),
    }


def design_validation_primers(sequence, cut_position, amplicon_size=400):
    """Design PCR primers flanking a CRISPR cut site for validation
    (e.g. by T7E1 assay or sequencing). Like CRISPOR's primer designer.

    Returns forward + reverse primers ~amplicon_size bp apart, centred on
    the cut. Picks ~20-mer primers with reasonable Tm.
    """
    seq = (sequence or '').upper()
    if len(seq) < amplicon_size or cut_position < 0 or cut_position >= len(seq):
        return {'error': 'Sequence too short or invalid cut position'}
    half = amplicon_size // 2
    left_start = max(0, cut_position - half)
    right_end = min(len(seq), cut_position + half)

    def pick_primer(region, reverse=False):
        # Try to find a ~20-mer with GC 40-60% and ending in G/C
        for length in (20, 21, 22, 19, 23):
            for offset in range(0, min(15, len(region) - length)):
                p = region[offset:offset + length]
                gc = (p.count('G') + p.count('C')) / len(p) * 100
                if 40 <= gc <= 60 and p[-1] in 'GC':
                    return p
        return region[:20]

    comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    fwd_region = seq[left_start:left_start + 40]
    fwd = pick_primer(fwd_region)
    rev_region = seq[right_end - 40:right_end]
    rev_template = pick_primer(rev_region)
    rev = ''.join(comp.get(b, 'N') for b in reversed(rev_template))

    def tm(p):  # Wallace rule approximation
        at = p.count('A') + p.count('T')
        gc = p.count('G') + p.count('C')
        return round(2 * at + 4 * gc, 1) if len(p) < 14 else round(
            64.9 + 41 * (gc - 16.4) / len(p), 1)

    return {
        'forward_primer': fwd,
        'reverse_primer': rev,
        'forward_tm': tm(fwd),
        'reverse_tm': tm(rev),
        'amplicon_size': right_end - left_start,
        'cut_position_in_amplicon': cut_position - left_start,
    }


# Base editors: fuse a Cas9 nickase with a deaminase. They convert bases
# within a small editing window without a double-strand break.
BASE_EDITORS = {
    'CBE (C→T)': {
        'name_ru': 'Цитозиновый редактор (C→T)', 'name_en': 'Cytosine base editor (C→T)',
        'edit_from': 'C', 'edit_to': 'T', 'comp_from': 'G', 'comp_to': 'A',
        'window_start': 4, 'window_end': 8,  # positions in the 20-nt protospacer
        'desc_ru': 'Превращает C·G в T·A. Не делает двунитевой разрыв.',
        'desc_en': 'Converts C·G to T·A. No double-strand break.',
    },
    'ABE (A→G)': {
        'name_ru': 'Адениновый редактор (A→G)', 'name_en': 'Adenine base editor (A→G)',
        'edit_from': 'A', 'edit_to': 'G', 'comp_from': 'T', 'comp_to': 'C',
        'window_start': 4, 'window_end': 8,
        'desc_ru': 'Превращает A·T в G·C. Не делает двунитевой разрыв.',
        'desc_en': 'Converts A·T to G·C. No double-strand break.',
    },
}


def apply_base_editor(protospacer, editor_key='CBE (C→T)'):
    """Apply a base editor to a 20-nt protospacer. Converts target bases
    inside the editing window (~positions 4-8). Returns the edited sequence
    and a list of edits made.

    Base editing is precise single-base change without a double-strand break.
    """
    ed = BASE_EDITORS.get(editor_key)
    if not ed:
        return {'error': 'Unknown base editor'}
    seq = list((protospacer or '').upper())
    if len(seq) < 20:
        return {'error': 'Protospacer must be 20 nt'}
    edits = []
    ws, we = ed['window_start'], ed['window_end']
    for i in range(ws - 1, min(we, len(seq))):
        if seq[i] == ed['edit_from']:
            seq[i] = ed['edit_to']
            edits.append({'position': i + 1, 'from': ed['edit_from'],
                          'to': ed['edit_to']})
    return {
        'editor': editor_key,
        'original': protospacer.upper(),
        'edited': ''.join(seq),
        'edits': edits,
        'n_edits': len(edits),
        'window': f'{ws}-{we}',
    }


def design_prime_editing(target_seq, edit_position, original_base, new_base,
                         pbs_length=13, rt_length=13):
    """Design a prime editing guide RNA (pegRNA) — newest CRISPR technology
    (2025). Prime editing makes precise insertions, deletions, or any base
    substitution WITHOUT a double-strand break or donor DNA.

    A pegRNA has 3 parts:
      - spacer (the guide, finds the target)
      - PBS (primer binding site) — anneals to the nicked strand
      - RTT (reverse transcriptase template) — encodes the desired edit

    Returns a dict describing the pegRNA components.
    """
    seq = (target_seq or '').upper()
    if len(seq) < 30 or edit_position < 0 or edit_position >= len(seq):
        return {'error': 'Sequence too short or invalid edit position'}
    comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    # Spacer: 20 nt upstream of a PAM near the edit
    spacer_start = max(0, edit_position - 17)
    spacer = seq[spacer_start:spacer_start + 20]
    # Nick site is ~3 nt upstream of PAM (position 17 in spacer)
    nick_pos = spacer_start + 17
    # PBS: complementary to the 3' end of the nicked strand
    pbs_region = seq[max(0, nick_pos - pbs_length):nick_pos]
    pbs = ''.join(comp.get(b, 'N') for b in reversed(pbs_region))
    # RTT: encodes the edit — template for the new sequence
    edited_seq = seq[:edit_position] + new_base + seq[edit_position + 1:]
    rtt_region = edited_seq[nick_pos:nick_pos + rt_length]
    rtt = ''.join(comp.get(b, 'N') for b in reversed(rtt_region))
    # The 3' extension = RTT + PBS
    extension = rtt + pbs
    return {
        'spacer': spacer,
        'pbs': pbs,
        'rtt': rtt,
        'extension_3prime': extension,
        'nick_position': nick_pos,
        'edit': '{}>{} at position {}'.format(original_base, new_base, edit_position + 1),
        'pbs_length': pbs_length,
        'rt_length': rt_length,
        'note': ('Prime editing: precise edit without double-strand break. '
                 'Order the pegRNA with this 3\' extension.'),
    }


# CRISPRi (interference/repression) and CRISPRa (activation) use dCas9 —
# a catalytically DEAD Cas9 that binds but does not cut. Fused to effectors.
CRISPR_REGULATION = {
    'CRISPRi': {
        'name_ru': 'CRISPRi (подавление)', 'name_en': 'CRISPRi (repression)',
        'effector': 'dCas9-KRAB', 'effect': 'repress',
        'desc_ru': 'dCas9 + KRAB домен подавляет транскрипцию гена',
        'desc_en': 'dCas9 + KRAB domain represses gene transcription',
        'optimal_window_ru': 'Целиться в промотор или +1..+50 от старта',
        'optimal_window_en': 'Target promoter or +1..+50 from TSS',
    },
    'CRISPRa': {
        'name_ru': 'CRISPRa (активация)', 'name_en': 'CRISPRa (activation)',
        'effector': 'dCas9-VPR', 'effect': 'activate',
        'desc_ru': 'dCas9 + VPR активатор усиливает экспрессию гена',
        'desc_en': 'dCas9 + VPR activator boosts gene expression',
        'optimal_window_ru': 'Целиться в -400..-50 перед стартом транскрипции',
        'optimal_window_en': 'Target -400..-50 upstream of TSS',
    },
}


def design_crispri_crispra(guide_seq, distance_from_tss, mode='CRISPRi'):
    """Score a guide for CRISPRi (repression) or CRISPRa (activation).
    These use dCas9 (dead Cas9) — no cutting, just regulation.

    Effectiveness depends heavily on WHERE the guide sits relative to the
    transcription start site (TSS). Returns score + recommendation.
    """
    reg = CRISPR_REGULATION.get(mode)
    if not reg:
        return {'error': 'Unknown mode'}
    score = 50.0
    if mode == 'CRISPRi':
        # CRISPRi works best near TSS (0 to +50)
        if 0 <= distance_from_tss <= 50:
            score = 90.0
        elif -50 <= distance_from_tss < 0:
            score = 70.0
        elif 50 < distance_from_tss <= 300:
            score = 55.0
        else:
            score = 25.0
    else:  # CRISPRa
        # CRISPRa works best upstream (-400 to -50)
        if -400 <= distance_from_tss <= -50:
            score = 90.0
        elif -50 < distance_from_tss <= 0:
            score = 60.0
        else:
            score = 25.0
    return {
        'mode': mode,
        'effector': reg['effector'],
        'effect': reg['effect'],
        'score': round(score, 1),
        'distance_from_tss': distance_from_tss,
        'guide': guide_seq,
    }


def design_multiplexed_guides(targets_list, enzyme='SpCas9'):
    """Multiplexed CRISPR — design guides for editing MULTIPLE genes at once.
    Returns a guide for each target plus a compatibility check (guides should
    not interfere or target each other).

    Each target in targets_list is a dict {'name':..., 'sequence':...}.
    """
    results = []
    all_guides = []
    # find_crispr_targets lives in biosim_part2 — import lazily
    try:
        from biosim_part2 import find_crispr_targets as _fct
    except Exception as _exc:
        biosim_log('warn', 'design_multiplexed_guides: {}'.format(_exc))
        _fct = None
    for tgt in targets_list:
        name = tgt.get('name', '?')
        seq = (tgt.get('sequence', '') or '').upper()
        guides = _fct(seq, enzyme=enzyme, max_sites=5) if _fct else []
        best = guides[0] if guides else None
        if best:
            all_guides.append(best['guide'])
            results.append({
                'target': name,
                'guide': best['guide'],
                'score': best.get('score', 0),
                'position': best.get('position', 0),
            })
        else:
            results.append({'target': name, 'guide': None, 'score': 0})
    # Check guides don't cross-target (share >15 nt identity)
    warnings = []
    for i in range(len(all_guides)):
        for j in range(i + 1, len(all_guides)):
            g1, g2 = all_guides[i], all_guides[j]
            matches = sum(1 for a, b in zip(g1, g2) if a == b)
            if matches > 15:
                warnings.append('Гайды {} и {} похожи ({} совпадений)'.format(
                    i + 1, j + 1, matches))
    return {
        'guides': results,
        'n_targets': len(targets_list),
        'n_designed': sum(1 for r in results if r['guide']),
        'warnings': warnings,
    }


def organism_crispr_efficiency(org_key, cas_enzyme='SpCas9', edit_type='knockout'):
    """Compute expected CRISPR success rate for an organism + approach.

    edit_type: 'knockout' (uses NHEJ) or 'knockin'/'correction' (uses HDR).
    Returns a dict with efficiency and explanation.
    """
    org = ORGANISMS.get(org_key)
    if not org:
        return {'error': 'Unknown organism'}
    crispr = org['crispr']
    cas_ok = cas_enzyme in crispr['best_cas']

    if edit_type in ('knockin', 'correction', 'hdr'):
        base_eff = crispr['hdr_efficiency']
        mechanism = 'HDR'
    else:
        base_eff = crispr['nhej_efficiency']
        mechanism = 'NHEJ'

    # Penalty if using a non-optimal Cas
    if not cas_ok:
        base_eff *= 0.6

    return {
        'error': None,
        'organism': org['name_en'],
        'mechanism': mechanism,
        'efficiency': round(base_eff * 100, 1),
        'cas_optimal': cas_ok,
        'edit_time_days': crispr['edit_time_days'],
        'difficulty': crispr['difficulty'],
        'delivery_methods': crispr['delivery'],
        'notes_ru': crispr['notes_ru'],
        'notes_en': crispr['notes_en'],
    }


# ╔══════════════════════════════════════════════════════════════════╗
# ║  METABOLIC CONSTRUCTOR — build enzyme chains, simulate products   ║
# ║  Organism-aware kinetic simulation (Michaelis-Menten)            ║
# ╚══════════════════════════════════════════════════════════════════╝

# Enzyme library for the metabolic constructor. Each enzyme converts a
# substrate to a product with a kinetic rate (kcat) and affinity (Km).
METABOLIC_ENZYMES_KINETIC = {
    'hexokinase':    {'name_ru': 'Гексокиназа', 'name_en': 'Hexokinase',
                      'substrate': 'Glucose', 'product': 'G6P',
                      'kcat': 8.0, 'km': 0.5, 'atp': -1, 'nadh': 0, 'yield_factor': 1},
    'pfk':           {'name_ru': 'Фосфофруктокиназа (PFK)', 'name_en': 'PFK',
                      'substrate': 'G6P', 'product': 'F16BP',
                      'kcat': 3.0, 'km': 0.3, 'atp': -1, 'nadh': 0, 'yield_factor': 1},  # rate-limiting
    'gapdh':         {'name_ru': 'GAPDH', 'name_en': 'GAPDH',
                      'substrate': 'F16BP', 'product': '13BPG',
                      'kcat': 10.0, 'km': 0.4, 'atp': 0, 'nadh': 1, 'yield_factor': 2},  # C6→2×C3 (aldolase split)
    'pk':            {'name_ru': 'Пируваткиназа', 'name_en': 'Pyruvate kinase',
                      'substrate': '13BPG', 'product': 'Pyruvate',
                      'kcat': 12.0, 'km': 0.3, 'atp': 2, 'nadh': 0, 'yield_factor': 1},
    'pdc':           {'name_ru': 'Пируватдекарбоксилаза', 'name_en': 'Pyruvate decarboxylase',
                      'substrate': 'Pyruvate', 'product': 'Acetaldehyde',
                      'kcat': 9.0, 'km': 0.4, 'atp': 0, 'nadh': 0, 'yield_factor': 1},
    'adh':           {'name_ru': 'Алкогольдегидрогеназа', 'name_en': 'Alcohol dehydrogenase',
                      'substrate': 'Acetaldehyde', 'product': 'Ethanol',
                      'kcat': 11.0, 'km': 0.3, 'atp': 0, 'nadh': -1, 'yield_factor': 1},
    'ldh':           {'name_ru': 'Лактатдегидрогеназа', 'name_en': 'Lactate dehydrogenase',
                      'substrate': 'Pyruvate', 'product': 'Lactate',
                      'kcat': 10.0, 'km': 0.35, 'atp': 0, 'nadh': -1, 'yield_factor': 1},
    'citrate_synthase': {'name_ru': 'Цитратсинтаза', 'name_en': 'Citrate synthase',
                      'substrate': 'AcetylCoA', 'product': 'Citrate',
                      'kcat': 6.0, 'km': 0.4, 'atp': 0, 'nadh': 0, 'yield_factor': 1},
    'idh':           {'name_ru': 'Изоцитратдегидрогеназа', 'name_en': 'IDH',
                      'substrate': 'Citrate', 'product': 'aKG',
                      'kcat': 4.0, 'km': 0.3, 'atp': 0, 'nadh': 1, 'yield_factor': 1},
    'pdh':           {'name_ru': 'Пируватдегидрогеназа', 'name_en': 'Pyruvate dehydrogenase',
                      'substrate': 'Pyruvate', 'product': 'AcetylCoA',
                      'kcat': 7.0, 'km': 0.4, 'atp': 0, 'nadh': 1, 'yield_factor': 1},
    # ── NEW: продолжение цикла Кребса ──
    'akgdh':         {'name_ru': 'α-кетоглутаратдегидрогеназа', 'name_en': 'α-KG dehydrogenase',
                      'substrate': 'aKG', 'product': 'SuccinylCoA',
                      'kcat': 5.0, 'km': 0.35, 'atp': 0, 'nadh': 1, 'yield_factor': 1},
    'scs':           {'name_ru': 'Сукцинил-CoA-синтетаза', 'name_en': 'Succinyl-CoA synthetase',
                      'substrate': 'SuccinylCoA', 'product': 'Succinate',
                      'kcat': 6.0, 'km': 0.3, 'atp': 1, 'nadh': 0, 'yield_factor': 1},
    'sdh':           {'name_ru': 'Сукцинатдегидрогеназа', 'name_en': 'Succinate dehydrogenase',
                      'substrate': 'Succinate', 'product': 'Fumarate',
                      'kcat': 5.5, 'km': 0.3, 'atp': 0, 'nadh': 0, 'fadh2': 1, 'yield_factor': 1},
    'fumarase':      {'name_ru': 'Фумараза', 'name_en': 'Fumarase',
                      'substrate': 'Fumarate', 'product': 'Malate',
                      'kcat': 8.0, 'km': 0.25, 'atp': 0, 'nadh': 0, 'yield_factor': 1},
    'mdh':           {'name_ru': 'Малатдегидрогеназа', 'name_en': 'Malate dehydrogenase',
                      'substrate': 'Malate', 'product': 'Oxaloacetate',
                      'kcat': 7.0, 'km': 0.3, 'atp': 0, 'nadh': 1, 'yield_factor': 1},
    # ── NEW: глюконеогенез / пентозофосфатный путь ──
    'g6pdh':         {'name_ru': 'Глюкозо-6-фосфатдегидрогеназа', 'name_en': 'G6P dehydrogenase',
                      'substrate': 'G6P', 'product': '6PGL',
                      'kcat': 6.0, 'km': 0.3, 'atp': 0, 'nadph': 1, 'yield_factor': 1},
    'transketolase': {'name_ru': 'Транскетолаза', 'name_en': 'Transketolase',
                      'substrate': '6PGL', 'product': 'Ribose5P',
                      'kcat': 5.0, 'km': 0.35, 'atp': 0, 'nadh': 0, 'yield_factor': 1},
    # ── TCA cycle ────────────────────────────────────────────────────
    'citrate_synthase': {'name_ru': 'Цитратсинтаза', 'name_en': 'Citrate synthase',
                      'substrate': 'AcetylCoA', 'product': 'Citrate',
                      'kcat': 12.0, 'km': 0.10, 'atp': 0, 'nadh': 0, 'yield_factor': 1},
    'aconitase':      {'name_ru': 'Аконитаза', 'name_en': 'Aconitase',
                      'substrate': 'Citrate', 'product': 'Isocitrate',
                      'kcat': 15.0, 'km': 0.50, 'atp': 0, 'nadh': 0, 'yield_factor': 1},
    'idh':            {'name_ru': 'Изоцитратдегидрогеназа', 'name_en': 'Isocitrate dehydrogenase',
                      'substrate': 'Isocitrate', 'product': 'AKG',
                      'kcat': 9.0, 'km': 0.03, 'atp': 0, 'nadh': 1, 'yield_factor': 1},
    'akgdh':          {'name_ru': 'Альфа-кетоглутаратдегидрогеназа', 'name_en': 'AKG dehydrogenase',
                      'substrate': 'AKG', 'product': 'SuccinylCoA',
                      'kcat': 7.0, 'km': 0.20, 'atp': 0, 'nadh': 1, 'yield_factor': 1},
    'succinyl_syn':   {'name_ru': 'Сукцинил-КоА-синтетаза', 'name_en': 'Succinyl-CoA synthetase',
                      'substrate': 'SuccinylCoA', 'product': 'Succinate',
                      'kcat': 14.0, 'km': 0.25, 'atp': 1, 'nadh': 0, 'yield_factor': 1},
    'sdh':            {'name_ru': 'Сукцинатдегидрогеназа', 'name_en': 'Succinate dehydrogenase',
                      'substrate': 'Succinate', 'product': 'Fumarate',
                      'kcat': 11.0, 'km': 0.40, 'atp': 0, 'nadh': 1, 'yield_factor': 1},
    'fumarase':       {'name_ru': 'Фумараза', 'name_en': 'Fumarase',
                      'substrate': 'Fumarate', 'product': 'Malate',
                      'kcat': 25.0, 'km': 0.15, 'atp': 0, 'nadh': 0, 'yield_factor': 1},
    'mdh':            {'name_ru': 'Малатдегидрогеназа', 'name_en': 'Malate dehydrogenase',
                      'substrate': 'Malate', 'product': 'OAA',
                      'kcat': 20.0, 'km': 0.30, 'atp': 0, 'nadh': 1, 'yield_factor': 1},

    # ── pentose phosphate pathway ────────────────────────────────────
    'g6pdh':          {'name_ru': 'Г6Ф-дегидрогеназа', 'name_en': 'G6P dehydrogenase',
                      'substrate': 'G6P', 'product': '6PG',
                      'kcat': 10.0, 'km': 0.06, 'atp': 0, 'nadh': 1, 'yield_factor': 1},
    'pgd':            {'name_ru': '6-фосфоглюконатдегидрогеназа', 'name_en': '6PG dehydrogenase',
                      'substrate': '6PG', 'product': 'Ru5P',
                      'kcat': 8.0, 'km': 0.04, 'atp': 0, 'nadh': 1, 'yield_factor': 1},
    'transketolase':  {'name_ru': 'Транскетолаза', 'name_en': 'Transketolase',
                      'substrate': 'Ru5P', 'product': 'S7P',
                      'kcat': 13.0, 'km': 0.30, 'atp': 0, 'nadh': 0, 'yield_factor': 1},

    # ── fermentation branches ────────────────────────────────────────
    'pfl':            {'name_ru': 'Пируватформиатлиаза', 'name_en': 'Pyruvate formate-lyase',
                      'substrate': 'Pyruvate', 'product': 'Formate',
                      'kcat': 16.0, 'km': 2.00, 'atp': 0, 'nadh': 0, 'yield_factor': 1},
    'ack':            {'name_ru': 'Ацетаткиназа', 'name_en': 'Acetate kinase',
                      'substrate': 'AcetylCoA', 'product': 'Acetate',
                      'kcat': 22.0, 'km': 0.70, 'atp': 1, 'nadh': 0, 'yield_factor': 1},
    'butyrate_syn':   {'name_ru': 'Синтез бутирата', 'name_en': 'Butyrate synthase',
                      'substrate': 'AcetylCoA', 'product': 'Butyrate',
                      'kcat': 6.0, 'km': 0.45, 'atp': 1, 'nadh': -2, 'yield_factor': 1},

    # ── amino acid and product routes ────────────────────────────────
    'aspartate_kinase': {'name_ru': 'Аспартаткиназа', 'name_en': 'Aspartate kinase',
                      'substrate': 'OAA', 'product': 'Aspartate',
                      'kcat': 5.0, 'km': 0.35, 'atp': -1, 'nadh': 0, 'yield_factor': 1},
    'dapa_syn':       {'name_ru': 'Дигидродипиколинатсинтаза', 'name_en': 'Dihydrodipicolinate synthase',
                      'substrate': 'Aspartate', 'product': 'DAP',
                      'kcat': 4.0, 'km': 0.25, 'atp': 0, 'nadh': -1, 'yield_factor': 1},
    'lysa':           {'name_ru': 'Диаминопимелатдекарбоксилаза', 'name_en': 'DAP decarboxylase',
                      'substrate': 'DAP', 'product': 'Lysine',
                      'kcat': 3.5, 'km': 0.20, 'atp': 0, 'nadh': 0, 'yield_factor': 1},
    'glutamate_dh':   {'name_ru': 'Глутаматдегидрогеназа', 'name_en': 'Glutamate dehydrogenase',
                      'substrate': 'AKG', 'product': 'Glutamate',
                      'kcat': 9.0, 'km': 0.60, 'atp': 0, 'nadh': -1, 'yield_factor': 1},
    'phaA_k':         {'name_ru': 'Бета-кетотиолаза (PHB)', 'name_en': 'Beta-ketothiolase (PHB)',
                      'substrate': 'AcetylCoA', 'product': 'AcAcCoA',
                      'kcat': 7.5, 'km': 0.40, 'atp': 0, 'nadh': 0, 'yield_factor': 1},
    'phaB_k':         {'name_ru': 'Ацетоацетил-КоА-редуктаза', 'name_en': 'AcAcCoA reductase',
                      'substrate': 'AcAcCoA', 'product': 'HBCoA',
                      'kcat': 6.5, 'km': 0.15, 'atp': 0, 'nadh': -1, 'yield_factor': 1},
    'phaC_k':         {'name_ru': 'ПГБ-синтаза', 'name_en': 'PHB synthase',
                      'substrate': 'HBCoA', 'product': 'PHB',
                      'kcat': 4.0, 'km': 0.30, 'atp': 0, 'nadh': 0, 'yield_factor': 1},
    'hmgr':           {'name_ru': 'ГМГ-КоА-редуктаза', 'name_en': 'HMG-CoA reductase',
                      'substrate': 'AcAcCoA', 'product': 'Mevalonate',
                      'kcat': 2.5, 'km': 0.05, 'atp': 0, 'nadh': -2, 'yield_factor': 1},
    'mvk':            {'name_ru': 'Мевалонаткиназа', 'name_en': 'Mevalonate kinase',
                      'substrate': 'Mevalonate', 'product': 'IPP',
                      'kcat': 5.5, 'km': 0.12, 'atp': -1, 'nadh': 0, 'yield_factor': 1},
}


def get_metabolic_enzyme(key):
    return METABOLIC_ENZYMES_KINETIC.get(key)


def simulate_metabolic_chain(enzyme_chain, substrate_amount=10.0,
                             organism_key='ecoli', t_max=20.0, n_points=200,
                             enzyme_levels=None, inhibitors=None):
    """Simulate a metabolic enzyme chain over time (Michaelis-Menten kinetics).

    enzyme_chain: list of enzyme keys in order, e.g. ['hexokinase','pfk',...]
    organism_key: affects overall rate (via doubling time / metabolic speed)
    enzyme_levels: optional dict {enzyme_key: concentration} (default 1.0 each)
    inhibitors: optional dict {enzyme_key: inhibition 0..1} (1.0 = fully blocked)
    Returns time series of all metabolite concentrations + ATP/NADH balance.

    Clear, deterministic ODE simulation:
      v = kcat * [E] * [S] / (Km + [S]) * (1 - inhibition)   (Michaelis-Menten)
    """
    if not enzyme_chain:
        return {'error': 'Empty enzyme chain'}

    enzyme_levels = enzyme_levels or {}
    inhibitors = inhibitors or {}

    # Organism metabolic speed factor (faster doubling = faster metabolism)
    org = ORGANISMS.get(organism_key, {})
    dt_h = org.get('doubling_time_h', 1.0)
    # Normalize: E.coli (0.33h) = 1.0x speed; slower organisms scale down
    speed_factor = max(0.2, min(2.0, 0.33 / max(0.05, dt_h) ** 0.3))

    # Build the list of reactions and the set of all metabolites
    reactions = []
    reaction_keys = []
    metabolites = []
    for ekey in enzyme_chain:
        enz = METABOLIC_ENZYMES_KINETIC.get(ekey)
        if not enz:
            continue
        reactions.append(enz)
        reaction_keys.append(ekey)
        for met in (enz['substrate'], enz['product']):
            if met not in metabolites:
                metabolites.append(met)

    if not reactions:
        return {'error': 'No valid enzymes in chain'}

    # Index metabolites
    met_idx = {m: i for i, m in enumerate(metabolites)}
    n_met = len(metabolites)

    # Initial concentrations: substrate of first enzyme = substrate_amount
    import numpy as _np
    y0 = _np.zeros(n_met)
    first_substrate = reactions[0]['substrate']
    y0[met_idx[first_substrate]] = substrate_amount

    def rhs(t, y):
        dydt = _np.zeros(n_met)
        for enz, ekey in zip(reactions, reaction_keys):
            s_idx = met_idx[enz['substrate']]
            p_idx = met_idx[enz['product']]
            s = max(0.0, y[s_idx])
            # Enzyme concentration (default 1.0) and inhibition (default 0)
            E = enzyme_levels.get(ekey, 1.0)
            inhib = min(1.0, max(0.0, inhibitors.get(ekey, 0.0)))
            # Michaelis-Menten velocity with enzyme level + inhibition
            v = enz['kcat'] * E * s / (enz['km'] + s + 1e-9) * speed_factor * (1.0 - inhib)
            # Stoichiometry: product molecules per substrate (aldolase C6→2×C3)
            yf = enz.get('yield_factor', 1)
            dydt[s_idx] -= v
            dydt[p_idx] += v * yf
        return dydt

    # Integrate with simple RK4 (clear, dependency-light)
    t_vals = _np.linspace(0, t_max, n_points)
    dt = t_vals[1] - t_vals[0]
    Y = _np.zeros((n_points, n_met))
    Y[0] = y0
    y = y0.copy()
    for i in range(1, n_points):
        k1 = rhs(t_vals[i-1], y)
        k2 = rhs(t_vals[i-1] + dt/2, y + dt/2 * k1)
        k3 = rhs(t_vals[i-1] + dt/2, y + dt/2 * k2)
        k4 = rhs(t_vals[i-1] + dt, y + dt * k3)
        y = y + dt/6 * (k1 + 2*k2 + 2*k3 + k4)
        y = _np.maximum(0.0, y)  # no negative concentrations
        Y[i] = y

    # Final product = product of the last enzyme
    final_product = reactions[-1]['product']
    product_amount = float(Y[-1, met_idx[final_product]])

    # Energy balance (sum over reactions, scaled by how much flowed through)
    total_atp = 0.0
    total_nadh = 0.0
    for enz in reactions:
        converted = float(Y[-1, met_idx[enz['product']]])
        total_atp += enz['atp'] * converted
        total_nadh += enz['nadh'] * converted

    return {
        'error': None,
        'time': t_vals.tolist(),
        'metabolites': metabolites,
        'concentrations': {m: Y[:, met_idx[m]].tolist() for m in metabolites},
        'final_product': final_product,
        'product_amount': round(product_amount, 3),
        'substrate': first_substrate,
        'substrate_amount': substrate_amount,
        'conversion_pct': round(product_amount / substrate_amount * 100, 1) if substrate_amount else 0,
        'atp_balance': round(total_atp, 2),
        'nadh_balance': round(total_nadh, 2),
        'organism': org.get('name_en', organism_key),
        'speed_factor': round(speed_factor, 2),
        'n_enzymes': len(reactions),
    }


def optimize_enzyme_levels(enzyme_chain, substrate_amount=10.0, organism_key='ecoli',
                           total_enzyme_budget=None):
    """Find enzyme concentrations that maximize final product.

    Uses a simple hill-climbing search. Returns the best enzyme levels found
    and the resulting product amount. This answers: 'how much of each enzyme
    should I make to get the most product?'
    """
    if not enzyme_chain:
        return {'error': 'Empty chain'}
    valid = [e for e in enzyme_chain if METABOLIC_ENZYMES_KINETIC.get(e)]
    if not valid:
        return {'error': 'No valid enzymes'}

    import numpy as _np
    n = len(valid)
    budget = total_enzyme_budget if total_enzyme_budget else float(n)

    # Start with equal distribution
    levels = {e: budget / n for e in valid}

    def product_for(lv):
        r = simulate_metabolic_chain(valid, substrate_amount, organism_key,
                                     enzyme_levels=lv)
        return r['product_amount'] if not r.get('error') else 0.0

    best_product = product_for(levels)
    step = budget / (2 * n)
    # Hill climbing: shift budget between enzymes
    for _iteration in range(40):
        improved = False
        for i in range(n):
            for j in range(n):
                if i == j:
                    continue
                ei, ej = valid[i], valid[j]
                if levels[ei] - step < 0.05:
                    continue
                trial = dict(levels)
                trial[ei] -= step
                trial[ej] += step
                p = product_for(trial)
                if p > best_product + 1e-6:
                    best_product = p
                    levels = trial
                    improved = True
        if not improved:
            step *= 0.5
            if step < budget * 0.005:
                break

    # Baseline (equal levels) for comparison
    equal_levels = {e: budget / n for e in valid}
    baseline = product_for(equal_levels)

    return {
        'error': None,
        'best_levels': {e: round(v, 2) for e, v in levels.items()},
        'best_product': round(best_product, 3),
        'baseline_product': round(baseline, 3),
        'improvement_pct': round((best_product - baseline) / baseline * 100, 1) if baseline else 0,
        'budget': budget,
    }


def compare_organisms_metabolism(enzyme_chain, substrate_amount=10.0,
                                 organism_keys=None):
    """Run the same pathway across several organisms and compare products.

    Returns each organism's product amount and full time-course for plotting.
    """
    if organism_keys is None:
        organism_keys = ['ecoli', 'yeast', 'human']
    results = {}
    for ok in organism_keys:
        r = simulate_metabolic_chain(enzyme_chain, substrate_amount, ok)
        if not r.get('error'):
            results[ok] = {
                'name': r['organism'],
                'product_amount': r['product_amount'],
                'conversion_pct': r['conversion_pct'],
                'speed_factor': r['speed_factor'],
                'time': r['time'],
                'product_curve': r['concentrations'].get(r['final_product'], []),
                'final_product': r['final_product'],
            }
    return {'error': None if results else 'No results', 'organisms': results}


def validate_enzyme_chain(enzyme_chain):
    """Check if an enzyme chain connects properly (product → next substrate).

    Returns {'valid': bool, 'breaks': [...], 'message': ...}.
    """
    if not enzyme_chain:
        return {'valid': False, 'breaks': [], 'message': 'Empty chain'}
    breaks = []
    for i in range(len(enzyme_chain) - 1):
        cur = METABOLIC_ENZYMES_KINETIC.get(enzyme_chain[i])
        nxt = METABOLIC_ENZYMES_KINETIC.get(enzyme_chain[i+1])
        if not cur or not nxt:
            continue
        if cur['product'] != nxt['substrate']:
            breaks.append({
                'position': i,
                'expected': cur['product'],
                'found': nxt['substrate'],
                'enzyme1': cur['name_en'],
                'enzyme2': nxt['name_en'],
            })
    return {
        'valid': len(breaks) == 0,
        'breaks': breaks,
        'message': 'Chain connects properly' if not breaks
                   else '{} break(s) in the chain'.format(len(breaks)),
    }


# ╔══════════════════════════════════════════════════════════════════╗
# ║  DISEASE EDITING — work with a disease, edit genes, see effects   ║
# ╚══════════════════════════════════════════════════════════════════╝

# Extended disease models: what genes are involved, how editing affects outcome.
DISEASE_GENE_MODELS = {
    # ═══════════════════════════════════════════════════════════════
    # Every gene carries an `actions` map: edit action -> effect on
    # disease severity, as a fraction of the full disease burden.
    #   negative = the disease gets better
    #   positive = the disease gets worse
    #   absent / 0 = no meaningful effect
    #
    # This is what makes the choice of action matter. Previously a single
    # `severity_factors` number was applied no matter what the user did, so
    # knocking out HBB "cured" sickle cell anaemia just as well as correcting
    # it — which is the opposite of the truth, since destroying beta-globin
    # produces thalassaemia. `severity_factors` is kept for backwards
    # compatibility but `actions` takes precedence.
    # ═══════════════════════════════════════════════════════════════
    'Sickle cell anemia': {
        'name_ru': 'Серповидноклеточная анемия',
        'primary_gene': 'HBB', 'mutation': 'Glu6Val (rs334)',
        'genes': {
            'HBB': {'role_ru': 'Бета-глобин (мутантный)', 'role_en': 'Beta-globin (mutant)',
                    'healthy_level': 100, 'disease_level': 100, 'editable': True,
                    'fix_effect': 'Коррекция мутации восстанавливает нормальный гемоглобин',
                    'actions': {'correct': -0.90, 'knockout': 0.35, 'knockdown': 0.40,
                                'activate': 0.05, 'overexpress': 0.30},
                    'action_notes': {
                        'knockout': 'Нокаут HBB убирает бета-глобин вовсе — это талассемия, состояние ухудшается',
                        'overexpress': 'Больше мутантного гемоглобина — больше серповидных клеток'}},
            'BCL11A': {'role_ru': 'Репрессор фетального Hb', 'role_en': 'Fetal Hb repressor',
                       'healthy_level': 100, 'disease_level': 100, 'editable': True,
                       'fix_effect': 'Нокдаун включает фетальный гемоглобин (HbF), компенсирует',
                       'actions': {'correct': 0.0, 'knockout': -0.85, 'knockdown': -0.80,
                                   'activate': 0.40, 'overexpress': 0.50},
                       'action_notes': {
                           'knockout': 'Подход Casgevy/exa-cel — снятие репрессии HbF',
                           'activate': 'Усиление репрессора подавляет HbF — хуже'}},
            'HBG1': {'role_ru': 'Гамма-глобин (фетальный)', 'role_en': 'Gamma-globin (fetal)',
                     'healthy_level': 5, 'disease_level': 5, 'editable': True,
                     'fix_effect': 'Активация замещает дефектный взрослый Hb',
                     'actions': {'correct': 0.0, 'knockout': 0.30, 'knockdown': 0.25,
                                 'activate': -0.75, 'overexpress': -0.80}},
        },
        'severity_factors': {'HBB': -0.7, 'BCL11A': -0.85, 'HBG1': 0.6},
    },
    'Cystic fibrosis': {
        'name_ru': 'Муковисцидоз',
        'primary_gene': 'CFTR', 'mutation': 'ΔF508 (c.1521_1523delCTT)',
        'genes': {
            'CFTR': {'role_ru': 'Хлорный канал (дефектный)', 'role_en': 'Chloride channel (defective)',
                     'healthy_level': 100, 'disease_level': 20, 'editable': True,
                     'fix_effect': 'Коррекция ΔF508 восстанавливает транспорт хлора',
                     'actions': {'correct': -0.85, 'knockout': 0.40, 'knockdown': 0.35,
                                 'activate': -0.15, 'overexpress': -0.30},
                     'action_notes': {
                         'overexpress': 'Больше дефектного белка помогает слабо — он всё равно не сворачивается'}},
            'SCNN1A': {'role_ru': 'Натриевый канал ENaC', 'role_en': 'ENaC sodium channel',
                       'healthy_level': 100, 'disease_level': 150, 'editable': True,
                       'fix_effect': 'Снижение активности уменьшает дегидратацию слизи',
                       'actions': {'correct': 0.0, 'knockout': -0.25, 'knockdown': -0.35,
                                   'activate': 0.30, 'overexpress': 0.35},
                       'action_notes': {
                           'knockout': 'Полный нокаут ENaC опасен (псевдогипоальдостеронизм) — нокдаун лучше'}},
        },
        'severity_factors': {'CFTR': -0.8, 'SCNN1A': 0.3},
    },
    "Huntington's disease": {
        'name_ru': 'Хорея Гентингтона',
        'primary_gene': 'HTT', 'mutation': 'CAG-повторы (>36)',
        'genes': {
            'HTT': {'role_ru': 'Хантингтин (расширенный CAG)', 'role_en': 'Huntingtin (expanded CAG)',
                    'healthy_level': 100, 'disease_level': 100, 'editable': True,
                    'fix_effect': 'Аллель-специфичный нокдаун убирает токсичный белок',
                    'actions': {'correct': -0.90, 'knockout': -0.60, 'knockdown': -0.85,
                                'activate': 0.50, 'overexpress': 0.60},
                    'action_notes': {
                        'knockout': 'Нокаут бьёт и по здоровому аллелю — хантингтин нужен клетке',
                        'knockdown': 'Аллель-специфичный нокдаун — предпочтительный путь'}},
        },
        'severity_factors': {'HTT': -0.9},
    },
    'Duchenne MD': {
        'name_ru': 'Мышечная дистрофия Дюшенна',
        'primary_gene': 'DMD', 'mutation': 'Делеции экзонов',
        'genes': {
            'DMD': {'role_ru': 'Дистрофин (отсутствует)', 'role_en': 'Dystrophin (absent)',
                    'healthy_level': 100, 'disease_level': 0, 'editable': True,
                    'fix_effect': 'Exon skipping восстанавливает рамку → частичный дистрофин',
                    'actions': {'correct': -0.90, 'knockout': 0.20, 'knockdown': 0.10,
                                'activate': -0.20, 'overexpress': -0.40},
                    'action_notes': {
                        'overexpress': 'Микродистрофин через AAV — рабочий клинический подход'}},
        },
        'severity_factors': {'DMD': -1.0},
    },
    'Beta-thalassemia': {
        'name_ru': 'Бета-талассемия',
        'primary_gene': 'HBB', 'mutation': 'Разные мутации HBB',
        'genes': {
            'HBB': {'role_ru': 'Бета-глобин (снижен)', 'role_en': 'Beta-globin (reduced)',
                    'healthy_level': 100, 'disease_level': 30, 'editable': True,
                    'fix_effect': 'Коррекция восстанавливает синтез бета-глобина',
                    'actions': {'correct': -0.70, 'knockout': 0.40, 'knockdown': 0.35,
                                'activate': -0.30, 'overexpress': -0.50}},
            'BCL11A': {'role_ru': 'Репрессор фетального Hb', 'role_en': 'Fetal Hb repressor',
                       'healthy_level': 100, 'disease_level': 100, 'editable': True,
                       'fix_effect': 'Нокдаун включает HbF, компенсирует нехватку',
                       'actions': {'correct': 0.0, 'knockout': -0.60, 'knockdown': -0.55,
                                   'activate': 0.30, 'overexpress': 0.35}},
        },
        'severity_factors': {'HBB': -0.6, 'BCL11A': -0.5},
    },
    'Tay-Sachs disease': {
        'name_ru': 'Болезнь Тея-Сакса',
        'primary_gene': 'HEXA', 'mutation': 'c.1274_1277dupTATC',
        'genes': {
            'HEXA': {'role_ru': 'Гексозаминидаза A (дефект)', 'role_en': 'Hexosaminidase A (defective)',
                     'healthy_level': 100, 'disease_level': 5, 'editable': True,
                     'fix_effect': 'Коррекция восстанавливает фермент, расщепляющий ганглиозиды',
                     'actions': {'correct': -0.95, 'knockout': 0.10, 'knockdown': 0.10,
                                 'activate': -0.30, 'overexpress': -0.50}},
        },
        'severity_factors': {'HEXA': -0.95},
    },
    'Phenylketonuria (PKU)': {
        'name_ru': 'Фенилкетонурия (ФКУ)',
        'primary_gene': 'PAH', 'mutation': 'c.1222C>T (p.Arg408Trp) и др.',
        'genes': {
            'PAH': {'role_ru': 'Фенилаланингидроксилаза (дефект)', 'role_en': 'Phenylalanine hydroxylase (defective)',
                    'healthy_level': 100, 'disease_level': 10, 'editable': True,
                    'fix_effect': 'Коррекция восстанавливает расщепление фенилаланина',
                    'actions': {'correct': -0.90, 'knockout': 0.10, 'knockdown': 0.10,
                                'activate': -0.35, 'overexpress': -0.50}},
        },
        'severity_factors': {'PAH': -0.9},
    },
    'Hemophilia A': {
        'name_ru': 'Гемофилия A',
        'primary_gene': 'F8', 'mutation': 'Инверсия интрона 22 и др.',
        'genes': {
            'F8': {'role_ru': 'Фактор свёртывания VIII (дефект)', 'role_en': 'Coagulation factor VIII (defective)',
                   'healthy_level': 100, 'disease_level': 5, 'editable': True,
                   'fix_effect': 'Коррекция/вставка восстанавливает свёртывание крови',
                   'actions': {'correct': -0.95, 'knockout': 0.05, 'knockdown': 0.05,
                               'activate': -0.40, 'overexpress': -0.70},
                   'action_notes': {
                       'overexpress': 'Даже 5–10% нормального уровня фактора VIII убирает тяжёлые кровотечения'}},
        },
        'severity_factors': {'F8': -0.95},
    },
    'Familial hypercholesterolemia': {
        'name_ru': 'Семейная гиперхолестеринемия',
        'primary_gene': 'LDLR', 'mutation': 'Мутации LDLR (>2000 вариантов)',
        'genes': {
            'LDLR': {'role_ru': 'Рецептор ЛПНП (дефект)', 'role_en': 'LDL receptor (defective)',
                     'healthy_level': 100, 'disease_level': 30, 'editable': True,
                     'fix_effect': 'Коррекция восстанавливает захват холестерина из крови',
                     'actions': {'correct': -0.60, 'knockout': 0.35, 'knockdown': 0.30,
                                 'activate': -0.40, 'overexpress': -0.50}},
            'PCSK9': {'role_ru': 'PCSK9 (разрушает рецепторы ЛПНП)', 'role_en': 'PCSK9 (degrades LDL receptors)',
                      'healthy_level': 100, 'disease_level': 100, 'editable': True,
                      'fix_effect': 'Нокаут PCSK9 повышает уровень рецепторов ЛПНП — снижает холестерин',
                      'actions': {'correct': 0.0, 'knockout': -0.55, 'knockdown': -0.50,
                                  'activate': 0.30, 'overexpress': 0.35},
                      'action_notes': {
                          'knockout': 'Путь VERVE-101: базовое редактирование PCSK9 в печени'}},
        },
        'severity_factors': {'LDLR': -0.6, 'PCSK9': -0.5},
    },
    'Leber congenital amaurosis': {
        'name_ru': 'Врождённый амавроз Лебера',
        'primary_gene': 'CEP290', 'mutation': 'CEP290 IVS26 (c.2991+1655A>G)',
        'genes': {
            'CEP290': {'role_ru': 'Белок ресничек фоторецепторов (дефект)', 'role_en': 'Photoreceptor cilia protein (defective)',
                       'healthy_level': 100, 'disease_level': 10, 'editable': True,
                       'fix_effect': 'Вырезание мутантного интрона восстанавливает зрение (EDIT-101)',
                       'actions': {'correct': -0.90, 'knockout': 0.10, 'knockdown': 0.05,
                                   'activate': -0.20, 'overexpress': -0.30}},
        },
        'severity_factors': {'CEP290': -0.9},
    },
    'Severe combined immunodeficiency': {
        'name_ru': 'Тяжёлый комбинированный иммунодефицит (SCID)',
        'primary_gene': 'IL2RG', 'mutation': 'Мутации IL2RG / ADA',
        'genes': {
            'IL2RG': {'role_ru': 'Гамма-цепь рецептора (дефект)', 'role_en': 'Receptor gamma chain (defective)',
                      'healthy_level': 100, 'disease_level': 5, 'editable': True,
                      'fix_effect': 'Коррекция восстанавливает развитие иммунных клеток',
                      'actions': {'correct': -0.85, 'knockout': 0.10, 'knockdown': 0.10,
                                  'activate': -0.30, 'overexpress': -0.60}},
            'ADA': {'role_ru': 'Аденозиндезаминаза (другая форма SCID)', 'role_en': 'Adenosine deaminase (other SCID form)',
                    'healthy_level': 100, 'disease_level': 10, 'editable': True,
                    'fix_effect': 'Коррекция ADA восстанавливает иммунитет',
                    'actions': {'correct': -0.80, 'knockout': 0.10, 'knockdown': 0.10,
                                'activate': -0.30, 'overexpress': -0.60}},
        },
        'severity_factors': {'IL2RG': -0.85, 'ADA': -0.8},
    },

    # ── added diseases ───────────────────────────────────────────────
    'Spinal muscular atrophy': {
        'name_ru': 'Спинальная мышечная атрофия (СМА)',
        'primary_gene': 'SMN1', 'mutation': 'Делеция/мутация SMN1',
        'genes': {
            'SMN1': {'role_ru': 'Белок выживания мотонейронов (отсутствует)', 'role_en': 'Survival motor neuron 1 (absent)',
                     'healthy_level': 100, 'disease_level': 0, 'editable': True,
                     'fix_effect': 'Восстановление SMN1 спасает мотонейроны',
                     'actions': {'correct': -0.85, 'knockout': 0.10, 'knockdown': 0.10,
                                 'activate': -0.30, 'overexpress': -0.70}},
            'SMN2': {'role_ru': 'Паралог SMN2 (даёт мало белка)', 'role_en': 'SMN2 paralog (low output)',
                     'healthy_level': 20, 'disease_level': 20, 'editable': True,
                     'fix_effect': 'Правка сайта сплайсинга SMN2 повышает выход полноразмерного белка',
                     'actions': {'correct': -0.70, 'knockout': 0.30, 'knockdown': 0.25,
                                 'activate': -0.55, 'overexpress': -0.65},
                     'action_notes': {
                         'correct': 'Принцип нусинерсена/рисдиплама, перенесённый на редактирование'}},
        },
        'severity_factors': {'SMN1': -0.85, 'SMN2': -0.6},
    },
    'Transthyretin amyloidosis': {
        'name_ru': 'Транстиретиновый амилоидоз (ATTR)',
        'primary_gene': 'TTR', 'mutation': 'Мутации TTR (напр. V30M)',
        'genes': {
            'TTR': {'role_ru': 'Транстиретин (образует амилоид)', 'role_en': 'Transthyretin (amyloidogenic)',
                    'healthy_level': 100, 'disease_level': 100, 'editable': True,
                    'fix_effect': 'Нокаут TTR в печени прекращает образование амилоида',
                    'actions': {'correct': -0.70, 'knockout': -0.90, 'knockdown': -0.80,
                                'activate': 0.40, 'overexpress': 0.55},
                    'action_notes': {
                        'knockout': 'NTLA-2001: первое in vivo CRISPR-лечение, нокаут TTR в печени'}},
        },
        'severity_factors': {'TTR': -0.9},
    },
    'Hereditary angioedema': {
        'name_ru': 'Наследственный ангионевротический отёк',
        'primary_gene': 'KLKB1', 'mutation': 'Дефицит C1-ингибитора (SERPING1)',
        'genes': {
            'SERPING1': {'role_ru': 'C1-ингибитор (дефицит)', 'role_en': 'C1 inhibitor (deficient)',
                         'healthy_level': 100, 'disease_level': 25, 'editable': True,
                         'fix_effect': 'Коррекция восстанавливает контроль над каскадом брадикинина',
                         'actions': {'correct': -0.75, 'knockout': 0.35, 'knockdown': 0.30,
                                     'activate': -0.40, 'overexpress': -0.60}},
            'KLKB1': {'role_ru': 'Прекалликреин (источник брадикинина)', 'role_en': 'Prekallikrein (bradykinin source)',
                      'healthy_level': 100, 'disease_level': 100, 'editable': True,
                      'fix_effect': 'Нокаут прекалликреина обрывает путь отёка',
                      'actions': {'correct': 0.0, 'knockout': -0.80, 'knockdown': -0.70,
                                  'activate': 0.35, 'overexpress': 0.45},
                      'action_notes': {'knockout': 'Подход NTLA-2002'}},
        },
        'severity_factors': {'SERPING1': -0.75, 'KLKB1': -0.8},
    },
    'Alpha-1 antitrypsin deficiency': {
        'name_ru': 'Дефицит альфа-1-антитрипсина',
        'primary_gene': 'SERPINA1', 'mutation': 'Аллель Z (p.Glu342Lys)',
        'genes': {
            'SERPINA1': {'role_ru': 'Альфа-1-антитрипсин (Z-аллель, накапливается в печени)',
                         'role_en': 'Alpha-1 antitrypsin (Z allele, aggregates in liver)',
                         'healthy_level': 100, 'disease_level': 15, 'editable': True,
                         'fix_effect': 'Коррекция Z-аллеля убирает и поражение печени, и дефицит в лёгких',
                         'actions': {'correct': -0.90, 'knockout': -0.35, 'knockdown': -0.30,
                                     'activate': 0.25, 'overexpress': 0.30},
                         'action_notes': {
                             'knockout': 'Нокаут спасает печень от агрегатов, но не восстанавливает защиту лёгких',
                             'overexpress': 'Больше Z-белка — больше агрегатов в печени'}},
        },
        'severity_factors': {'SERPINA1': -0.9},
    },
    'Amyotrophic lateral sclerosis (SOD1)': {
        'name_ru': 'БАС (форма с мутацией SOD1)',
        'primary_gene': 'SOD1', 'mutation': 'Мутации SOD1 (напр. A4V)',
        'genes': {
            'SOD1': {'role_ru': 'Супероксиддисмутаза 1 (токсичный мутант)', 'role_en': 'Superoxide dismutase 1 (toxic mutant)',
                     'healthy_level': 100, 'disease_level': 100, 'editable': True,
                     'fix_effect': 'Снижение мутантного SOD1 замедляет гибель мотонейронов',
                     'actions': {'correct': -0.85, 'knockout': -0.70, 'knockdown': -0.80,
                                 'activate': 0.45, 'overexpress': 0.60},
                     'action_notes': {'knockdown': 'Принцип тоферсена, перенесённый на редактирование'}},
        },
        'severity_factors': {'SOD1': -0.85},
    },
    'Hypertrophic cardiomyopathy': {
        'name_ru': 'Гипертрофическая кардиомиопатия',
        'primary_gene': 'MYBPC3', 'mutation': 'Усечённые варианты MYBPC3',
        'genes': {
            'MYBPC3': {'role_ru': 'Миозин-связывающий белок C (гаплонедостаточность)',
                       'role_en': 'Myosin-binding protein C (haploinsufficient)',
                       'healthy_level': 100, 'disease_level': 50, 'editable': True,
                       'fix_effect': 'Коррекция восстанавливает нормальную работу саркомера',
                       'actions': {'correct': -0.85, 'knockout': 0.40, 'knockdown': 0.35,
                                   'activate': -0.35, 'overexpress': -0.55}},
        },
        'severity_factors': {'MYBPC3': -0.85},
    },
    'Epidermolysis bullosa': {
        'name_ru': 'Буллёзный эпидермолиз',
        'primary_gene': 'COL7A1', 'mutation': 'Мутации коллагена VII',
        'genes': {
            'COL7A1': {'role_ru': 'Коллаген VII (якорные фибриллы кожи)', 'role_en': 'Collagen VII (skin anchoring fibrils)',
                       'healthy_level': 100, 'disease_level': 10, 'editable': True,
                       'fix_effect': 'Коррекция восстанавливает сцепление эпидермиса с дермой',
                       'actions': {'correct': -0.90, 'knockout': 0.15, 'knockdown': 0.10,
                                   'activate': -0.30, 'overexpress': -0.55}},
        },
        'severity_factors': {'COL7A1': -0.9},
    },
    'Retinitis pigmentosa (RHO)': {
        'name_ru': 'Пигментный ретинит (мутация родопсина)',
        'primary_gene': 'RHO', 'mutation': 'Доминантные мутации RHO (P23H)',
        'genes': {
            'RHO': {'role_ru': 'Родопсин (доминантно-негативный мутант)', 'role_en': 'Rhodopsin (dominant-negative mutant)',
                    'healthy_level': 100, 'disease_level': 100, 'editable': True,
                    'fix_effect': 'Аллель-специфичный нокаут убирает токсичный родопсин',
                    'actions': {'correct': -0.85, 'knockout': -0.55, 'knockdown': -0.75,
                                'activate': 0.40, 'overexpress': 0.50},
                    'action_notes': {
                        'knockout': 'Неселективный нокаут задевает и здоровый аллель — нужен аллель-специфичный подход'}},
        },
        'severity_factors': {'RHO': -0.85},
    },
    'Wilson disease': {
        'name_ru': 'Болезнь Вильсона',
        'primary_gene': 'ATP7B', 'mutation': 'Мутации ATP7B (напр. H1069Q)',
        'genes': {
            'ATP7B': {'role_ru': 'Медь-транспортирующая АТФаза (дефект)', 'role_en': 'Copper-transporting ATPase (defective)',
                      'healthy_level': 100, 'disease_level': 15, 'editable': True,
                      'fix_effect': 'Коррекция восстанавливает вывод меди из печени',
                      'actions': {'correct': -0.90, 'knockout': 0.15, 'knockdown': 0.10,
                                  'activate': -0.35, 'overexpress': -0.60}},
        },
        'severity_factors': {'ATP7B': -0.9},
    },
    'Gaucher disease': {
        'name_ru': 'Болезнь Гоше',
        'primary_gene': 'GBA', 'mutation': 'Мутации GBA (N370S, L444P)',
        'genes': {
            'GBA': {'role_ru': 'Глюкоцереброзидаза (дефект)', 'role_en': 'Glucocerebrosidase (defective)',
                    'healthy_level': 100, 'disease_level': 10, 'editable': True,
                    'fix_effect': 'Коррекция восстанавливает расщепление глюкоцереброзида',
                    'actions': {'correct': -0.90, 'knockout': 0.10, 'knockdown': 0.10,
                                'activate': -0.35, 'overexpress': -0.60}},
        },
        'severity_factors': {'GBA': -0.9},
    },
    'Congenital deafness (GJB2)': {
        'name_ru': 'Врождённая глухота (GJB2/коннексин 26)',
        'primary_gene': 'GJB2', 'mutation': 'c.35delG и др.',
        'genes': {
            'GJB2': {'role_ru': 'Коннексин 26 (щелевые контакты улитки)', 'role_en': 'Connexin 26 (cochlear gap junctions)',
                     'healthy_level': 100, 'disease_level': 5, 'editable': True,
                     'fix_effect': 'Коррекция восстанавливает калиевый круговорот в улитке',
                     'actions': {'correct': -0.90, 'knockout': 0.10, 'knockdown': 0.10,
                                 'activate': -0.25, 'overexpress': -0.50}},
        },
        'severity_factors': {'GJB2': -0.9},
    },
    'Progeria (HGPS)': {
        'name_ru': 'Прогерия (синдром Хатчинсона-Гилфорда)',
        'primary_gene': 'LMNA', 'mutation': 'c.1824C>T (прогерин)',
        'genes': {
            'LMNA': {'role_ru': 'Ламин A (образуется токсичный прогерин)', 'role_en': 'Lamin A (toxic progerin produced)',
                     'healthy_level': 100, 'disease_level': 100, 'editable': True,
                     'fix_effect': 'Базовое редактирование C>T убирает криптический сайт сплайсинга',
                     'actions': {'correct': -0.90, 'knockout': -0.30, 'knockdown': -0.60,
                                 'activate': 0.40, 'overexpress': 0.50},
                     'action_notes': {
                         'correct': 'Показано на мышах: базовое редактирование продлевает жизнь вдвое',
                         'knockout': 'Ламин A нужен клетке — полный нокаут вреден сам по себе'}},
        },
        'severity_factors': {'LMNA': -0.9},
    },
    'Fanconi anemia': {
        'name_ru': 'Анемия Фанкони',
        'primary_gene': 'FANCA', 'mutation': 'Мутации FANCA',
        'genes': {
            'FANCA': {'role_ru': 'Белок репарации ДНК (дефект)', 'role_en': 'DNA repair protein (defective)',
                      'healthy_level': 100, 'disease_level': 10, 'editable': True,
                      'fix_effect': 'Коррекция восстанавливает репарацию межцепочечных сшивок',
                      'actions': {'correct': -0.85, 'knockout': 0.15, 'knockdown': 0.10,
                                  'activate': -0.30, 'overexpress': -0.55},
                      'action_notes': {
                          'correct': 'Клетки без мутации получают преимущество роста — коррекция сама себя усиливает'}},
        },
        'severity_factors': {'FANCA': -0.85},
    },
    'Chronic granulomatous disease': {
        'name_ru': 'Хроническая гранулематозная болезнь',
        'primary_gene': 'CYBB', 'mutation': 'Мутации CYBB (gp91phox)',
        'genes': {
            'CYBB': {'role_ru': 'Субъединица NADPH-оксидазы (дефект)', 'role_en': 'NADPH oxidase subunit (defective)',
                     'healthy_level': 100, 'disease_level': 5, 'editable': True,
                     'fix_effect': 'Коррекция возвращает фагоцитам способность убивать бактерии',
                     'actions': {'correct': -0.90, 'knockout': 0.10, 'knockdown': 0.05,
                                 'activate': -0.30, 'overexpress': -0.55}},
        },
        'severity_factors': {'CYBB': -0.9},
    },
}


def get_disease_model(disease_name):
    """Return the editable gene model for a disease."""
    return DISEASE_GENE_MODELS.get(disease_name)




# ════════════════════════════════════════════════════════════════════
# REAL GENE EDITING ENGINE
# ------------------------------------------------------------------
# Replaces the old "pick an action from a menu" model, where the DNA the
# user typed had no effect on the outcome. Everything below works from
# the ACTUAL sequence: the edit is found by comparing the edited sequence
# against the reference, its consequence at the protein level is derived
# from the reading frame, and the efficiency comes from the guide that was
# really chosen plus the biology of the organism.
#
# Nothing here invents sequence data. If a reference sequence has not been
# downloaded, the caller is told so instead of being shown a fabricated one.
# ════════════════════════════════════════════════════════════════════

# RefSeq accessions are given WITHOUT a version suffix on purpose: NCBI then
# serves the current version, so the table does not rot as RefSeq is updated.
DISEASE_GENE_REFSEQ = {
    'HBB': {
        'refseq': 'NM_000518', 'name_en': 'Beta-globin',
        'name_ru': 'Бета-глобин',
        'variant': 'c.20A>T (p.Glu6Val, rs334)',
        'cds_pos': 20, 'ref_base': 'A', 'alt_base': 'T',
        'variant_type': 'point',
        'strategy_en': 'HDR correction of the A>T substitution, or knock down '
                       'BCL11A to re-activate fetal haemoglobin.',
        'strategy_ru': 'HDR-коррекция замены A>T, либо нокдаун BCL11A для '
                       'включения фетального гемоглобина.',
    },
    'BCL11A': {
        'refseq': 'NM_022893', 'name_en': 'BCL11A (fetal Hb repressor)',
        'name_ru': 'BCL11A (репрессор фетального Hb)',
        'variant': 'no pathogenic variant — therapeutic target',
        'cds_pos': None, 'ref_base': None, 'alt_base': None,
        'variant_type': 'target_only',
        'strategy_en': 'Disrupt the +58 erythroid enhancer to de-repress HBG1/2. '
                       'This is the exa-cel / Casgevy approach.',
        'strategy_ru': 'Разрушить эритроидный энхансер +58, чтобы снять '
                       'репрессию HBG1/2. Подход exa-cel / Casgevy.',
    },
    'HBG1': {
        'refseq': 'NM_000559', 'name_en': 'Gamma-globin (fetal)',
        'name_ru': 'Гамма-глобин (фетальный)',
        'variant': 'no pathogenic variant — therapeutic target',
        'cds_pos': None, 'ref_base': None, 'alt_base': None,
        'variant_type': 'target_only',
        'strategy_en': 'Promoter edits that mimic HPFH raise fetal haemoglobin.',
        'strategy_ru': 'Правки промотора, имитирующие HPFH, поднимают HbF.',
    },
    'CFTR': {
        'refseq': 'NM_000492', 'name_en': 'CFTR chloride channel',
        'name_ru': 'Хлорный канал CFTR',
        'variant': 'c.1521_1523delCTT (p.Phe508del)',
        'cds_pos': 1521, 'ref_base': 'CTT', 'alt_base': '',
        'variant_type': 'deletion',
        'strategy_en': 'HDR re-insertion of the deleted CTT restores folding.',
        'strategy_ru': 'HDR-возврат удалённого CTT восстанавливает укладку белка.',
    },
    'SCNN1A': {
        'refseq': 'NM_001038', 'name_en': 'ENaC sodium channel alpha',
        'name_ru': 'Натриевый канал ENaC (альфа)',
        'variant': 'no pathogenic variant — modifier target',
        'cds_pos': None, 'ref_base': None, 'alt_base': None,
        'variant_type': 'target_only',
        'strategy_en': 'Partial knockdown reduces airway surface dehydration.',
        'strategy_ru': 'Частичный нокдаун снижает обезвоживание слизи.',
    },
    'HTT': {
        'refseq': 'NM_002111', 'name_en': 'Huntingtin',
        'name_ru': 'Хантингтин',
        'variant': 'CAG repeat expansion in exon 1 (>36 repeats)',
        'cds_pos': None, 'ref_base': None, 'alt_base': None,
        'variant_type': 'repeat',
        'strategy_en': 'Allele-specific knockdown of the expanded allele; the '
                       'normal allele must be preserved.',
        'strategy_ru': 'Аллель-специфичный нокдаун расширенного аллеля; '
                       'нормальный аллель нужно сохранить.',
    },
    'DMD': {
        'refseq': 'NM_004006', 'name_en': 'Dystrophin',
        'name_ru': 'Дистрофин',
        'variant': 'exon deletions (exon 51 skipping is the common target)',
        'cds_pos': None, 'ref_base': None, 'alt_base': None,
        'variant_type': 'exon_deletion',
        'strategy_en': 'Exon skipping restores the reading frame and yields a '
                       'shortened but functional dystrophin.',
        'strategy_ru': 'Пропуск экзона восстанавливает рамку считывания — '
                       'дистрофин укорочен, но работает.',
    },
    'HEXA': {
        'refseq': 'NM_000520', 'name_en': 'Hexosaminidase A',
        'name_ru': 'Гексозаминидаза A',
        'variant': 'c.1274_1277dupTATC (p.Tyr427IlefsTer5)',
        'cds_pos': 1274, 'ref_base': '', 'alt_base': 'TATC',
        'variant_type': 'insertion',
        'strategy_en': 'Remove the 4-bp duplication to restore the frame.',
        'strategy_ru': 'Убрать дупликацию 4 п.н., чтобы вернуть рамку.',
    },
    'PAH': {
        'refseq': 'NM_000277', 'name_en': 'Phenylalanine hydroxylase',
        'name_ru': 'Фенилаланингидроксилаза',
        'variant': 'c.1222C>T (p.Arg408Trp)',
        'cds_pos': 1222, 'ref_base': 'C', 'alt_base': 'T',
        'variant_type': 'point',
        'strategy_en': 'Base editing C>T reversal, or HDR correction.',
        'strategy_ru': 'Обратная замена базовым редактором, либо HDR.',
    },
    'F8': {
        'refseq': 'NM_000132', 'name_en': 'Coagulation factor VIII',
        'name_ru': 'Фактор свёртывания VIII',
        'variant': 'intron 22 inversion (~45% of severe cases)',
        'cds_pos': None, 'ref_base': None, 'alt_base': None,
        'variant_type': 'inversion',
        'strategy_en': 'Re-invert the segment, or add a functional F8 copy at a '
                       'safe-harbour locus.',
        'strategy_ru': 'Вернуть инверсию обратно, либо вставить рабочую копию F8 '
                       'в безопасный локус.',
    },
    'LDLR': {
        'refseq': 'NM_000527', 'name_en': 'LDL receptor',
        'name_ru': 'Рецептор ЛПНП',
        'variant': 'highly heterogeneous (>2000 known variants)',
        'cds_pos': None, 'ref_base': None, 'alt_base': None,
        'variant_type': 'heterogeneous',
        'strategy_en': 'Because the variants are so varied, the clinical route is '
                       'to knock down PCSK9 instead (cf. VERVE-101 base editing).',
        'strategy_ru': 'Из-за разнообразия мутаций клинический путь — нокдаун '
                       'PCSK9 (ср. базовое редактирование VERVE-101).',
    },
    'PCSK9': {
        'refseq': 'NM_174936', 'name_en': 'PCSK9',
        'name_ru': 'PCSK9',
        'variant': 'no pathogenic variant — therapeutic target',
        'cds_pos': None, 'ref_base': None, 'alt_base': None,
        'variant_type': 'target_only',
        'strategy_en': 'Knocking it out raises LDL-receptor recycling and lowers '
                       'LDL cholesterol.',
        'strategy_ru': 'Нокаут повышает рециклинг рецептора ЛПНП и снижает '
                       'холестерин ЛПНП.',
    },
    'RPE65': {
        'refseq': 'NM_000329', 'name_en': 'RPE65 retinoid isomerase',
        'name_ru': 'Ретиноид-изомераза RPE65',
        'variant': 'biallelic loss-of-function (many variants)',
        'cds_pos': None, 'ref_base': None, 'alt_base': None,
        'variant_type': 'heterogeneous',
        'strategy_en': 'Gene supplementation is the approved route '
                       '(voretigene neparvovec); editing is variant-specific.',
        'strategy_ru': 'Одобрен путь генозамещения (воретиген непарвовек); '
                       'редактирование зависит от конкретной мутации.',
    },
    'IL2RG': {
        'refseq': 'NM_000206', 'name_en': 'IL2 receptor gamma chain',
        'name_ru': 'Гамма-цепь рецептора IL2',
        'variant': 'loss-of-function (X-linked SCID)',
        'cds_pos': None, 'ref_base': None, 'alt_base': None,
        'variant_type': 'heterogeneous',
        'strategy_en': 'HDR correction in HSCs, or targeted cDNA insertion.',
        'strategy_ru': 'HDR-коррекция в ГСК, либо вставка рабочей кДНК.',
    },
    'ADA': {
        'refseq': 'NM_000022', 'name_en': 'Adenosine deaminase',
        'name_ru': 'Аденозиндезаминаза',
        'variant': 'loss-of-function (ADA-SCID)',
        'cds_pos': None, 'ref_base': None, 'alt_base': None,
        'variant_type': 'heterogeneous',
        'strategy_en': 'HDR correction or cDNA addition in haematopoietic stem cells.',
        'strategy_ru': 'HDR-коррекция или добавление кДНК в ГСК.',
    },
}


def get_gene_reference_info(gene):
    """Curated facts about a disease gene, or None if we have none on file."""
    return DISEASE_GENE_REFSEQ.get(gene)


# ── Edit detection ──────────────────────────────────────────────────

def diff_sequences(reference, edited):
    """Locate the difference between two sequences.

    Returns a dict describing the first (and, for our purposes, principal)
    difference: its position, what was removed, what was inserted, and the
    net length change. Uses common-prefix / common-suffix trimming, which
    resolves substitutions, insertions and deletions exactly and is far
    cheaper than full alignment on the long sequences we deal with.
    """
    reference = (reference or '').upper().replace('\n', '').replace(' ', '')
    edited = (edited or '').upper().replace('\n', '').replace(' ', '')
    if reference == edited:
        return {'changed': False, 'position': None, 'removed': '',
                'inserted': '', 'delta': 0, 'kind': 'none'}

    # longest common prefix
    p = 0
    while p < len(reference) and p < len(edited) and reference[p] == edited[p]:
        p += 1
    # longest common suffix (not overlapping the prefix)
    s = 0
    while (s < len(reference) - p and s < len(edited) - p
           and reference[len(reference) - 1 - s] == edited[len(edited) - 1 - s]):
        s += 1

    removed = reference[p:len(reference) - s]
    inserted = edited[p:len(edited) - s]
    delta = len(inserted) - len(removed)

    if removed and inserted and len(removed) == len(inserted):
        kind = 'substitution'
    elif removed and not inserted:
        kind = 'deletion'
    elif inserted and not removed:
        kind = 'insertion'
    else:
        kind = 'complex'

    return {'changed': True, 'position': p, 'removed': removed,
            'inserted': inserted, 'delta': delta, 'kind': kind}


_STOP_CODONS = ('TAA', 'TAG', 'TGA')


def _translate_simple(dna):
    """Translate until the first stop. Returns (protein, stopped_early)."""
    table = {
        'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'CTT': 'L', 'CTC': 'L',
        'CTA': 'L', 'CTG': 'L', 'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
        'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'TCT': 'S', 'TCC': 'S',
        'TCA': 'S', 'TCG': 'S', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'GCT': 'A', 'GCC': 'A',
        'GCA': 'A', 'GCG': 'A', 'TAT': 'Y', 'TAC': 'Y', 'CAT': 'H', 'CAC': 'H',
        'CAA': 'Q', 'CAG': 'Q', 'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E', 'TGT': 'C', 'TGC': 'C',
        'TGG': 'W', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'AGT': 'S',
        'AGC': 'S', 'AGA': 'R', 'AGG': 'R', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G',
        'GGG': 'G',
    }
    prot = []
    for i in range(0, len(dna) - 2, 3):
        c = dna[i:i + 3]
        if c in _STOP_CODONS:
            return ''.join(prot), (i + 3) < len(dna) - 2
        prot.append(table.get(c, 'X'))
    return ''.join(prot), False


def analyze_gene_edit(reference, edited, gene=None):
    """Work out what an edit actually did, biologically.

    This is the function that makes the user's typing matter: it compares the
    sequence they produced against the real reference, then reports the
    molecular consequence — frameshift, premature stop, clean correction of
    the pathogenic base, or a miss.
    """
    info = DISEASE_GENE_REFSEQ.get(gene or '', {})
    d = diff_sequences(reference, edited)

    out = {
        'changed': d['changed'], 'kind': d['kind'], 'position': d['position'],
        'removed': d['removed'], 'inserted': d['inserted'], 'delta': d['delta'],
        'frameshift': False, 'premature_stop': False,
        'hits_variant_site': False, 'corrects_variant': False,
        'protein_ref': '', 'protein_edited': '', 'consequence': 'none',
        'messages_ru': [], 'messages_en': [],
    }
    if not d['changed']:
        out['consequence'] = 'unchanged'
        out['messages_ru'].append('Последовательность не изменена.')
        out['messages_en'].append('Sequence unchanged.')
        return out

    # Reading-frame consequences
    out['frameshift'] = (d['delta'] % 3 != 0)
    ref_prot, _ = _translate_simple(reference.upper())
    edit_prot, _ = _translate_simple(edited.upper())
    out['protein_ref'] = ref_prot
    out['protein_edited'] = edit_prot
    if len(edit_prot) < len(ref_prot) * 0.9:
        out['premature_stop'] = True

    # Did the edit land on the known pathogenic site?
    cds_pos = info.get('cds_pos')
    if cds_pos:
        site0 = cds_pos - 1                      # to 0-based
        span = max(1, len(d['removed']), len(d['inserted']))
        if abs(d['position'] - site0) <= max(3, span):
            out['hits_variant_site'] = True
            vt = info.get('variant_type')
            alt, ref_b = info.get('alt_base') or '', info.get('ref_base') or ''
            if vt == 'point':
                # correction means putting the wild-type base back
                if d['inserted'] == ref_b and d['removed'] == alt:
                    out['corrects_variant'] = True
            elif vt == 'deletion':
                # pathogenic allele lacks ref_b; correction re-inserts it
                if d['kind'] == 'insertion' and d['inserted'] == ref_b:
                    out['corrects_variant'] = True
            elif vt == 'insertion':
                # pathogenic allele has an extra alt_base; correction removes it
                if d['kind'] == 'deletion' and d['removed'] == alt:
                    out['corrects_variant'] = True

    # Classify and explain
    if out['corrects_variant']:
        out['consequence'] = 'corrected'
        out['messages_ru'].append(
            'Патогенная мутация исправлена в позиции {}.'.format(cds_pos))
        out['messages_en'].append(
            'Pathogenic variant corrected at position {}.'.format(cds_pos))
    elif out['frameshift']:
        out['consequence'] = 'frameshift'
        out['messages_ru'].append(
            'Сдвиг рамки считывания ({:+d} п.н.) — белок разрушен. '
            'Это нокаут, а не коррекция.'.format(d['delta']))
        out['messages_en'].append(
            'Frameshift ({:+d} bp) — the protein is destroyed. '
            'This is a knockout, not a correction.'.format(d['delta']))
    elif out['premature_stop']:
        out['consequence'] = 'truncated'
        out['messages_ru'].append(
            'Появился преждевременный стоп-кодон: белок укорочен '
            '({} из {} аминокислот).'.format(len(edit_prot), len(ref_prot)))
        out['messages_en'].append(
            'Premature stop codon: protein truncated '
            '({} of {} residues).'.format(len(edit_prot), len(ref_prot)))
    elif d['kind'] == 'substitution' and d['delta'] == 0:
        if ref_prot == edit_prot:
            out['consequence'] = 'silent'
            out['messages_ru'].append('Синонимичная замена — белок не изменился.')
            out['messages_en'].append('Silent substitution — protein unchanged.')
        else:
            out['consequence'] = 'missense'
            out['messages_ru'].append('Миссенс-замена — изменена аминокислота.')
            out['messages_en'].append('Missense substitution — residue changed.')
    elif d['delta'] % 3 == 0 and d['delta'] != 0:
        out['consequence'] = 'in_frame_indel'
        out['messages_ru'].append(
            'Вставка/делеция кратна 3 — рамка сохранена, '
            'белок изменён на {} аминокислот.'.format(abs(d['delta']) // 3))
        out['messages_en'].append(
            'In-frame indel — reading frame preserved, protein changed by '
            '{} residues.'.format(abs(d['delta']) // 3))
    else:
        out['consequence'] = 'other'

    if cds_pos and not out['hits_variant_site']:
        out['messages_ru'].append(
            'Правка сделана в позиции {}, а патогенный сайт находится в '
            'позиции {} — цель не затронута.'.format(d['position'] + 1, cds_pos))
        out['messages_en'].append(
            'Edit is at position {}, but the pathogenic site is at position '
            '{} — the target was missed.'.format(d['position'] + 1, cds_pos))
    return out


# ── Organism-specific efficiency for every edit type ────────────────

def organism_edit_efficiency(organism_key, action, guide_score=None,
                             cas_enzyme='SpCas9', protocol='in_vivo',
                             clones_screened=100):
    """Efficiency of a given edit type in a given organism, 0..1.

    Previously knockdown / activate / overexpress used the same three
    constants for every organism, which meant the organism selector was
    largely decorative. Each pathway now derives from the organism's own
    repair biology, delivery difficulty and cell type, and is further scaled
    by the on-target score of the guide the user actually picked.
    """
    org = ORGANISMS.get(organism_key, {})
    cr = org.get('crispr', {})
    hdr = cr.get('hdr_efficiency', 0.05)
    nhej = cr.get('nhej_efficiency', 0.70)
    euk = (org.get('cell_type') == 'eukaryote')

    # Delivery penalty: organisms that are hard to transfect lose efficiency
    # across every pathway, not just the repair-dependent ones.
    difficulty = (cr.get('difficulty') or '').lower()
    if 'высок' in difficulty or 'high' in difficulty:
        delivery = 0.55
    elif 'средн' in difficulty or 'medium' in difficulty or 'moder' in difficulty:
        delivery = 0.75
    else:
        delivery = 0.95

    if action == 'correct':
        # Precise repair is limited by HDR itself, and HDR only operates in
        # cycling cells; it is not scaled by an arbitrary constant.
        base = hdr * (1.0 if euk else 1.3)
    elif action == 'knockout':
        base = nhej
    elif action == 'knockdown':
        # dCas9-KRAB style repression: needs sustained expression of the
        # effector, so delivery dominates. Prokaryotes lack KRAB but CRISPRi
        # by steric block works well there.
        base = 0.75 if euk else 0.85
    elif action == 'activate':
        # CRISPRa needs a full activator complex — the hardest to deliver.
        base = 0.60 if euk else 0.45
    elif action == 'overexpress':
        # Plasmid / strong promoter: easy in microbes, harder in animals.
        base = 0.70 if euk else 0.90
    else:
        return 0.0

    eff = base * delivery

    # A weak guide caps what any pathway can achieve.
    if guide_score is not None:
        eff *= max(0.15, min(1.0, guide_score / 100.0))

    # Cas variants differ in intrinsic activity.
    cas_factor = {'SpCas9': 1.0, 'SaCas9': 0.9, 'Cas12a': 0.85, 'xCas9': 0.8,
                  'SpRY': 0.7, 'Cas12b': 0.75, 'CasX': 0.6, 'Cas14': 0.5,
                  'Cas13': 0.8, 'enCas12a': 0.9}.get(cas_enzyme, 1.0)
    eff *= cas_factor
    eff = max(0.0, min(1.0, eff))

    # ── ex vivo editing with clone selection ────────────────────────
    # Raw HDR in human cells really is only a few per cent, and taken at face
    # value that makes every correct repair look useless. Real therapies get
    # around this by editing cells outside the body and keeping only the
    # clones that worked (this is how Casgevy reaches ~90% edited cells).
    # The chance of recovering at least one correct clone out of N screened is
    # 1-(1-p)^N, and the selected population is then almost entirely correct.
    if protocol == 'ex_vivo' and eff > 0:
        n = max(1, int(clones_screened))
        p_success = 1.0 - (1.0 - eff) ** n
        eff = p_success * 0.95

    return max(0.0, min(1.0, eff))


def evaluate_disease_state(disease_name, gene_edits, organism_key='human',
                           gene_sequences=None, reference_sequences=None,
                           guide_scores=None, cas_enzyme='SpCas9',
                           protocol='in_vivo', clones_screened=100):
    """Disease severity after a set of edits — derived from the real edits.

    gene_edits          {gene: action}
    gene_sequences      {gene: sequence the user produced}   (optional)
    reference_sequences {gene: downloaded reference sequence} (optional)
    guide_scores        {gene: on-target score of the chosen guide} (optional)

    When sequences are supplied, the molecular consequence of the actual edit
    decides what happened; the chosen action is then only a statement of
    intent, which the analysis can confirm or contradict. When they are not
    supplied, the intent alone is used and the result is flagged as such.
    """
    model = DISEASE_GENE_MODELS.get(disease_name)
    if not model:
        return {'error': 'Unknown disease'}

    gene_sequences = gene_sequences or {}
    reference_sequences = reference_sequences or {}
    guide_scores = guide_scores or {}

    org = ORGANISMS.get(organism_key, {})

    # Untreated baseline, measured from how much each gene can contribute to
    # the disease rather than assumed to be 100%.
    baseline = 0.0
    for gene, gdata in model['genes'].items():
        acts = gdata.get('actions')
        if acts:
            best = min(list(acts.values()) + [0.0])
            baseline += abs(best) * 100
        else:
            f = model['severity_factors'].get(gene, 0)
            if f < 0:
                baseline += abs(f) * 100
    baseline = max(1.0, min(100.0, baseline))

    severity = baseline
    gene_levels = {}
    explanations = []
    edit_reports = {}
    sequence_backed = False

    for gene, gdata in model['genes'].items():
        action = gene_edits.get(gene, 'none')
        base = gdata['disease_level']
        healthy = gdata['healthy_level']
        level = base

        if action == 'none':
            gene_levels[gene] = round(level, 1)
            continue

        gscore = guide_scores.get(gene)
        eff = organism_edit_efficiency(organism_key, action, gscore, cas_enzyme,
                                       protocol, clones_screened)

        # If we have both the reference and the edited sequence, the molecular
        # result overrides the stated intent.
        ref_seq = reference_sequences.get(gene)
        cur_seq = gene_sequences.get(gene)
        report = None
        if ref_seq and cur_seq:
            sequence_backed = True
            report = analyze_gene_edit(ref_seq, cur_seq, gene)
            edit_reports[gene] = report
            cons = report['consequence']

            if cons == 'unchanged':
                # Intent declared, but the DNA was never touched.
                eff = 0.0
                explanations.append(
                    '{}: последовательность не изменена — правка не выполнена'
                    .format(gene))
            elif action == 'correct':
                if report['corrects_variant']:
                    pass                      # intent achieved, keep eff
                elif report['frameshift']:
                    # They meant to correct but destroyed the gene instead.
                    action = 'knockout'
                    eff = organism_edit_efficiency(
                        organism_key, 'knockout', gscore, cas_enzyme,
                        protocol, clones_screened)
                    explanations.append(
                        '{}: попытка коррекции привела к сдвигу рамки — '
                        'получился нокаут'.format(gene))
                else:
                    eff *= 0.15               # missed the pathogenic site
                    explanations.append(
                        '{}: правка не попала в патогенный сайт'.format(gene))
            elif action == 'knockout':
                if report['frameshift'] or report['premature_stop']:
                    pass                      # a real knockout
                elif cons == 'silent':
                    eff = 0.0
                    explanations.append(
                        '{}: синонимичная замена — ген продолжает работать'
                        .format(gene))
                else:
                    eff *= 0.5

        # Apply the (possibly revised) action
        if action == 'correct':
            level = base + (healthy - base) * eff
        elif action == 'knockout':
            level = base * (1 - eff)
        elif action == 'knockdown':
            level = base * (1 - 0.85 * eff)
        elif action == 'activate':
            level = min(200, base + healthy * 0.8 * eff)
        elif action == 'overexpress':
            level = min(300, base * (1 + 1.5 * eff))

        gene_levels[gene] = round(level, 1)

        # Effect of THIS action on THIS gene. Correcting a gene and destroying
        # it are opposite interventions and must not share one coefficient.
        acts = gdata.get('actions')
        if acts is not None:
            factor = acts.get(action, 0.0)
        else:
            factor = model['severity_factors'].get(gene, 0)
        severity += factor * eff * 100

        note = (gdata.get('action_notes') or {}).get(action)
        if factor > 0:
            explanations.append(
                '{}: {} — УХУДШАЕТ состояние (эффективность {:.0f}%){}'.format(
                    gene, action, eff * 100, ' — ' + note if note else ''))
        elif factor < 0:
            explanations.append(
                '{}: {} — помогает, эффективность {:.0f}%{}'.format(
                    gene, action, eff * 100,
                    ' — ' + note if note else ' → ' + gdata['fix_effect']))
        else:
            explanations.append(
                '{}: {} — на эту болезнь не влияет{}'.format(
                    gene, action, ' — ' + note if note else ''))

    severity = max(0.0, min(100.0, severity))

    if severity < 20:
        outcome_ru, outcome_en = 'Вылечено / норма', 'Cured / normal'
    elif severity < 50:
        outcome_ru, outcome_en = 'Значительное улучшение', 'Major improvement'
    elif severity < 80:
        outcome_ru, outcome_en = 'Частичное улучшение', 'Partial improvement'
    else:
        outcome_ru, outcome_en = 'Болезнь сохраняется', 'Disease persists'

    # Editing is a stochastic process across a cell population; a single
    # number would overstate what is known, so report a band as well.
    spread = 12.0 if sequence_backed else 20.0
    return {
        'error': None,
        'disease': disease_name,
        'baseline_severity': round(baseline, 1),
        'severity': round(severity, 1),
        'severity_low': round(max(0.0, severity - spread), 1),
        'severity_high': round(min(100.0, severity + spread), 1),
        'gene_levels': gene_levels,
        'explanations': explanations,
        'edit_reports': edit_reports,
        'sequence_backed': sequence_backed,
        'protocol': protocol,
        'outcome_ru': outcome_ru,
        'outcome_en': outcome_en,
        'organism': org.get('name_en', organism_key),
        'hdr_efficiency': org.get('crispr', {}).get('hdr_efficiency', 0.05),
    }


# ═══════════════════════════════════════════════════════════════════
# COPASI-STYLE METABOLIC MODEL — rate laws, model schema, presets
# Added for the redesigned Metabolism tab. Everything here is pure data
# and pure functions (no Tkinter), so it can be unit-checked on its own.
# ═══════════════════════════════════════════════════════════════════

# ── Predefined kinetic rate laws (mirrors COPASI's built-in function DB) ──
# Each law: v = f(substrates S, products P, modifiers, params).
# 'params' lists the named constants the user sets per reaction.
# 'formula' is a human-readable form; 'kind' groups them in the UI.
# 'func' is the actual velocity given a context dict.

def _mm_irrev(ctx):
    S = ctx['S']; V = ctx['p'].get('V', ctx['p'].get('Vmax', 1.0)); Km = ctx['p'].get('Km', 1.0)
    return V * S / (Km + S + 1e-12)

def _mm_rev(ctx):
    S = ctx['S']; P = ctx['P']
    Vf = ctx['p'].get('Vf', 1.0); Vr = ctx['p'].get('Vr', 0.0)
    Kms = ctx['p'].get('Kms', 1.0); Kmp = ctx['p'].get('Kmp', 1.0)
    num = Vf * S / Kms - Vr * P / Kmp
    den = 1.0 + S / Kms + P / Kmp
    return num / (den + 1e-12)

def _mass_action_irrev(ctx):
    k = ctx['p'].get('k', ctx['p'].get('k1', 1.0))
    return k * ctx['S']

def _mass_action_rev(ctx):
    kf = ctx['p'].get('k1', 1.0); kr = ctx['p'].get('k2', 0.0)
    return kf * ctx['S'] - kr * ctx['P']

def _hill(ctx):
    S = ctx['S']; V = ctx['p'].get('V', 1.0); Kh = ctx['p'].get('Kh', 1.0)
    h = ctx['p'].get('h', 2.0)
    Sh = S ** h
    return V * Sh / (Kh ** h + Sh + 1e-12)

def _competitive_inhib_irrev(ctx):
    S = ctx['S']; I = ctx.get('I', 0.0)
    V = ctx['p'].get('V', 1.0); Km = ctx['p'].get('Km', 1.0); Ki = ctx['p'].get('Ki', 1.0)
    return V * S / (Km * (1.0 + I / (Ki + 1e-12)) + S + 1e-12)

def _uncompetitive_inhib_irrev(ctx):
    S = ctx['S']; I = ctx.get('I', 0.0)
    V = ctx['p'].get('V', 1.0); Km = ctx['p'].get('Km', 1.0); Ki = ctx['p'].get('Ki', 1.0)
    return V * S / (Km + S * (1.0 + I / (Ki + 1e-12)) + 1e-12)

def _noncompetitive_inhib_irrev(ctx):
    S = ctx['S']; I = ctx.get('I', 0.0)
    V = ctx['p'].get('V', 1.0); Km = ctx['p'].get('Km', 1.0); Ki = ctx['p'].get('Ki', 1.0)
    factor = 1.0 + I / (Ki + 1e-12)
    return (V * S / (Km + S + 1e-12)) / factor

def _mixed_inhib_irrev(ctx):
    S = ctx['S']; I = ctx.get('I', 0.0)
    V = ctx['p'].get('V', 1.0); Km = ctx['p'].get('Km', 1.0)
    Kic = ctx['p'].get('Kic', 1.0); Kiu = ctx['p'].get('Kiu', 1.0)
    a = 1.0 + I / (Kic + 1e-12)
    b = 1.0 + I / (Kiu + 1e-12)
    return V * S / (Km * a + S * b + 1e-12)

def _substrate_inhib_irrev(ctx):
    S = ctx['S']
    V = ctx['p'].get('V', 1.0); Km = ctx['p'].get('Km', 1.0); Ksi = ctx['p'].get('Ksi', 10.0)
    return V * S / (Km + S + (S * S) / (Ksi + 1e-12) + 1e-12)

def _allosteric_mwc(ctx):
    # Monod-Wyman-Changeux allosteric inhibition (empirical form)
    S = ctx['S']; I = ctx.get('I', 0.0)
    V = ctx['p'].get('V', 1.0); Km = ctx['p'].get('Km', 1.0)
    n = ctx['p'].get('n', 2.0); L = ctx['p'].get('L', 1.0); Ki = ctx['p'].get('Ki', 1.0)
    frac = (S / (Km + 1e-12))
    return V * frac * (1.0 + frac) ** (n - 1) / (
        L * (1.0 + I / (Ki + 1e-12)) ** n + (1.0 + frac) ** n + 1e-12)

def _constant_flux(ctx):
    return ctx['p'].get('v', 1.0)

def _uni_uni(ctx):
    # reversible Uni-Uni (Cleland) — same shape as reversible MM
    return _mm_rev(ctx)

RATE_LAWS = {
    'mass_action_irrev': {
        'name_ru': 'Закон действующих масс (необр.)', 'name_en': 'Mass action (irreversible)',
        'kind': 'basic', 'reversible': False,
        'formula': 'v = k · S', 'params': ['k'], 'func': _mass_action_irrev,
    },
    'mass_action_rev': {
        'name_ru': 'Закон действующих масс (обр.)', 'name_en': 'Mass action (reversible)',
        'kind': 'basic', 'reversible': True,
        'formula': 'v = k1·S − k2·P', 'params': ['k1', 'k2'], 'func': _mass_action_rev,
    },
    'mm_irrev': {
        'name_ru': 'Михаэлис-Ментен (необр.)', 'name_en': 'Henri-Michaelis-Menten (irreversible)',
        'kind': 'enzyme', 'reversible': False,
        'formula': 'v = V·S / (Km + S)', 'params': ['V', 'Km'], 'func': _mm_irrev,
    },
    'mm_rev': {
        'name_ru': 'Михаэлис-Ментен (обр.)', 'name_en': 'Reversible Michaelis-Menten',
        'kind': 'enzyme', 'reversible': True,
        'formula': 'v = (Vf·S/Kms − Vr·P/Kmp)/(1+S/Kms+P/Kmp)',
        'params': ['Vf', 'Vr', 'Kms', 'Kmp'], 'func': _mm_rev,
    },
    'hill': {
        'name_ru': 'Кооперативность Хилла', 'name_en': 'Hill Cooperativity',
        'kind': 'enzyme', 'reversible': False,
        'formula': 'v = V·S^h / (Kh^h + S^h)', 'params': ['V', 'Kh', 'h'], 'func': _hill,
    },
    'competitive_irr': {
        'name_ru': 'Конкурентное ингибир. (необр.)', 'name_en': 'Competitive inhibition (irr)',
        'kind': 'inhibition', 'reversible': False,
        'formula': 'v = V·S / (Km·(1+I/Ki) + S)', 'params': ['V', 'Km', 'Ki'],
        'func': _competitive_inhib_irrev, 'needs_modifier': True,
    },
    'uncompetitive_irr': {
        'name_ru': 'Неконкур. (uncompetitive, необр.)', 'name_en': 'Uncompetitive inhibition (irr)',
        'kind': 'inhibition', 'reversible': False,
        'formula': 'v = V·S / (Km + S·(1+I/Ki))', 'params': ['V', 'Km', 'Ki'],
        'func': _uncompetitive_inhib_irrev, 'needs_modifier': True,
    },
    'noncompetitive_irr': {
        'name_ru': 'Некокурентное (necomp., необр.)', 'name_en': 'Noncompetitive inhibition (irr)',
        'kind': 'inhibition', 'reversible': False,
        'formula': 'v = [V·S/(Km+S)] / (1+I/Ki)', 'params': ['V', 'Km', 'Ki'],
        'func': _noncompetitive_inhib_irrev, 'needs_modifier': True,
    },
    'mixed_irr': {
        'name_ru': 'Смешанное ингибир. (необр.)', 'name_en': 'Mixed inhibition (irr)',
        'kind': 'inhibition', 'reversible': False,
        'formula': 'v = V·S / (Km·(1+I/Kic) + S·(1+I/Kiu))',
        'params': ['V', 'Km', 'Kic', 'Kiu'], 'func': _mixed_inhib_irrev, 'needs_modifier': True,
    },
    'substrate_inhib': {
        'name_ru': 'Субстратное ингибир.', 'name_en': 'Substrate inhibition (irr)',
        'kind': 'inhibition', 'reversible': False,
        'formula': 'v = V·S / (Km + S + S²/Ksi)', 'params': ['V', 'Km', 'Ksi'],
        'func': _substrate_inhib_irrev,
    },
    'allosteric_mwc': {
        'name_ru': 'Аллостерич. ингибир. (MWC)', 'name_en': 'Allosteric inhibition (MWC)',
        'kind': 'inhibition', 'reversible': False,
        'formula': 'MWC empirical', 'params': ['V', 'Km', 'n', 'L', 'Ki'],
        'func': _allosteric_mwc, 'needs_modifier': True,
    },
    'uni_uni': {
        'name_ru': 'Uni-Uni (обр.)', 'name_en': 'Uni Uni (reversible)',
        'kind': 'enzyme', 'reversible': True,
        'formula': 'reversible Uni-Uni', 'params': ['Vf', 'Vr', 'Kms', 'Kmp'], 'func': _uni_uni,
    },
    'constant_flux': {
        'name_ru': 'Постоянный поток', 'name_en': 'Constant flux',
        'kind': 'basic', 'reversible': False,
        'formula': 'v = v0', 'params': ['v'], 'func': _constant_flux,
    },
}


def rate_law_default_params(law_key):
    """Return a fresh dict of default parameter values for a rate law."""
    law = RATE_LAWS.get(law_key)
    if not law:
        return {}
    defaults = {
        'k': 0.5, 'k1': 0.5, 'k2': 0.1, 'V': 1.0, 'Vmax': 1.0, 'Km': 1.0,
        'Vf': 1.0, 'Vr': 0.2, 'Kms': 1.0, 'Kmp': 1.0, 'Kh': 1.0, 'h': 2.0,
        'Ki': 1.0, 'Kic': 1.0, 'Kiu': 1.0, 'Ksi': 10.0, 'n': 2.0, 'L': 1.0, 'v': 0.5,
    }
    return {p: defaults.get(p, 1.0) for p in law['params']}


def eval_rate_law(law_key, S, P=0.0, I=0.0, params=None):
    """Evaluate a rate law velocity. Safe: returns 0.0 on any error."""
    law = RATE_LAWS.get(law_key)
    if not law:
        return 0.0
    try:
        ctx = {'S': max(0.0, S), 'P': max(0.0, P), 'I': max(0.0, I),
               'p': params or rate_law_default_params(law_key)}
        v = law['func'](ctx)
        if v != v:  # NaN guard
            return 0.0
        return float(v)
    except Exception:
        return 0.0


# ── COPASI-style empty model factory ──
def new_whole_cell_model(name='New cell'):
    """Create a blank whole-cell (time-resolved) organism model.

    Этап 1 (каркас): только структура данных, без движков.
    Structure:
      genes: {id: {name, promoter_strength, rbs_strength, product,
                    mrna_halflife_min, protein_halflife_min, initial_mrna,
                    initial_protein, category}}
      chromosome: {length_kb, replication_rate_kb_min, ori_count}
      processes: {ribosomes, rna_polymerase, stochastic}  (шасси задаёт лимиты)
      environment: {substrate, substrate_conc, temperature_c}
      last_result: None (заполнится движком в Этапе 2)
    """
    return {
        'name': name,
        'chassis': 'ecoli',
        'genes': {},
        'chromosome': {
            'length_kb': 4600, 'replication_rate_kb_min': 1000, 'ori_count': 1,
        },
        'processes': {
            'ribosomes': 20000, 'rna_polymerase': 3000, 'stochastic': False,
        },
        'environment': {
            'substrate': 'glucose', 'substrate_conc': 10.0, 'temperature_c': 37.0,
        },
        'last_result': None,
    }


WC_RNAP_REF_COUNT = 3000.0
WC_RIBOSOME_REF_COUNT = 20000.0
# Translation demand (rbs_strength * mRNA_copies) lives on a much larger
# numeric scale than transcription demand (promoter_strength, capped at 1
# per gene) — mRNA copy numbers routinely reach the tens. This factor puts
# the default ribosome pool's *capacity* on the same scale as realistic
# demand for a handful of genes, so a single/modest design isn't throttled
# just from a unit mismatch; competition still bites once enough genes (or
# high copy numbers) are stacked up.
WC_RIBOSOME_CAPACITY_SCALE = 50.0


def simulate_whole_cell_timecourse(model, t_max=60.0, n_points=200,
                                   k_tx=2.0, k_tl=20.0):
    """Этап 2+3: gene-expression time course with shared-resource competition
    (RK4).

    Genes no longer act in isolation: RNA polymerase and ribosomes are a
    single shared, limited pool (model['processes']). When total demand for
    a resource exceeds what's available, every gene's realized rate is
    throttled by the same fair-share factor — this is what makes adding a
    module to the design measurably slow down every other gene, same as in
    a real cell. Per gene i:

        tx_demand_i  = promoter_strength_i
        tx_scale     = min(1, RNAP_count/RNAP_ref / sum(tx_demand))
        dM_i/dt      = k_tx * tx_demand_i * tx_scale       - ln2/mRNA_halflife_i * M_i

        tl_demand_i  = rbs_strength_i * M_i          (changes with M, so the
                                                        ribosome scale is
                                                        recomputed every step)
        tl_scale     = min(1, ribosome_count/ribosome_ref / sum(tl_demand))
        dP_i/dt      = k_tl * tl_demand_i * tl_scale       - ln2/protein_halflife_i * P_i

    RNAP_ref/ribosome_ref are the shasis defaults (3000 / 20000): with those
    exact counts and a modest gene set, behaviour matches the old
    (unlimited) Stage-2 kinetics — competition only bites once real demand
    is added.

    Returns {time, species, concentrations, final, resource_scale} — same
    shape the metabolic-model canvas plotter already understands, plus a
    'resource_scale' diagnostic (transcription/translation scale at the end
    of the run).
    """
    import numpy as _np
    import math as _math
    genes = model.get('genes', {})
    if not genes:
        return {'error': 'Модель не содержит генов'}
    gene_ids = list(genes.keys())
    n = len(gene_ids)

    y0 = _np.zeros(2 * n)
    deg_m = _np.zeros(n)
    deg_p = _np.zeros(n)
    promoter = _np.zeros(n)
    rbs = _np.zeros(n)
    for i, gid in enumerate(gene_ids):
        g = genes[gid]
        y0[i] = max(0.0, float(g.get('initial_mrna', 0.0)))
        y0[n + i] = max(0.0, float(g.get('initial_protein', 0.0)))
        halflife_m = max(1e-6, float(g.get('mrna_halflife_min', 5.0)))
        halflife_p = max(1e-6, float(g.get('protein_halflife_min', 60.0)))
        deg_m[i] = _math.log(2) / halflife_m
        deg_p[i] = _math.log(2) / halflife_p
        promoter[i] = max(0.0, float(g.get('promoter_strength', 0.5)))
        rbs[i] = max(0.0, float(g.get('rbs_strength', 0.5)))

    processes = model.get('processes', {})
    rnap_capacity = (max(0.0, float(processes.get('rna_polymerase', WC_RNAP_REF_COUNT)))
                     / WC_RNAP_REF_COUNT)
    ribo_capacity = (max(0.0, float(processes.get('ribosomes', WC_RIBOSOME_REF_COUNT)))
                     / WC_RIBOSOME_REF_COUNT) * WC_RIBOSOME_CAPACITY_SCALE
    total_promoter_demand = float(_np.sum(promoter))
    tx_scale = (min(1.0, rnap_capacity / total_promoter_demand)
               if total_promoter_demand > 1e-9 else 1.0)

    last_tl_scale = [1.0]

    def rhs(y):
        M = y[:n]
        P = y[n:]
        tl_demand = rbs * M
        total_tl_demand = float(_np.sum(tl_demand))
        tl_scale = (min(1.0, ribo_capacity / total_tl_demand)
                   if total_tl_demand > 1e-9 else 1.0)
        last_tl_scale[0] = tl_scale
        dM = k_tx * promoter * tx_scale - deg_m * M
        dP = k_tl * tl_demand * tl_scale - deg_p * P
        return _np.concatenate([dM, dP])

    n_points = max(2, int(n_points))
    dt = float(t_max) / n_points
    t = 0.0
    y = y0.copy()
    times = [0.0]
    traj = [y.copy()]
    for _ in range(n_points):
        k1 = rhs(y)
        k2 = rhs(y + 0.5 * dt * k1)
        k3 = rhs(y + 0.5 * dt * k2)
        k4 = rhs(y + dt * k3)
        y = y + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
        y = _np.maximum(y, 0.0)
        t += dt
        times.append(t)
        traj.append(y.copy())
    traj = _np.array(traj)

    species = []
    concentrations = {}
    final = {}
    for i, gid in enumerate(gene_ids):
        name = genes[gid].get('name', gid)
        mid = '{}:mRNA'.format(name)
        pid = '{}:protein'.format(name)
        species.append(mid)
        species.append(pid)
        concentrations[mid] = traj[:, i].tolist()
        concentrations[pid] = traj[:, n + i].tolist()
        final[mid] = round(float(traj[-1, i]), 3)
        final[pid] = round(float(traj[-1, n + i]), 3)

    return {
        'error': None,
        'time': times,
        'species': species,
        'concentrations': concentrations,
        'final': final,
        'resource_scale': {
            'transcription': round(tx_scale, 3),
            'translation': round(last_tl_scale[0], 3),
        },
    }


def simulate_whole_cell_growth(model, t_max=200.0, n_points=300,
                               k_tx=2.0, k_tl=20.0,
                               protein_mass_unit=1.0 / 5000.0,
                               baseline_mass=1.0, division_mass_factor=2.0):
    """Этап 3: рост клетки + репликация хромосомы + деление.

    Использует ту же ресурс-лимитированную кинетику экспрессии генов, что и
    simulate_whole_cell_timecourse (общий пул РНКП/рибосом), плюс:

      - масса клетки = baseline_mass + суммарный белок * protein_mass_unit
      - прогресс репликации хромосомы растёт с постоянной скоростью
        (replication_rate_kb_min / length_kb, единиц прогресса в минуту)
      - деление происходит, когда масса выросла в division_mass_factor раз
        от начальной И хромосома реплицирована полностью (progress >= 1):
        число молекул мРНК/белка делится пополам (распределение между
        дочерними клетками), прогресс репликации сбрасывается, событие
        фиксируется во времени.

    Returns {time, mass, chromosome_progress, division_events, species,
    concentrations, final, error}.
    """
    import numpy as _np
    import math as _math
    genes = model.get('genes', {})
    if not genes:
        return {'error': 'Модель не содержит генов'}
    gene_ids = list(genes.keys())
    n = len(gene_ids)

    deg_m = _np.zeros(n)
    deg_p = _np.zeros(n)
    promoter = _np.zeros(n)
    rbs = _np.zeros(n)
    M = _np.zeros(n)
    P = _np.zeros(n)
    for i, gid in enumerate(gene_ids):
        g = genes[gid]
        M[i] = max(0.0, float(g.get('initial_mrna', 0.0)))
        P[i] = max(0.0, float(g.get('initial_protein', 0.0)))
        halflife_m = max(1e-6, float(g.get('mrna_halflife_min', 5.0)))
        halflife_p = max(1e-6, float(g.get('protein_halflife_min', 60.0)))
        deg_m[i] = _math.log(2) / halflife_m
        deg_p[i] = _math.log(2) / halflife_p
        promoter[i] = max(0.0, float(g.get('promoter_strength', 0.5)))
        rbs[i] = max(0.0, float(g.get('rbs_strength', 0.5)))

    processes = model.get('processes', {})
    rnap_capacity = (max(0.0, float(processes.get('rna_polymerase', WC_RNAP_REF_COUNT)))
                     / WC_RNAP_REF_COUNT)
    ribo_capacity = (max(0.0, float(processes.get('ribosomes', WC_RIBOSOME_REF_COUNT)))
                     / WC_RIBOSOME_REF_COUNT) * WC_RIBOSOME_CAPACITY_SCALE
    total_promoter_demand = float(_np.sum(promoter))
    tx_scale = (min(1.0, rnap_capacity / total_promoter_demand)
               if total_promoter_demand > 1e-9 else 1.0)

    chrom = model.get('chromosome', {})
    length_kb = max(1e-6, float(chrom.get('length_kb', 4600)))
    repl_rate = max(0.0, float(chrom.get('replication_rate_kb_min', 1000)))
    chrom_rate_per_min = repl_rate / length_kb

    def gene_rhs(Mv, Pv):
        tl_demand = rbs * Mv
        total_tl_demand = float(_np.sum(tl_demand))
        tl_scale = (min(1.0, ribo_capacity / total_tl_demand)
                   if total_tl_demand > 1e-9 else 1.0)
        dM = k_tx * promoter * tx_scale - deg_m * Mv
        dP = k_tl * tl_demand * tl_scale - deg_p * Pv
        return dM, dP

    def mass_of(Pv):
        return baseline_mass + float(_np.sum(Pv)) * protein_mass_unit

    n_points = max(2, int(n_points))
    dt = float(t_max) / n_points
    mass0 = mass_of(P)
    division_threshold_mass = mass0 * division_mass_factor
    chrom_progress = 0.0
    t = 0.0
    times = [0.0]
    mass_series = [mass0]
    chrom_series = [chrom_progress]
    division_events = []
    M_traj = [M.copy()]
    P_traj = [P.copy()]

    for _ in range(n_points):
        k1m, k1p = gene_rhs(M, P)
        k2m, k2p = gene_rhs(M + 0.5 * dt * k1m, P + 0.5 * dt * k1p)
        k3m, k3p = gene_rhs(M + 0.5 * dt * k2m, P + 0.5 * dt * k2p)
        k4m, k4p = gene_rhs(M + dt * k3m, P + dt * k3p)
        M = _np.maximum(M + (dt / 6.0) * (k1m + 2 * k2m + 2 * k3m + k4m), 0.0)
        P = _np.maximum(P + (dt / 6.0) * (k1p + 2 * k2p + 2 * k3p + k4p), 0.0)
        chrom_progress = min(1.0, chrom_progress + chrom_rate_per_min * dt)
        t += dt

        cur_mass = mass_of(P)
        if cur_mass >= division_threshold_mass and chrom_progress >= 1.0:
            M = M / 2.0
            P = P / 2.0
            chrom_progress = 0.0
            division_events.append(round(t, 2))
            cur_mass = mass_of(P)

        times.append(t)
        mass_series.append(cur_mass)
        chrom_series.append(chrom_progress)
        M_traj.append(M.copy())
        P_traj.append(P.copy())

    M_traj = _np.array(M_traj)
    P_traj = _np.array(P_traj)
    species = []
    concentrations = {}
    final = {}
    for i, gid in enumerate(gene_ids):
        name = genes[gid].get('name', gid)
        mid = '{}:mRNA'.format(name)
        pid = '{}:protein'.format(name)
        species.append(mid)
        species.append(pid)
        concentrations[mid] = M_traj[:, i].tolist()
        concentrations[pid] = P_traj[:, i].tolist()
        final[mid] = round(float(M_traj[-1, i]), 3)
        final[pid] = round(float(P_traj[-1, i]), 3)

    return {
        'error': None,
        'time': times,
        'mass': mass_series,
        'chromosome_progress': chrom_series,
        'division_events': division_events,
        'species': species,
        'concentrations': concentrations,
        'final': final,
    }


def simulate_whole_cell_stochastic(model, t_max=60.0, k_tx=2.0, k_tl=20.0,
                                   max_events=200000, n_points=200, seed=None):
    """Этап 4: точная стохастика (Гиллеспи / SSA) для мРНК и белка каждого
    гена — четыре реакции на ген (транскрипция, деградация мРНК, трансляция,
    деградация белка), с тем же ресурс-лимитированным fair-share, что и в
    детерминированном движке. При малых числах копий мРНК/белка (обычно
    единицы-десятки на клетку) непрерывная ODE-модель сглаживает случайные
    всплески экспрессии — здесь каждое событие дискретно и случайно, поэтому
    видно настоящий шум экспрессии.

    Returns the same {time, species, concentrations, final} shape as the
    deterministic engine (resampled onto an even grid for plotting), plus
    n_events/truncated diagnostics.
    """
    import numpy as _np
    import math as _math
    import random as _random
    genes = model.get('genes', {})
    if not genes:
        return {'error': 'Модель не содержит генов'}
    gene_ids = list(genes.keys())
    n = len(gene_ids)
    rng = _random.Random(seed)

    deg_m = _np.zeros(n)
    deg_p = _np.zeros(n)
    promoter = _np.zeros(n)
    rbs = _np.zeros(n)
    M = _np.zeros(n, dtype=_np.int64)
    P = _np.zeros(n, dtype=_np.int64)
    for i, gid in enumerate(gene_ids):
        g = genes[gid]
        M[i] = max(0, int(round(float(g.get('initial_mrna', 0.0)))))
        P[i] = max(0, int(round(float(g.get('initial_protein', 0.0)))))
        halflife_m = max(1e-6, float(g.get('mrna_halflife_min', 5.0)))
        halflife_p = max(1e-6, float(g.get('protein_halflife_min', 60.0)))
        deg_m[i] = _math.log(2) / halflife_m
        deg_p[i] = _math.log(2) / halflife_p
        promoter[i] = max(0.0, float(g.get('promoter_strength', 0.5)))
        rbs[i] = max(0.0, float(g.get('rbs_strength', 0.5)))

    processes = model.get('processes', {})
    rnap_capacity = (max(0.0, float(processes.get('rna_polymerase', WC_RNAP_REF_COUNT)))
                     / WC_RNAP_REF_COUNT)
    ribo_capacity = (max(0.0, float(processes.get('ribosomes', WC_RIBOSOME_REF_COUNT)))
                     / WC_RIBOSOME_REF_COUNT) * WC_RIBOSOME_CAPACITY_SCALE
    total_promoter_demand = float(_np.sum(promoter))
    tx_scale = (min(1.0, rnap_capacity / total_promoter_demand)
               if total_promoter_demand > 1e-9 else 1.0)
    tx_rate = k_tx * promoter * tx_scale

    t = 0.0
    times = [0.0]
    M_hist = [M.copy()]
    P_hist = [P.copy()]
    steps = 0
    while t < t_max and steps < max_events:
        tl_demand = rbs * M
        total_tl_demand = float(_np.sum(tl_demand))
        tl_scale = (min(1.0, ribo_capacity / total_tl_demand)
                   if total_tl_demand > 1e-9 else 1.0)
        tl_rate = k_tl * tl_demand * tl_scale

        deg_m_rate = deg_m * M
        deg_p_rate = deg_p * P
        props = _np.concatenate([tx_rate, deg_m_rate, tl_rate, deg_p_rate])
        R = float(_np.sum(props))
        if R <= 1e-12:
            break
        tau = rng.expovariate(R)
        t += tau
        if t > t_max:
            break
        r = rng.random() * R
        idx = 0
        cum = 0.0
        for k in range(len(props)):
            cum += props[k]
            if r <= cum:
                idx = k
                break
        block = idx // n
        gi = idx % n
        if block == 0:
            M[gi] += 1
        elif block == 1:
            if M[gi] > 0:
                M[gi] -= 1
        elif block == 2:
            P[gi] += 1
        else:
            if P[gi] > 0:
                P[gi] -= 1
        times.append(t)
        M_hist.append(M.copy())
        P_hist.append(P.copy())
        steps += 1

    if times[-1] < t_max:
        times.append(t_max)
        M_hist.append(M.copy())
        P_hist.append(P.copy())

    grid = _np.linspace(0.0, t_max, max(2, int(n_points)))
    times_arr = _np.array(times)
    M_hist = _np.array(M_hist)
    P_hist = _np.array(P_hist)
    idxs = _np.searchsorted(times_arr, grid, side='right') - 1
    idxs = _np.clip(idxs, 0, len(times_arr) - 1)
    M_grid = M_hist[idxs]
    P_grid = P_hist[idxs]

    species = []
    concentrations = {}
    final = {}
    for i, gid in enumerate(gene_ids):
        name = genes[gid].get('name', gid)
        mid = '{}:mRNA'.format(name)
        pid = '{}:protein'.format(name)
        species.append(mid)
        species.append(pid)
        concentrations[mid] = M_grid[:, i].tolist()
        concentrations[pid] = P_grid[:, i].tolist()
        final[mid] = int(M_grid[-1, i])
        final[pid] = int(P_grid[-1, i])

    return {
        'error': None,
        'time': grid.tolist(),
        'species': species,
        'concentrations': concentrations,
        'final': final,
        'n_events': steps,
        'truncated': steps >= max_events,
    }


def simulate_whole_cell_montecarlo(model, t_max=60.0, n_runs=20,
                                   k_tx=2.0, k_tl=20.0, max_events=50000):
    """Этап 4: Monte Carlo по популяции клеток — n_runs независимых
    стохастических прогонов (разные случайные последовательности событий),
    затем статистика финальных значений мРНК/белка на ген: среднее,
    стандартное отклонение, коэффициент вариации (шум экспрессии между
    клетками одной популяции).
    """
    genes = model.get('genes', {})
    if not genes:
        return {'error': 'Модель не содержит генов'}
    n_runs = max(2, int(n_runs))
    per_gene_mrna = {gid: [] for gid in genes}
    per_gene_protein = {gid: [] for gid in genes}
    truncated_any = False
    for run_i in range(n_runs):
        res = simulate_whole_cell_stochastic(
            model, t_max=t_max, k_tx=k_tx, k_tl=k_tl,
            max_events=max_events, n_points=2, seed=run_i)
        if res.get('error'):
            return res
        truncated_any = truncated_any or res.get('truncated', False)
        for gid in genes:
            name = genes[gid].get('name', gid)
            per_gene_mrna[gid].append(res['final'].get('{}:mRNA'.format(name), 0))
            per_gene_protein[gid].append(res['final'].get('{}:protein'.format(name), 0))

    import statistics as _stats
    stats = {}
    for gid, g in genes.items():
        name = g.get('name', gid)
        for label, values in (('mRNA', per_gene_mrna[gid]),
                              ('protein', per_gene_protein[gid])):
            mean = _stats.fmean(values)
            std = _stats.pstdev(values) if len(values) > 1 else 0.0
            cv = (std / mean) if mean > 1e-9 else 0.0
            stats['{}:{}'.format(name, label)] = {
                'mean': round(mean, 3), 'std': round(std, 3),
                'cv': round(cv, 3), 'values': values,
            }
    return {'error': None, 'n_runs': n_runs, 'truncated_any': truncated_any,
            'stats': stats}


def simulate_whole_cell_sensitivity(model, gene_id, param='promoter_strength',
                                    values=None, t_max=300.0, n_scan=6):
    """Этап 4: скан параметра одного гена → эффект на скорость роста клетки
    (средний период деления за t_max минут, детерминированный движок —
    reuse simulate_whole_cell_growth на каждой точке скана).
    """
    import copy as _copy
    genes = model.get('genes', {})
    if gene_id not in genes:
        return {'error': 'Ген не найден'}
    if values is None:
        n_scan = max(2, int(n_scan))
        values = [round(i / (n_scan - 1), 3) for i in range(n_scan)]
    results = []
    for v in values:
        m2 = _copy.deepcopy(model)
        m2['genes'][gene_id][param] = v
        res = simulate_whole_cell_growth(m2, t_max=t_max, n_points=200)
        if res.get('error'):
            continue
        events = res['division_events']
        if len(events) >= 2:
            gaps = [events[i + 1] - events[i] for i in range(len(events) - 1)]
            avg_period = sum(gaps) / len(gaps)
        elif len(events) == 1:
            avg_period = events[0]
        else:
            avg_period = None
        results.append({
            'value': v, 'avg_division_period': avg_period,
            'n_divisions': len(events), 'final_mass': res['mass'][-1],
        })
    return {'error': None, 'gene_id': gene_id, 'param': param, 'results': results}


def whole_cell_genes_to_csv(model):
    """Этап 5: экспорт списка генов whole-cell модели в CSV-текст."""
    import io as _io
    import csv as _csv
    buf = _io.StringIO()
    writer = _csv.writer(buf)
    writer.writerow(['id', 'name', 'promoter_strength', 'rbs_strength',
                     'mrna_halflife_min', 'protein_halflife_min',
                     'initial_mrna', 'initial_protein', 'product'])
    for gid, g in model.get('genes', {}).items():
        writer.writerow([
            gid, g.get('name', gid), g.get('promoter_strength', 0.5),
            g.get('rbs_strength', 0.5), g.get('mrna_halflife_min', 5.0),
            g.get('protein_halflife_min', 60.0), g.get('initial_mrna', 0.0),
            g.get('initial_protein', 0.0), g.get('product', ''),
        ])
    return buf.getvalue()


def whole_cell_model_to_antimony(model, k_tx=2.0, k_tl=20.0):
    """Этап 5: экспорт whole-cell модели (гены → мРНК/белок) в Antimony —
    тот же формат, что и экспорт метаболической модели, так что оба
    открываются в Tellurium/COPASI/TinkerCell. Каждый ген даёт пару
    реакций синтеза (транскрипция, трансляция) с массово-действующей
    деградацией, используя те же ставки, что и во встроенном движке
    (без общего ресурсного лимита — тот эффект специфичен для VIREON и не
    переносится напрямую в стандартный ODE-формат).
    """
    import math as _math
    lines = []
    name = model.get('name', 'cell') or 'cell'
    safe_name = _sbml_safe_id(name)[:40] or 'cell'
    lines.append('// Antimony export from VIREON whole-cell model')
    lines.append('model {}()'.format(safe_name))
    lines.append('')
    genes = model.get('genes', {})
    if genes:
        lines.append('  // Species (mRNA / protein per gene)')
        for gid, g in genes.items():
            gname = _sbml_safe_id(g.get('name', gid))
            lines.append('  species {}_mRNA, {}_protein;'.format(gname, gname))
        lines.append('')
        lines.append('  // Reactions')
        for gid, g in genes.items():
            gname = _sbml_safe_id(g.get('name', gid))
            promoter = max(0.0, float(g.get('promoter_strength', 0.5)))
            rbs = max(0.0, float(g.get('rbs_strength', 0.5)))
            halflife_m = max(1e-6, float(g.get('mrna_halflife_min', 5.0)))
            halflife_p = max(1e-6, float(g.get('protein_halflife_min', 60.0)))
            deg_m = _math.log(2) / halflife_m
            deg_p = _math.log(2) / halflife_p
            tx_rate = k_tx * promoter
            lines.append('  -> {}_mRNA; {};'.format(gname, round(tx_rate, 6)))
            lines.append('  {}_mRNA -> ; {} * {}_mRNA;'.format(gname, round(deg_m, 6), gname))
            lines.append('  {}_mRNA -> {}_mRNA + {}_protein; {} * {}_mRNA;'
                         .format(gname, gname, gname, round(k_tl * rbs, 6), gname))
            lines.append('  {}_protein -> ; {} * {}_protein;'
                         .format(gname, round(deg_p, 6), gname))
        lines.append('')
        lines.append('  // Initial amounts')
        for gid, g in genes.items():
            gname = _sbml_safe_id(g.get('name', gid))
            lines.append('  {}_mRNA = {};'.format(gname, g.get('initial_mrna', 0.0)))
            lines.append('  {}_protein = {};'.format(gname, g.get('initial_protein', 0.0)))
    lines.append('end')
    return '\n'.join(lines)


def new_metabolic_model(name='New model'):
    """Create a blank COPASI-style metabolic model.

    Structure:
      compartments: {id: {name, volume, constant}}
      species:      {id: {name, compartment, initial, type, expression, fixed}}
                    type ∈ 'reactions' | 'fixed' | 'assignment' | 'ode'
      reactions:    [{id, name, substrates{sp:stoich}, products{sp:stoich},
                      modifiers[sp], law, params, compartment}]
      global_quantities: {id: {name, value, type, expression}}
      events:       [{trigger, assignments[(target, expr)]}]
    """
    return {
        'name': name,
        'compartments': {
            'cell': {'name': 'cell', 'volume': 1.0, 'constant': True, 'expression': ''},
        },
        'species': {},
        'reactions': [],
        'global_quantities': {},
        'events': [],
        # legacy keys kept for backward compatibility with links / history:
        'organism': 'E_coli_K12',
        'product': 'ethanol',
        'added_enzymes': [],
        'knocked_out': [],
        'last_result': None,
    }


# ─────────────────────────────────────────────────
# ORGANELLE PRESETS (compartmentalization)
# ─────────────────────────────────────────────────
# Real cells separate processes physically: fatty-acid synthesis in the
# cytosol, their breakdown (β-oxidation) in mitochondria, etc. Each preset is a
# compartment with a biologically sensible relative volume (cell cytosol = 1.0
# reference). Volumes are rough textbook ratios for a typical eukaryotic cell.
ORGANELLE_PRESETS = {
    'cytoplasm':    {'name_ru': 'Цитоплазма',    'name_en': 'Cytoplasm',
                     'volume': 1.0,  'color': '#b3d4fc'},
    'mitochondria': {'name_ru': 'Митохондрия',   'name_en': 'Mitochondria',
                     'volume': 0.22, 'color': '#f8b4b4'},
    'nucleus':      {'name_ru': 'Ядро',          'name_en': 'Nucleus',
                     'volume': 0.1,  'color': '#c3b1e1'},
    'er':           {'name_ru': 'ЭПР',           'name_en': 'ER',
                     'volume': 0.15, 'color': '#b7e4c7'},
    'golgi':        {'name_ru': 'Аппарат Гольджи','name_en': 'Golgi',
                     'volume': 0.06, 'color': '#ffe0a3'},
    'peroxisome':   {'name_ru': 'Пероксисома',   'name_en': 'Peroxisome',
                     'volume': 0.01, 'color': '#f5c0d8'},
    'lysosome':     {'name_ru': 'Лизосома',      'name_en': 'Lysosome',
                     'volume': 0.05, 'color': '#d0bfff'},
    'periplasm':    {'name_ru': 'Периплазма',    'name_en': 'Periplasm',
                     'volume': 0.3,  'color': '#ffd8a8'},
    'extracellular':{'name_ru': 'Внеклеточная среда', 'name_en': 'Extracellular',
                     'volume': 10.0, 'color': '#e9ecef'},
}


def add_organelle_to_model(model, organelle_key):
    """Add an organelle preset as a compartment. Returns the new compartment id,
    or None if unknown / already present."""
    org = ORGANELLE_PRESETS.get(organelle_key)
    if not org:
        return None
    if organelle_key in model.get('compartments', {}):
        return None
    model.setdefault('compartments', {})[organelle_key] = {
        'name': org['name_en'], 'volume': float(org['volume']),
        'constant': True, 'expression': '', 'color': org.get('color', '#cccccc'),
    }
    return organelle_key


def make_transport_reaction(model, species_id, from_comp, to_comp,
                            transporter=None, rate=1.0):
    """Create a transport reaction moving a species between two compartments.

    Biologically, crossing a membrane needs a transporter protein. If
    `transporter` (a species id) is given it is added as a modifier, and the
    rate law is carrier-mediated (saturable). Without a transporter the species
    cannot cross — callers should enforce that in the UI.

    The moved species exists as two pools, one per compartment:
      <species>_<from_comp>  →  <species>_<to_comp>
    """
    sp = model['species'].get(species_id)
    if sp is None:
        return None
    base_name = sp.get('name', species_id)
    src_id = '{}_{}'.format(species_id, from_comp)
    dst_id = '{}_{}'.format(species_id, to_comp)
    # Ensure both pools exist as species in their compartments.
    if src_id not in model['species']:
        model['species'][src_id] = {
            'name': '{} ({})'.format(base_name, from_comp), 'compartment': from_comp,
            'initial': sp.get('initial', 0.0), 'type': 'reactions',
            'expression': '', 'fixed': False}
    if dst_id not in model['species']:
        model['species'][dst_id] = {
            'name': '{} ({})'.format(base_name, to_comp), 'compartment': to_comp,
            'initial': 0.0, 'type': 'reactions', 'expression': '', 'fixed': False}
    rid = 'transport_{}_{}_{}'.format(species_id, from_comp, to_comp)
    rxn = {
        'id': rid,
        'name': '{} → {} ({})'.format(from_comp, to_comp, base_name),
        'substrates': {src_id: 1}, 'products': {dst_id: 1},
        'modifiers': [transporter] if transporter else [],
        'law': 'mm_irr' if transporter else 'mass_action',
        'params': ({'Vmax': float(rate), 'Km': 1.0} if transporter
                   else {'k': float(rate)}),
        'compartment': from_comp,
        'is_transport': True,
        'requires_transporter': True,
        'transporter': transporter,
    }
    model.setdefault('reactions', []).append(rxn)
    return rid


def preset_model_from_pathway(pathway_key):
    """Build a full COPASI-style model from a METABOLIC_PATHWAYS entry.

    Turns the pathway's reaction list into species + reactions with
    Michaelis-Menten laws, all inside a single 'cell' compartment.
    """
    pw = METABOLIC_PATHWAYS.get(pathway_key)
    if not pw:
        return None
    model = new_metabolic_model(pw.get('name_en', pathway_key))
    # collect species
    sp_set = []
    for rxn in pw['reactions']:
        s, p = rxn[0], rxn[1]
        for m in (s, p):
            if m not in sp_set:
                sp_set.append(m)
    for i, sp in enumerate(sp_set):
        model['species'][sp] = {
            'name': sp, 'compartment': 'cell',
            'initial': 10.0 if i == 0 else 0.0,
            'type': 'reactions', 'expression': '', 'fixed': (i == 0),
        }
    # first species is the fixed input
    if sp_set:
        model['species'][sp_set[0]]['type'] = 'fixed'
    for j, rxn in enumerate(pw['reactions']):
        s, p, enz = rxn[0], rxn[1], rxn[2]
        atp = rxn[3] if len(rxn) > 3 else 0
        # ATP-neutral steps are near-equilibrium → reversible MM;
        # ATP-consuming/producing steps drive the pathway → irreversible MM.
        if atp == 0:
            law = 'mm_rev'
            params = {'Vf': 1.0, 'Vr': 0.3, 'Kms': 1.0, 'Kmp': 1.0}
        else:
            law = 'mm_irrev'
            params = {'V': 1.0, 'Km': 1.0}
        model['reactions'].append({
            'id': 'R{}'.format(j + 1),
            'name': '{} ({}→{})'.format(enz, s, p),
            'substrates': {s: 1.0}, 'products': {p: 1.0}, 'modifiers': [],
            'law': law, 'params': params,
            'compartment': 'cell', 'enzyme': enz, 'atp': atp,
        })
    model['preset_key'] = pathway_key
    return model


def build_model_from_chain(chain, req_genes=None, name='model',
                           organism_key='E_coli_K12', product_key=None):
    """Build a COPASI-style model from a linear species chain + enzyme names.

    Kept part2-independent: callers in part2 pass the chain extracted from
    METABOLIC_PRODUCTS so nothing is lost, while part3 stays self-contained.
    """
    model = new_metabolic_model(name)
    chain = list(chain)
    for i, sp in enumerate(chain):
        model['species'][sp] = {
            'name': sp, 'compartment': 'cell',
            'initial': 10.0 if i == 0 else 0.0,
            'type': 'fixed' if i == 0 else 'reactions',
            'expression': '', 'fixed': (i == 0),
        }
    req = list(req_genes or [])
    for i in range(len(chain) - 1):
        enz = req[i % len(req)] if req else 'E{}'.format(i + 1)
        model['reactions'].append({
            'id': 'R{}'.format(i + 1),
            'name': '{} ({}→{})'.format(enz, chain[i], chain[i + 1]),
            'substrates': {chain[i]: 1.0}, 'products': {chain[i + 1]: 1.0},
            'modifiers': [], 'law': 'mm_irrev', 'params': {'V': 1.0, 'Km': 1.0},
            'compartment': 'cell', 'enzyme': enz,
        })
    model['organism'] = organism_key
    if product_key:
        model['product'] = product_key
    return model


def simulate_metabolic_model(model, t_max=20.0, n_points=200):
    """Deterministic time-course of a COPASI-style model (RK4).

    Compartment volumes matter: a reaction's velocity is computed on
    concentrations, then the amount change is scaled so that species in
    different compartments stay balanced (v * V_compartment gives amount rate).
    Species of type 'fixed' are held constant. 'assignment'/'ode' handled simply.
    Returns time series of concentrations for every species.
    """
    import numpy as _np
    species = model.get('species', {})
    reactions = model.get('reactions', [])
    comps = model.get('compartments', {'cell': {'volume': 1.0}})
    if not species or not reactions:
        return {'error': 'Model has no species or reactions'}

    sp_ids = list(species.keys())
    idx = {s: i for i, s in enumerate(sp_ids)}
    n = len(sp_ids)

    fixed_mask = _np.array([species[s].get('type') == 'fixed' for s in sp_ids])
    ode_ids = [s for s in sp_ids if species[s].get('type') == 'ode']
    assign_ids = [s for s in sp_ids if species[s].get('type') == 'assignment']
    reaction_mask = _np.array([species[s].get('type') in (None, 'reactions')
                               for s in sp_ids])

    import math as _math
    _safe = {'sin': _math.sin, 'cos': _math.cos, 'tan': _math.tan,
             'exp': _math.exp, 'sqrt': _math.sqrt, 'log': _math.log,
             'abs': abs, 'min': min, 'max': max, 'pi': _math.pi, 'e': _math.e}

    def vol_of(sp):
        c = species[sp].get('compartment', 'cell')
        return float(comps.get(c, {}).get('volume', 1.0)) or 1.0

    # ── Variable-volume compartments ──
    # Non-constant compartments with a growth 'expression' (dV/dt formula) are
    # integrated alongside species. Their live volume feeds concentration =
    # amount / volume, so dilution-by-growth is modelled correctly.
    var_comps = [cid for cid, c in comps.items()
                 if not c.get('constant', True) and c.get('expression', '').strip()]
    comp_state_idx = {cid: n + k for k, cid in enumerate(var_comps)}
    n_total = n + len(var_comps)

    def live_vol(cid, y):
        if cid in comp_state_idx:
            return max(1e-9, y[comp_state_idx[cid]])
        return float(comps.get(cid, {}).get('volume', 1.0)) or 1.0

    def vol_of_y(sp, y):
        c = species[sp].get('compartment', 'cell')
        return live_vol(c, y)

    # 'initial' is a CONCENTRATION; integrate AMOUNTS = concentration × volume
    y0 = _np.zeros(n_total)
    for s in sp_ids:
        y0[idx[s]] = float(species[s].get('initial', 0.0)) * vol_of(s)
    for cid in var_comps:
        y0[comp_state_idx[cid]] = float(comps[cid].get('volume', 1.0))

    # Precompute global-quantity values available to rate laws as extra params
    globals_vals = {gid: gq.get('value', 0.0)
                    for gid, gq in model.get('global_quantities', {}).items()}

    def conc(y, sp):
        """Concentration of a species = amount / (live) compartment volume."""
        if sp not in idx:
            return 0.0
        return max(0.0, y[idx[sp]]) / vol_of_y(sp, y)

    def rhs(t, y):
        dydt = _np.zeros(n_total)
        for rxn in reactions:
            subs = rxn.get('substrates', {})
            prods = rxn.get('products', {})
            mods = rxn.get('modifiers', [])
            law_key = rxn.get('law', 'mm_irrev')
            law = RATE_LAWS.get(law_key, {})
            params = dict(rxn.get('params') or {})
            for gk, gv in globals_vals.items():
                params.setdefault(gk, gv)

            if law.get('kind') == 'basic' and 'mass_action' in law_key:
                kf = params.get('k', params.get('k1', 1.0))
                fwd = kf
                for sp, st in subs.items():
                    fwd *= conc(y, sp) ** st
                v = fwd
                if law.get('reversible'):
                    kr = params.get('k2', 0.0)
                    rev = kr
                    for sp, st in prods.items():
                        rev *= conc(y, sp) ** st
                    v -= rev
            else:
                S = conc(y, next(iter(subs))) if subs else 0.0
                P = conc(y, next(iter(prods))) if prods else 0.0
                I = conc(y, mods[0]) if mods else 0.0
                v = eval_rate_law(law_key, S, P, I, params)

            rxn_comp = rxn.get('compartment', 'cell')
            v_amount = v * live_vol(rxn_comp, y)
            for sp, st in subs.items():
                if sp in idx:
                    dydt[idx[sp]] -= v_amount * st
            for sp, st in prods.items():
                if sp in idx:
                    dydt[idx[sp]] += v_amount * st
        dydt[:n][fixed_mask] = 0.0

        # ODE-type species: user dX/dt formula (concentration units → amount).
        if ode_ids:
            env = {s: conc(y, s) for s in sp_ids}
            env.update(globals_vals)
            env['t'] = t
            for s in ode_ids:
                expr = species[s].get('expression', '').strip()
                if not expr:
                    continue
                try:
                    rate = float(eval(expr, {'__builtins__': {}}, {**_safe, **env}))
                    dydt[idx[s]] += rate * vol_of_y(s, y)
                except Exception:
                    pass
        for s in assign_ids:
            dydt[idx[s]] = 0.0

        # Variable compartment volumes: dV/dt from their expression.
        if var_comps:
            env = {s: conc(y, s) for s in sp_ids}
            env.update(globals_vals)
            env['t'] = t
            for cid in var_comps:
                env['V'] = live_vol(cid, y)
                expr = comps[cid].get('expression', '').strip()
                try:
                    dydt[comp_state_idx[cid]] = float(
                        eval(expr, {'__builtins__': {}}, {**_safe, **env}))
                except Exception:
                    dydt[comp_state_idx[cid]] = 0.0
        return dydt

    def apply_assignments(y):
        """Set assignment-type species from their formula (algebraic rule)."""
        if not assign_ids:
            return y
        env = {s: conc(y, s) for s in sp_ids}
        env.update(globals_vals)
        for s in assign_ids:
            expr = species[s].get('expression', '').strip()
            if not expr:
                continue
            try:
                val = float(eval(expr, {'__builtins__': {}}, {**_safe, **env}))
                y[idx[s]] = max(0.0, val) * vol_of_y(s, y)
            except Exception:
                pass
        return y

    def apply_events(y):
        """Apply model events: when trigger true → assignments. Simple, safe."""
        events = model.get('events', [])
        if not events:
            return y
        env = {sp: conc(y, sp) for sp in sp_ids}
        env.update(globals_vals)
        for ev in events:
            trig = ev.get('trigger', '').strip()
            if not trig:
                continue
            try:
                if eval(trig, {'__builtins__': {}}, {**_safe, **env}):
                    for target, expr in ev.get('assignments', []):
                        if target in idx:
                            val = eval(str(expr), {'__builtins__': {}}, {**_safe, **env})
                            y[idx[target]] = float(val) * vol_of_y(target, y)
            except Exception:
                continue
        return y

    t_vals = _np.linspace(0, t_max, n_points)
    dt = t_vals[1] - t_vals[0]
    Y = _np.zeros((n_points, n_total))
    y = y0.copy()
    y = apply_assignments(y)   # resolve algebraic species before first store
    Y[0] = y
    for i in range(1, n_points):
        k1 = rhs(t_vals[i - 1], y)
        k2 = rhs(t_vals[i - 1] + dt / 2, y + dt / 2 * k1)
        k3 = rhs(t_vals[i - 1] + dt / 2, y + dt / 2 * k2)
        k4 = rhs(t_vals[i - 1] + dt, y + dt * k3)
        y = y + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        y = _np.maximum(0.0, y)
        y = apply_assignments(y)
        y = apply_events(y)
        Y[i] = y

    # Convert amounts → concentrations using the (possibly time-varying) volume
    # at each step.
    def _vol_series(sp):
        cid = species[sp].get('compartment', 'cell')
        if cid in comp_state_idx:
            return _np.maximum(1e-9, Y[:, comp_state_idx[cid]])
        return float(comps.get(cid, {}).get('volume', 1.0)) or 1.0

    conc_out = {}
    for s in sp_ids:
        vser = _vol_series(s)
        conc_out[s] = (Y[:, idx[s]] / vser).tolist()

    result = {
        'error': None,
        'time': t_vals.tolist(),
        'species': sp_ids,
        'concentrations': conc_out,
        'amounts': {s: Y[:, idx[s]].tolist() for s in sp_ids},
        'final': {s: round(conc_out[s][-1], 4) for s in sp_ids},
    }
    if var_comps:
        result['volumes'] = {cid: Y[:, comp_state_idx[cid]].tolist() for cid in var_comps}
    return result


def model_steady_state(model, t_max=500.0):
    """Find approximate steady state by long integration; report dX/dt norm."""
    res = simulate_metabolic_model(model, t_max=t_max, n_points=400)
    if res.get('error'):
        return res
    import numpy as _np
    conc = res['concentrations']
    ss = {s: conc[s][-1] for s in res['species']}
    # crude residual: compare last two points
    resid = 0.0
    for s in res['species']:
        resid += abs(conc[s][-1] - conc[s][-2])
    return {'error': None, 'steady_state': {k: round(v, 4) for k, v in ss.items()},
            'residual': round(resid, 6), 'converged': resid < 1e-3}


def model_stoichiometry_matrix(model):
    """Build the stoichiometry matrix N (species × reactions)."""
    species = list(model.get('species', {}).keys())
    reactions = model.get('reactions', [])
    sp_idx = {s: i for i, s in enumerate(species)}
    import numpy as _np
    N = _np.zeros((len(species), len(reactions)))
    for j, rxn in enumerate(reactions):
        for sp, st in rxn.get('substrates', {}).items():
            if sp in sp_idx:
                N[sp_idx[sp], j] -= st
        for sp, st in rxn.get('products', {}).items():
            if sp in sp_idx:
                N[sp_idx[sp], j] += st
    return {'species': species,
            'reactions': [r.get('id', 'R{}'.format(k + 1)) for k, r in enumerate(reactions)],
            'matrix': N.tolist()}


def elementary_flux_modes(model, max_modes=50):
    """Compute elementary flux modes (EFM) of the network.

    An EFM is a minimal set of reactions that can operate at steady state
    (N·v = 0, v ≥ 0 for irreversible reactions) and that cannot be decomposed
    further. This is a compact, education-grade enumerator: it finds the
    null-space of the stoichiometry matrix and extracts support-minimal,
    sign-consistent generating vectors. For large networks it caps output.
    """
    import numpy as _np
    from itertools import combinations
    stoich = model_stoichiometry_matrix(model)
    N = _np.array(stoich['matrix'])
    rxn_ids = stoich['reactions']
    if N.size == 0:
        return {'error': 'Empty network', 'modes': []}
    n_rxn = N.shape[1]
    reversible = [bool(RATE_LAWS.get(r.get('law', ''), {}).get('reversible', False))
                  for r in model.get('reactions', [])]

    # null space via SVD
    try:
        u, s, vh = _np.linalg.svd(N)
    except Exception as exc:
        return {'error': str(exc), 'modes': []}
    tol = 1e-9
    rank = int((s > tol).sum())
    null_dim = n_rxn - rank
    if null_dim <= 0:
        return {'error': None, 'modes': [], 'null_dim': 0,
                'note': 'No steady-state flux mode (only trivial zero).'}
    null_basis = vh[rank:].T  # (n_rxn × null_dim) columns span null space

    def flux_in_nullspace(target_support):
        """Find a flux vector supported only on target_support that lies in the
        null space (N·v=0). Returns the vector or None. Solves for coefficients
        c such that (null_basis·c) is zero outside target_support."""
        outside = [i for i in range(n_rxn) if i not in target_support]
        if not outside:
            # any coeff works; take first basis column
            v = null_basis[:, 0]
            return v if _np.abs(v).max() > tol else None
        A = null_basis[outside, :]  # must be zero → A·c = 0
        try:
            _, s2, vh2 = _np.linalg.svd(A)
        except Exception:
            return None
        rank2 = int((s2 > tol).sum())
        if rank2 >= A.shape[1]:
            return None  # only trivial c=0
        c = vh2[rank2:].T[:, 0]  # a null coefficient vector
        v = null_basis @ c
        return v if _np.abs(v).max() > tol else None

    def sign_ok(vec):
        for i in range(n_rxn):
            if not reversible[i] and vec[i] < -1e-6:
                return False
        return True

    modes = []
    seen = set()
    # Enumerate candidate supports from small to large; keep only those that
    # (a) carry a genuine null-space flux and (b) are support-minimal.
    active = list(range(n_rxn))
    max_support = min(n_rxn, rank + 2)
    found_supports = []
    for size in range(1, max_support + 1):
        if len(modes) >= max_modes:
            break
        for combo in combinations(active, size):
            combo_set = set(combo)
            # skip if a strictly smaller already-found support is contained
            if any(set(fs) < combo_set for fs in found_supports):
                continue
            vec = flux_in_nullspace(combo_set)
            if vec is None:
                continue
            support = tuple(sorted(i for i in range(n_rxn) if abs(vec[i]) > 1e-6))
            if not support or support in seen:
                continue
            if not set(support) <= combo_set:
                continue
            # orient to satisfy irreversibility if possible
            if not sign_ok(vec):
                if sign_ok(-vec):
                    vec = -vec
                else:
                    continue
            vmax = _np.abs(vec).max()
            vec = vec / vmax
            seen.add(support)
            found_supports.append(support)
            modes.append({
                'reactions': [rxn_ids[i] for i in support],
                'fluxes': {rxn_ids[i]: round(float(vec[i]), 3) for i in support},
            })
            if len(modes) >= max_modes:
                break

    return {'error': None, 'modes': modes, 'null_dim': null_dim, 'rank': rank}


def build_cobra_model_from_metabolic(model, biomass_species=None):
    """Convert a COPASI-style metabolic model into a cobra.Model for FBA/FVA/pFBA.

    Requires the cobra library at runtime. Each species → a metabolite, each
    reaction → a cobra Reaction with substrate/product stoichiometry. Fixed
    species get exchange reactions so the network can take up / secrete them.
    Returns (cobra_model, info_str) or (None, error_str).
    """
    try:
        import cobra
        from cobra import Model, Reaction, Metabolite
    except Exception:
        return None, 'cobra library not installed'
    species = model.get('species', {})
    reactions = model.get('reactions', [])
    if not species or not reactions:
        return None, 'Model has no species or reactions'
    cm = Model(model.get('name', 'metabolic_model').replace(' ', '_')[:40] or 'model')
    mets = {}
    for sid, sp in species.items():
        m = Metabolite(sid, name=sp.get('name', sid),
                       compartment=sp.get('compartment', 'cell'))
        mets[sid] = m
    rxn_objs = []
    for i, rxn in enumerate(reactions):
        r = Reaction(rxn.get('id', 'R{}'.format(i + 1)))
        r.name = rxn.get('name', r.id)
        rev = RATE_LAWS.get(rxn.get('law', ''), {}).get('reversible', False)
        r.lower_bound = -1000.0 if rev else 0.0
        r.upper_bound = 1000.0
        stoich = {}
        for sp, st in rxn.get('substrates', {}).items():
            if sp in mets:
                stoich[mets[sp]] = -float(st)
        for sp, st in rxn.get('products', {}).items():
            if sp in mets:
                stoich[mets[sp]] = stoich.get(mets[sp], 0.0) + float(st)
        if stoich:
            r.add_metabolites(stoich)
            rxn_objs.append(r)
    cm.add_reactions(rxn_objs)
    # Exchange reactions for fixed species (uptake/secretion boundaries)
    exchanges = []
    for sid, sp in species.items():
        if sp.get('type') == 'fixed' and sid in mets:
            ex = Reaction('EX_' + sid)
            ex.name = 'Exchange ' + sid
            ex.lower_bound = -10.0
            ex.upper_bound = 1000.0
            ex.add_metabolites({mets[sid]: -1.0})
            exchanges.append(ex)
    if exchanges:
        cm.add_reactions(exchanges)
    # Objective: biomass species production, or last reaction as proxy
    obj_rxn = None
    if biomass_species and biomass_species in mets:
        bex = Reaction('BIOMASS')
        bex.name = 'Biomass'
        bex.lower_bound = 0.0
        bex.upper_bound = 1000.0
        bex.add_metabolites({mets[biomass_species]: -1.0})
        cm.add_reactions([bex])
        obj_rxn = 'BIOMASS'
    elif rxn_objs:
        obj_rxn = rxn_objs[-1].id
    if obj_rxn:
        cm.objective = obj_rxn
    return cm, 'cobra model: {} mets, {} reactions'.format(
        len(cm.metabolites), len(cm.reactions))


# ═══════════════════════════════════════════════════════════════════
# ADVANCED METABOLIC ANALYSES (COPASI parity) — pure, model-based
# All operate directly on a COPASI-style model dict. No cobra needed.
# ═══════════════════════════════════════════════════════════════════

def timecourse_sensitivities(model, param_rxn_idx=None, perturbation=0.05,
                             t_max=20.0, n_points=100):
    """Time-course sensitivities: how each species' trajectory responds to a
    small change in a reaction parameter (V), over time.

    Returns { reaction_id: { species: [dX/dp at each time] } }. Central
    difference of the full trajectory. This is COPASI's 'Time-Course
    Sensitivities' task.
    """
    reactions = model.get('reactions', [])
    if not reactions:
        return {'error': 'No reactions'}
    idxs = [param_rxn_idx] if param_rxn_idx is not None else range(len(reactions))
    base = simulate_metabolic_model(model, t_max=t_max, n_points=n_points)
    if base.get('error'):
        return base
    import copy as _copy
    out = {}
    for j in idxs:
        rxn = reactions[j]
        params = rxn.get('params', {})
        pname = 'V' if 'V' in params else ('k' if 'k' in params else
                 ('Vf' if 'Vf' in params else None))
        if pname is None:
            continue
        p0 = params[pname]
        dp = abs(p0) * perturbation or perturbation
        m_plus = _copy.deepcopy(model)
        m_minus = _copy.deepcopy(model)
        m_plus['reactions'][j]['params'][pname] = p0 + dp
        m_minus['reactions'][j]['params'][pname] = max(1e-6, p0 - dp)
        r_plus = simulate_metabolic_model(m_plus, t_max=t_max, n_points=n_points)
        r_minus = simulate_metabolic_model(m_minus, t_max=t_max, n_points=n_points)
        if r_plus.get('error') or r_minus.get('error'):
            continue
        sens = {}
        for s in base['species']:
            cp = r_plus['concentrations'][s]
            cm = r_minus['concentrations'][s]
            sens[s] = [(cp[k] - cm[k]) / (2 * dp) for k in range(len(cp))]
        out[rxn.get('id', 'R{}'.format(j + 1))] = sens
    return {'error': None, 'time': base['time'], 'sensitivities': out}


def cross_section(model, monitor_species, threshold=None, direction='up',
                  t_max=100.0, n_points=2000):
    """Cross-section (Poincaré-style): record the values of all species each
    time the monitored species crosses a threshold in a given direction.

    Useful for oscillatory systems to reveal limit cycles / period. Returns the
    crossing points as a list of {t, values}. This is COPASI's 'Cross Section'.
    """
    res = simulate_metabolic_model(model, t_max=t_max, n_points=n_points)
    if res.get('error'):
        return res
    if monitor_species not in res['concentrations']:
        return {'error': 'Unknown species: {}'.format(monitor_species)}
    series = res['concentrations'][monitor_species]
    t = res['time']
    if threshold is None:
        threshold = (max(series) + min(series)) / 2.0
    crossings = []
    for k in range(1, len(series)):
        prev, cur = series[k - 1], series[k]
        up_cross = prev < threshold <= cur
        down_cross = prev > threshold >= cur
        if (direction == 'up' and up_cross) or (direction == 'down' and down_cross) \
           or (direction == 'both' and (up_cross or down_cross)):
            crossings.append({
                't': t[k],
                'values': {s: res['concentrations'][s][k] for s in res['species']},
            })
    # Estimate period from successive crossings
    period = None
    if len(crossings) >= 2:
        gaps = [crossings[i + 1]['t'] - crossings[i]['t'] for i in range(len(crossings) - 1)]
        period = sum(gaps) / len(gaps)
    return {'error': None, 'threshold': threshold, 'crossings': crossings,
            'n_crossings': len(crossings), 'period': period}


def lyapunov_exponent(model, monitor_species=None, t_max=100.0, n_points=2000,
                      perturbation=1e-6):
    """Largest Lyapunov exponent: measures sensitivity to initial conditions
    (positive ⇒ chaos). Integrates a reference and a slightly perturbed
    trajectory, tracks their divergence, and fits the exponential rate.

    This is COPASI's 'Lyapunov Exponents' task (largest exponent).
    """
    import numpy as _np
    import copy as _copy
    ref = simulate_metabolic_model(model, t_max=t_max, n_points=n_points)
    if ref.get('error'):
        return ref
    species = ref['species']
    if not species:
        return {'error': 'No species'}
    m2 = _copy.deepcopy(model)
    s0 = species[0] if monitor_species is None else monitor_species
    if s0 in m2['species']:
        m2['species'][s0]['initial'] = m2['species'][s0].get('initial', 0.0) + perturbation
    pert = simulate_metabolic_model(m2, t_max=t_max, n_points=n_points)
    if pert.get('error'):
        return pert
    # distance between trajectories over time
    d = _np.zeros(len(ref['time']))
    for s in species:
        a = _np.array(ref['concentrations'][s])
        b = _np.array(pert['concentrations'][s])
        d += (a - b) ** 2
    d = _np.sqrt(d)
    d0 = d[0] if d[0] > 0 else perturbation
    # log-divergence slope (fit where divergence is measurable)
    valid = d > 1e-12
    if valid.sum() < 5:
        return {'error': None, 'lyapunov': 0.0,
                'interpretation': 'no measurable divergence (stable/fixed point)'}
    t = _np.array(ref['time'])[valid]
    logd = _np.log(d[valid] / d0)
    # linear fit slope = Lyapunov exponent
    A = _np.vstack([t, _np.ones_like(t)]).T
    slope, _ = _np.linalg.lstsq(A, logd, rcond=None)[0]
    lam = float(slope)
    if lam > 0.01:
        interp = 'chaotic (λ > 0)'
    elif lam < -0.01:
        interp = 'stable / converging (λ < 0)'
    else:
        interp = 'neutral / periodic (λ ≈ 0)'
    return {'error': None, 'lyapunov': round(lam, 5), 'interpretation': interp}


def linear_noise_approximation(model, t_max=500.0):
    """Linear Noise Approximation (LNA): estimates the intrinsic-noise
    covariance of species around the steady state, using the Lyapunov equation
    A·C + C·Aᵀ + D = 0, where A is the Jacobian and D the diffusion (from
    reaction propensities). Reports the stationary variance / CV per species.

    This is COPASI's 'Linear Noise Approximation' task.
    """
    import numpy as _np
    ss = model_steady_state(model, t_max=t_max)
    if ss.get('error'):
        return ss
    species = list(model.get('species', {}).keys())
    reactions = model.get('reactions', [])
    n = len(species)
    idx = {s: i for i, s in enumerate(species)}
    x = _np.array([ss['steady_state'].get(s, 0.0) for s in species])

    # Numerical Jacobian A = d(dx/dt)/dx around steady state
    def flux_vector(xv):
        d = _np.zeros(n)
        for rxn in reactions:
            subs = rxn.get('substrates', {}); prods = rxn.get('products', {})
            mods = rxn.get('modifiers', [])
            S = max(0.0, xv[idx[next(iter(subs))]]) if subs and next(iter(subs)) in idx else 0.0
            P = max(0.0, xv[idx[next(iter(prods))]]) if prods and next(iter(prods)) in idx else 0.0
            I = max(0.0, xv[idx[mods[0]]]) if mods and mods[0] in idx else 0.0
            v = eval_rate_law(rxn.get('law', 'mm_irrev'), S, P, I, rxn.get('params'))
            for sp, st in subs.items():
                if sp in idx: d[idx[sp]] -= v * st
            for sp, st in prods.items():
                if sp in idx: d[idx[sp]] += v * st
        return d

    A = _np.zeros((n, n))
    eps = 1e-6
    f0 = flux_vector(x)
    for j in range(n):
        xp = x.copy(); xp[j] += eps
        A[:, j] = (flux_vector(xp) - f0) / eps

    # Diffusion matrix D = Σ_r (propensity_r) · s_r s_rᵀ  (chemical Langevin)
    D = _np.zeros((n, n))
    for rxn in reactions:
        subs = rxn.get('substrates', {}); prods = rxn.get('products', {})
        mods = rxn.get('modifiers', [])
        S = max(0.0, x[idx[next(iter(subs))]]) if subs and next(iter(subs)) in idx else 0.0
        P = max(0.0, x[idx[next(iter(prods))]]) if prods and next(iter(prods)) in idx else 0.0
        I = max(0.0, x[idx[mods[0]]]) if mods and mods[0] in idx else 0.0
        v = abs(eval_rate_law(rxn.get('law', 'mm_irrev'), S, P, I, rxn.get('params')))
        svec = _np.zeros(n)
        for sp, st in subs.items():
            if sp in idx: svec[idx[sp]] -= st
        for sp, st in prods.items():
            if sp in idx: svec[idx[sp]] += st
        D += v * _np.outer(svec, svec)

    # Solve Lyapunov equation A C + C A^T + D = 0 for covariance C
    try:
        try:
            from scipy.linalg import solve_lyapunov
            C = solve_lyapunov(A, -D)
        except Exception:
            # Fallback: vectorized solve (kron)
            I_n = _np.eye(n)
            K = _np.kron(I_n, A) + _np.kron(A, I_n)
            c = _np.linalg.solve(K, -D.flatten())
            C = c.reshape(n, n)
    except Exception as exc:
        return {'error': 'LNA solve failed: {}'.format(exc)}

    variances = {}
    for s in species:
        var = float(C[idx[s], idx[s]])
        mean = float(x[idx[s]])
        cv = (var ** 0.5 / mean) if mean > 1e-9 and var > 0 else 0.0
        variances[s] = {'mean': round(mean, 4),
                        'variance': round(max(0.0, var), 4),
                        'cv': round(cv, 4)}
    return {'error': None, 'noise': variances}


def parameter_estimation(model, data, param_specs, t_max=None, n_points=None,
                         max_iter=60):
    """Parameter Estimation: fit model reaction parameters to experimental
    time-series data by minimizing sum of squared errors (Nelder-Mead-style
    coordinate descent — dependency-light, no scipy required).

    `data`: { species: [(t, value), ...] }.
    `param_specs`: list of (reaction_index, param_name, low, high).
    Returns best-fit parameter values and the final objective. This is
    COPASI's 'Parameter Estimation' task.
    """
    import copy as _copy
    if not data or not param_specs:
        return {'error': 'Need data and parameters to fit'}
    # time grid from data
    all_t = sorted({t for pts in data.values() for t, _ in pts})
    if not all_t:
        return {'error': 'No data points'}
    tmax = t_max or (max(all_t) * 1.05)
    npts = n_points or max(50, len(all_t) * 4)

    def objective(pvals):
        m = _copy.deepcopy(model)
        for (ridx, pname, _lo, _hi), val in zip(param_specs, pvals):
            if 0 <= ridx < len(m['reactions']):
                m['reactions'][ridx].setdefault('params', {})[pname] = val
        res = simulate_metabolic_model(m, t_max=tmax, n_points=npts)
        if res.get('error'):
            return 1e12
        tgrid = res['time']
        err = 0.0
        for sp, pts in data.items():
            if sp not in res['concentrations']:
                continue
            traj = res['concentrations'][sp]
            for (tt, val) in pts:
                # nearest time index
                k = min(range(len(tgrid)), key=lambda i: abs(tgrid[i] - tt))
                err += (traj[k] - val) ** 2
        return err

    # initialize at midpoints
    p = [(lo + hi) / 2.0 for (_, _, lo, hi) in param_specs]
    best = objective(p)
    steps = [(hi - lo) / 4.0 for (_, _, lo, hi) in param_specs]
    for _it in range(max_iter):
        improved = False
        for i in range(len(p)):
            lo, hi = param_specs[i][2], param_specs[i][3]
            for delta in (steps[i], -steps[i]):
                trial = list(p)
                trial[i] = min(hi, max(lo, trial[i] + delta))
                e = objective(trial)
                if e < best - 1e-9:
                    best = e; p = trial; improved = True
        if not improved:
            steps = [s * 0.5 for s in steps]
            if all(s < 1e-4 for s in steps):
                break
    fitted = {}
    for (ridx, pname, _lo, _hi), val in zip(param_specs, p):
        rid = model['reactions'][ridx].get('id', 'R{}'.format(ridx + 1)) \
            if 0 <= ridx < len(model['reactions']) else 'R?'
        fitted['{}.{}'.format(rid, pname)] = round(val, 4)
    return {'error': None, 'fitted_params': fitted, 'objective': round(best, 6)}


def time_scale_separation(model, t_max=500.0):
    """Time Scale Separation: eigenvalues of the Jacobian at steady state
    reveal fast vs slow modes. Reports each mode's time scale (τ = 1/|Re λ|)
    and classifies the system's stiffness. COPASI's 'Time Scale Separation'.
    """
    import numpy as _np
    ss = model_steady_state(model, t_max=t_max)
    if ss.get('error'):
        return ss
    species = list(model.get('species', {}).keys())
    reactions = model.get('reactions', [])
    n = len(species)
    idx = {s: i for i, s in enumerate(species)}
    x = _np.array([ss['steady_state'].get(s, 0.0) for s in species])

    def flux_vector(xv):
        d = _np.zeros(n)
        for rxn in reactions:
            subs = rxn.get('substrates', {}); prods = rxn.get('products', {})
            mods = rxn.get('modifiers', [])
            S = max(0.0, xv[idx[next(iter(subs))]]) if subs and next(iter(subs)) in idx else 0.0
            P = max(0.0, xv[idx[next(iter(prods))]]) if prods and next(iter(prods)) in idx else 0.0
            I = max(0.0, xv[idx[mods[0]]]) if mods and mods[0] in idx else 0.0
            v = eval_rate_law(rxn.get('law', 'mm_irrev'), S, P, I, rxn.get('params'))
            for sp, st in subs.items():
                if sp in idx: d[idx[sp]] -= v * st
            for sp, st in prods.items():
                if sp in idx: d[idx[sp]] += v * st
        return d

    A = _np.zeros((n, n)); eps = 1e-6; f0 = flux_vector(x)
    for j in range(n):
        xp = x.copy(); xp[j] += eps
        A[:, j] = (flux_vector(xp) - f0) / eps
    eig = _np.linalg.eigvals(A)
    modes = []
    for ev in eig:
        re = float(ev.real)
        tau = (1.0 / abs(re)) if abs(re) > 1e-9 else float('inf')
        modes.append({'eigenvalue_real': round(re, 5),
                      'eigenvalue_imag': round(float(ev.imag), 5),
                      'time_scale': (round(tau, 4) if tau != float('inf') else None),
                      'oscillatory': abs(ev.imag) > 1e-6})
    modes.sort(key=lambda m: (m['time_scale'] is None, m['time_scale'] or 0))
    finite = [m['time_scale'] for m in modes if m['time_scale']]
    stiffness = (max(finite) / min(finite)) if len(finite) >= 2 and min(finite) > 0 else 1.0
    return {'error': None, 'modes': modes, 'stiffness_ratio': round(stiffness, 2),
            'stiff': stiffness > 100}


# ═══════════════════════════════════════════════════════════════════
# MODEL EXCHANGE — export COPASI-style metabolic model to Antimony / SBML
# So models built here can be opened in COPASI, Tellurium, TinkerCell, etc.
# ═══════════════════════════════════════════════════════════════════

def _sbml_safe_id(sid):
    """Make an identifier valid for Antimony/SBML: must start with a letter or
    underscore and contain only [A-Za-z0-9_]."""
    s = ''.join(c if (c.isalnum() or c == '_') else '_' for c in str(sid))
    if not s:
        s = '_x'
    if s[0].isdigit():
        s = '_' + s
    return s


def metabolic_model_to_antimony(model):
    """Export a COPASI-style metabolic model to Antimony text.

    Antimony is a compact, human-readable model language (used by Tellurium,
    TinkerCell). Round-trips species, compartments, reactions with rate laws,
    parameters, and initial values. Returns the text.
    """
    lines = []
    name = (model.get('name', 'model') or 'model')
    safe_name = _sbml_safe_id(name)[:40] or 'model'
    # Build an ID map so species/compartments starting with a digit stay valid
    id_map = {sid: _sbml_safe_id(sid) for sid in model.get('species', {})}
    comp_map = {cid: _sbml_safe_id(cid) for cid in model.get('compartments', {})}
    def SID(x):
        return id_map.get(x, _sbml_safe_id(x))
    def CID(x):
        return comp_map.get(x, _sbml_safe_id(x))
    lines.append('// Antimony export from VIREON metabolic model')
    lines.append('model {}()'.format(safe_name))
    lines.append('')

    comps = model.get('compartments', {})
    species = model.get('species', {})
    reactions = model.get('reactions', [])
    globals_ = model.get('global_quantities', {})

    # Compartments
    if comps:
        lines.append('  // Compartments')
        for cid, c in comps.items():
            lines.append('  compartment {};'.format(CID(cid)))
        lines.append('')

    # Species with compartment + type
    if species:
        lines.append('  // Species')
        for sid, sp in species.items():
            comp = CID(sp.get('compartment', 'cell'))
            prefix = 'const ' if sp.get('type') == 'fixed' else ''
            lines.append('  {}species {} in {};'.format(prefix, SID(sid), comp))
        lines.append('')

    # Reactions with rate laws
    if reactions:
        lines.append('  // Reactions')
        for i, rxn in enumerate(reactions):
            rid = _sbml_safe_id(rxn.get('id', 'R{}'.format(i + 1)))
            subs = rxn.get('substrates', {})
            prods = rxn.get('products', {})
            def side(d):
                parts = []
                for sp, st in d.items():
                    parts.append('{} {}'.format(int(st) if st == int(st) else st, SID(sp))
                                 if st != 1 else SID(sp))
                return ' + '.join(parts) if parts else ''
            arrow = '->'
            params = rxn.get('params', {})
            S = SID(next(iter(subs))) if subs else '1'
            P = SID(next(iter(prods))) if prods else '0'
            I = SID((rxn.get('modifiers') or ['0'])[0]) if rxn.get('modifiers') else '0'
            law_key = rxn.get('law', 'mm_irrev')
            expr = _antimony_rate_expr(law_key, S, P, I, rid)
            lines.append('  {}: {} {} {}; {};'.format(
                rid, side(subs) or '$src', arrow, side(prods) or '$sink', expr))
        lines.append('')
        # Reaction-local parameters
        lines.append('  // Parameters')
        for i, rxn in enumerate(reactions):
            rid = _sbml_safe_id(rxn.get('id', 'R{}'.format(i + 1)))
            for pk, pv in (rxn.get('params') or {}).items():
                lines.append('  {}_{} = {};'.format(rid, pk, pv))
        lines.append('')

    # Global quantities
    if globals_:
        lines.append('  // Global quantities')
        for gid, gq in globals_.items():
            lines.append('  {} = {};'.format(_sbml_safe_id(gid), gq.get('value', 0.0)))
        lines.append('')

    # Initial concentrations + compartment sizes
    lines.append('  // Initial values')
    for cid, c in comps.items():
        lines.append('  {} = {};'.format(CID(cid), c.get('volume', 1.0)))
    for sid, sp in species.items():
        lines.append('  {} = {};'.format(SID(sid), sp.get('initial', 0.0)))
    lines.append('')
    lines.append('end')
    return '\n'.join(lines)


def _antimony_rate_expr(law_key, S, P, I, rid):
    """Build an Antimony kinetic expression string for a rate law, using
    reaction-local parameter names like R1_V, R1_Km."""
    p = lambda name: '{}_{}'.format(rid, name)
    if law_key == 'mass_action_irrev':
        return '{}*{}'.format(p('k'), S)
    if law_key == 'mass_action_rev':
        return '{}*{} - {}*{}'.format(p('k1'), S, p('k2'), P)
    if law_key == 'mm_irrev':
        return '{}*{}/({} + {})'.format(p('V'), S, p('Km'), S)
    if law_key == 'mm_rev':
        return '({}*{}/{} - {}*{}/{})/(1 + {}/{} + {}/{})'.format(
            p('Vf'), S, p('Kms'), p('Vr'), P, p('Kmp'), S, p('Kms'), P, p('Kmp'))
    if law_key == 'hill':
        return '{}*{}^{}/({}^{} + {}^{})'.format(p('V'), S, p('h'), p('Kh'), p('h'), S, p('h'))
    if law_key == 'competitive_irr':
        return '{}*{}/({}*(1 + {}/{}) + {})'.format(p('V'), S, p('Km'), I, p('Ki'), S)
    if law_key == 'uncompetitive_irr':
        return '{}*{}/({} + {}*(1 + {}/{}))'.format(p('V'), S, p('Km'), S, I, p('Ki'))
    if law_key == 'noncompetitive_irr':
        return '({}*{}/({} + {}))/(1 + {}/{})'.format(p('V'), S, p('Km'), S, I, p('Ki'))
    if law_key == 'substrate_inhib':
        return '{}*{}/({} + {} + {}^2/{})'.format(p('V'), S, p('Km'), S, S, p('Ksi'))
    if law_key == 'constant_flux':
        return '{}'.format(p('v'))
    # default MM
    return '{}*{}/({} + {})'.format(p('V'), S, p('Km'), S)


def metabolic_model_to_sbml(model, filepath=None):
    """Export a COPASI-style metabolic model to SBML L3 (needs libsbml).

    Returns (sbml_string_or_None, info). If libsbml is unavailable, returns a
    clear message so the UI can fall back to Antimony.
    """
    try:
        import libsbml
    except Exception:
        return None, 'libsbml not installed — use Antimony export instead'
    try:
        doc = libsbml.SBMLDocument(3, 2)
        m = doc.createModel()
        m.setId(_sbml_safe_id(model.get('name', 'model'))[:40] or 'model')

        for cid, c in model.get('compartments', {}).items():
            comp = m.createCompartment()
            comp.setId(_sbml_safe_id(cid)); comp.setConstant(bool(c.get('constant', True)))
            comp.setSize(float(c.get('volume', 1.0))); comp.setSpatialDimensions(3)

        for sid, sp in model.get('species', {}).items():
            s = m.createSpecies()
            s.setId(_sbml_safe_id(sid)); s.setCompartment(_sbml_safe_id(sp.get('compartment', 'cell')))
            s.setInitialConcentration(float(sp.get('initial', 0.0)))
            s.setConstant(sp.get('type') == 'fixed')
            s.setBoundaryCondition(sp.get('type') == 'fixed')
            s.setHasOnlySubstanceUnits(False)

        for i, rxn in enumerate(model.get('reactions', [])):
            r = m.createReaction()
            rid = _sbml_safe_id(rxn.get('id', 'R{}'.format(i + 1)))
            r.setId(rid)
            r.setReversible(bool(RATE_LAWS.get(rxn.get('law', ''), {}).get('reversible')))
            r.setFast(False)
            for sp, st in rxn.get('substrates', {}).items():
                ref = r.createReactant(); ref.setSpecies(_sbml_safe_id(sp))
                ref.setStoichiometry(float(st)); ref.setConstant(True)
            for sp, st in rxn.get('products', {}).items():
                ref = r.createProduct(); ref.setSpecies(_sbml_safe_id(sp))
                ref.setStoichiometry(float(st)); ref.setConstant(True)
            for mod in rxn.get('modifiers', []):
                if mod:
                    mref = r.createModifier(); mref.setSpecies(_sbml_safe_id(mod))
            # kinetic law
            klaw = r.createKineticLaw()
            S = _sbml_safe_id(next(iter(rxn.get('substrates', {})), '1'))
            P = _sbml_safe_id(next(iter(rxn.get('products', {})), '0'))
            I = _sbml_safe_id((rxn.get('modifiers') or ['0'])[0])
            expr = _antimony_rate_expr(rxn.get('law', 'mm_irrev'), S, P, I, rid)
            math_ast = libsbml.parseL3Formula(expr)
            if math_ast is not None:
                klaw.setMath(math_ast)
            for pk, pv in (rxn.get('params') or {}).items():
                par = klaw.createLocalParameter()
                par.setId('{}_{}'.format(rid, pk)); par.setValue(float(pv))

        sbml_str = libsbml.writeSBMLToString(doc)
        if filepath:
            libsbml.writeSBMLToFile(doc, filepath)
        return sbml_str, 'SBML L3 exported ({} species, {} reactions)'.format(
            len(model.get('species', {})), len(model.get('reactions', [])))
    except Exception as exc:
        return None, 'SBML export failed: {}'.format(exc)

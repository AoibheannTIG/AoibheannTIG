#!/usr/bin/env python3
"""
r1cs_optimize.py

Pipeline:
  Step 1 — eliminate internal (non I/O) variables, iteratively
      Greedy safe linear elimination (mod p), then optional Groebner elimination over GF(p)
      to shrink to the smallest set of higher-order constraints in I/O (often one polynomial).

  Step 2 — algebraic factoring to simplify structure
      Factor each remaining constraint.

  Step 3 — reduce all constraints to the form a*b + c = 0
      Compile each polynomial into R1CS with only one product per constraint by introducing
      tmp variables with constraints of the form: (linA)*(linB) +linC (-tmp) = 0.

Outputs:
  A compact R1CS JSON with the same basic schema (prime, counts, constraints),
  where every constraint is a*b + c = 0 (a,b,c are linear).

Usage:
  python r1cs_optimize.py input_circuit.json -o final_r1cs.json
"""

from __future__ import annotations
import argparse
import json
from dataclasses import dataclass
from typing import Dict, List, Tuple, Set, Iterable, Optional

from sympy import (
    symbols, Symbol, Integer, Poly, Add, Mul, Pow, expand, factor, factor_list, groebner, factor_terms, collect, gcd
)
from sympy.polys.domains import GF


# -----------------------------
# Utilities (mod p arithmetic)
# -----------------------------

def modp(x: int, p: int) -> int:
    return int(x) % p


def addm(a: int, b: int, p: int) -> int:
    return (a + b) % p


def subm(a: int, b: int, p: int) -> int:
    return (a - b) % p


def mulm(a: int, b: int, p: int) -> int:
    return (a * b) % p


def invm(a: int, p: int) -> int:
    # Python's pow handles modular inverse when exponent = -1 (Python 3.8+)
    return pow(int(a) % p, -1, p)


# -----------------------------
# Circuit loading / symbol map
# -----------------------------

@dataclass
class CircuitMeta:
    prime: int
    n_outputs: int
    n_pub: int
    n_prv: int
    n_vars: int  # total variables including 1, outs, ins, internals


def load_circuit(path: str) -> Tuple[CircuitMeta, List[Tuple[Dict[int, int], Dict[int, int], Dict[int, int]]]]:
    with open(path, "r") as f:
        c = json.load(f)

    meta = CircuitMeta(
        prime=int(c["prime"]),
        n_outputs=int(c["nOutputs"]),
        n_pub=int(c["nPubInputs"]),
        n_prv=int(c["nPrvInputs"]),
        n_vars=int(c["nVars"])
    )

    constraints = []
    for A, B, C in c["constraints"]:
        # keys are var indices as strings in some dumps; normalize to int
        constraints.append((
            {int(k): int(v) for k, v in A.items()},
            {int(k): int(v) for k, v in B.items()},
            {int(k): int(v) for k, v in C.items()},
        ))
    return meta, constraints


def build_symbols(meta: CircuitMeta) -> Tuple[List[Symbol], List[Symbol], List[Symbol], List[Symbol]]:
    """
    Construct witness symbol order (by index):
      0: constant "1" (we represent by Integer(1) in expressions)
      1..n_outputs: out0..out{n_outputs-1}
      next n_inputs: in0..in{n_inputs-1}
      remaining:     var0.. (internal variables)
    """
    one = Integer(1)  # not a Symbol; used directly in expressions
    outs = [symbols(f"out{i}") for i in range(meta.n_outputs)]
    n_inputs = meta.n_pub + meta.n_prv
    ins = [symbols(f"in{i}") for i in range(n_inputs)]
    n_internals = meta.n_vars - (1 + meta.n_outputs + n_inputs)
    internals = [symbols(f"var{i}") for i in range(n_internals)]
    # Return as: [1], outs, ins, internals
    return [one], outs, ins, internals


def linear_form_from_dict(d: Dict[int, int], idx2sym: List, p: int):
    """
    Build a SymPy linear expression from a sparse dict {index: coeff}.
    Coefficients are reduced modulo p; index 0 (the constant slot) is added directly.
    """
    expr = Integer(0)
    n = len(idx2sym)
    for raw_idx, raw_coeff in d.items():
        # validate and normalize
        if not isinstance(raw_idx, int):
            raw_idx = int(raw_idx)
        if raw_idx < 0 or raw_idx >= n:
            raise ValueError(f"Variable index {raw_idx} out of range (0..{n-1}).")
        c = modp(int(raw_coeff), p)
        if c == 0:
            continue
        if raw_idx == 0:
            # constant slot: add directly, avoid big-int * 1 path
            expr += Integer(c)
        else:
            sym = idx2sym[raw_idx]
            expr += Integer(c) * sym
    return expr


def constraints_to_polys(raw_constraints, meta: CircuitMeta, one_out_in_var_lists):
    one, outs, ins, internals = one_out_in_var_lists
    # index to symbol, preserving the same layout
    idx2sym = [one] + outs + ins + internals
    polys = []
    p = meta.prime
    for A, B, C in raw_constraints:
        a = linear_form_from_dict(A, idx2sym, p)
        b = linear_form_from_dict(B, idx2sym, p)
        c = linear_form_from_dict(C, idx2sym, p)
        polys.append(expand(a * b - c))
    return polys, idx2sym


# -----------------------------
# Step 1: Eliminate internals
# -----------------------------

def is_pure_linear(expr) -> bool:
    """Return True if expr is affine-linear in all symbols (sum of numbers and 1-degree terms)."""
    e = expand(expr)
    if e.is_Number:
        return True
    if isinstance(e, Add):
        terms = e.args
    else:
        terms = (e,)

    for t in terms:
        if t.is_Number:
            continue
        if isinstance(t, Mul):
            # Check total symbolic degree of product term
            deg = 0
            for f in t.args:
                if f.is_Number:
                    continue
                if isinstance(f, Pow):
                    if not f.exp.is_Integer or int(f.exp) != 1:
                        return False
                    deg += 1
                elif f.is_Symbol:
                    deg += 1
                else:
                    return False
            if deg > 1:
                return False
        elif isinstance(t, Pow):
            if not t.exp.is_Integer or int(t.exp) != 1:
                return False
        elif t.is_Symbol:
            pass
        else:
            return False
    return True


def try_eliminate_one_linear(poly_list: List, internals: List[Symbol], p: int) -> Tuple[bool, List]:
    """
    Try to find a linear definition for an internal var: a*var + rest == 0 with 'a' constant in GF(p).
    Substitute and drop that equation.
    """
    for i, expr in enumerate(poly_list):
        e = expand(expr)
        if not is_pure_linear(e):
            continue
        # Try each internal variable
        for v in internals:
            if v not in e.free_symbols:
                continue
            coeff_v = expand(e).coeff(v)
            # Require constant coefficient (no symbols in coeff)
            if coeff_v.free_symbols:
                continue
            a = int(coeff_v)
            if a % p == 0:
                continue
            # rest = e - a*v
            rest = expand(e - coeff_v * v)
            # Solve v = - a^{-1} * rest  (mod p)
            inv = invm(a, p)
            # Build substitution expression modulo p
            # Note: reduce numeric coefficients in 'rest' modulo p
            rest_mod = _mod_reduce_expr(rest, p)
            rhs = expand((-inv) * rest_mod)
            # Substitute into all other polys
            new_list = []
            for j, ex in enumerate(poly_list):
                if j == i:
                    continue  # drop the defining equation
                new_list.append(expand(ex.subs({v: rhs})))
            return True, new_list
    return False, poly_list


def _mod_reduce_expr(expr, p: int):
    """
    Reduce the integer coefficients of an expression modulo p.
    Works structurally on Add/Mul/Pow/Symbol/Number.
    """
    e = expand(expr)
    if e.is_Number:
        return Integer(modp(int(e), p))
    if isinstance(e, Add):
        return Add(*[_mod_reduce_expr(a, p) for a in e.args])
    if isinstance(e, Mul):
        return Mul(*[_mod_reduce_expr(a, p) for a in e.args])
    if isinstance(e, Pow):
        base = _mod_reduce_expr(e.base, p)
        exp = e.exp
        return base ** exp
    # Symbol or other atom
    return e


def eliminate_internals_iterative(polys: List, internals: List[Symbol], p: int) -> List:
    changed = True
    current = polys[:]
    iteration = 0
    while changed:
        iteration += 1
        print(f"    Iteration {iteration}: trying to eliminate internal variables...")
        changed, current = try_eliminate_one_linear(current, internals, p)
        if changed:
            print(f"    -> Successfully eliminated a variable, {len(current)} equations remaining")
        else:
            print(f"    -> No more variables can be eliminated")
    return current


def groebner_eliminate(polys: List, internals: List[Symbol], ios: List[Symbol], p: int,
                       enable: bool = True) -> List:
    """
    Use Groebner basis over GF(p) with lex order to eliminate internals.
    Order: internals first, then IOs -> G ∩ k[IOs] generates the elimination ideal.
    Returns polynomials that do not contain internals (if any found), else original list.
    """
    if not enable or not polys:
        return polys
    gens = internals + ios
    try:
        print(f"    Computing Groebner basis with {len(polys)} polynomials...")
        G = groebner(polys, *gens, order='lex', domain=GF(p))
        print(f"    Groebner basis has {len(G.polys)} polynomials")
        elim = []
        internal_set = set(internals)
        for g in G.polys:
            ge = g.as_expr()
            if ge.free_symbols.isdisjoint(internal_set):
                elim.append(expand(ge))
        if elim:
            print(f"    Found {len(elim)} polynomials without internal variables")
            # Keep a minimal non-trivial subset (heuristic: shortest string repr first)
            elim.sort(key=lambda e: (len(str(e)), len(e.free_symbols)))
            return elim
        print(f"    No polynomials found without internal variables")
        return polys
    except Exception as e:
        print(f"    Groebner basis computation failed: {e}")
        return polys


# -----------------------------
# Step 2: Factoring (GF(p))
# -----------------------------

def selective_gcd_factor(expr, outputs: List[Symbol], p: int):
    """
    Factor out a common term only from the subset of terms that do NOT involve outputs.
    This creates a partially factored expression where output terms remain separate.
    
    Example: -in0**4*in1 - in0**2*in1 + out0  ->  -in0**2*in1*(in0**2 + 1) + out0
    """
    e = expand(expr)
    
    # Handle trivial cases
    if e.is_Number or e.is_Symbol:
        return e
    
    if not isinstance(e, Add):
        # Single term or already factored; try symbolic factoring first
        try:
            symbolic_fact = factor(e)
            # Only use modular factoring if symbolic fails
            if symbolic_fact == e:
                return factor(e, modulus=p)
            return symbolic_fact
        except:
            return e

    outs_set = set(outputs)
    factorable_terms = []
    output_terms = []
    
    # Separate terms into two groups: those with outputs and those without
    for t in e.args:
        if t.is_Number:
            # Keep constants with factorable terms (they contribute to GCD)
            factorable_terms.append(t)
            continue
        if t.free_symbols & outs_set:
            # Any term touching an output stays separate
            output_terms.append(t)
        else:
            # Terms without outputs can be factored together
            factorable_terms.append(t)

    # If no factorable terms, try symbolic factoring first
    if not factorable_terms:
        try:
            return factor(e)
        except:
            return e
    
    # If all terms are factorable (no outputs), try symbolic factoring first
    if not output_terms:
        try:
            factored = factor(e)
            # Apply R1CS-optimized factorization
            return optimize_for_r1cs(factored)
        except:
            return e

    # Factor only the non-output chunk
    chunk = Add(*factorable_terms)
    
    # Try to extract GCD and factor symbolically first
    try:
        # First try factor_terms to pull out monomial GCD
        chunk_fact = factor_terms(chunk)
        # Then try symbolic factoring
        chunk_fact = factor(chunk_fact)
        # Apply R1CS optimizations
        chunk_fact = optimize_for_r1cs(chunk_fact)
    except Exception:
        # If symbolic factoring fails, just expand
        chunk_fact = expand(chunk)

    # Combine factored chunk with output terms
    result = chunk_fact + Add(*output_terms)
    return result

def collect_desired_powers(polys: List) -> Dict[Symbol, Set[int]]:
    """
    Look through all polys and record powers base**exp that appear,
    so lowering can build intermediate powers that will be reused.
    """
    desired: Dict[Symbol, Set[int]] = {}
    def visit(e):
        if isinstance(e, Pow) and e.base.is_Symbol and e.exp.is_Integer and int(e.exp) >= 2:
            b, k = e.base, int(e.exp)
            desired.setdefault(b, set()).add(k)
        if hasattr(e, 'args'):
            for a in e.args:
                visit(a)
    for p in polys:
        visit(p)
    return desired
def optimize_for_r1cs(expr):
    """
    Optimize expression structure to minimize R1CS constraints by finding
    optimal factorizations that reduce the number of multiplications needed.
    
    Strategy: Analyze all monomials and find common factors that can be shared.
    Example: x^5*y^2 + y^2*z^5 + x^3 + z^3 -> y^2*(x^5 + z^5) + x^3 + z^3
    """
    # If the expression is already factored (like (a-b)*(c+d+e)), 
    # we need to check if the inner expressions can be optimized
    if isinstance(expr, Mul):
        optimized_args = []
        changed = False
        for arg in expr.args:
            if isinstance(arg, Add):
                # Try to optimize this additive sub-expression
                optimized_arg = optimize_for_r1cs(arg)
                if optimized_arg != arg:
                    changed = True
                optimized_args.append(optimized_arg)
            else:
                optimized_args.append(arg)
        
        if changed:
            if len(optimized_args) == 1:
                return optimized_args[0]
            else:
                return Mul(*optimized_args)
    
    # Handle additive expressions
    e = expand(expr)
    if not isinstance(e, Add):
        return expr  # Return original if not additive after expansion
    
    # Extract all monomials (terms)
    terms = e.args
    if len(terms) < 2:
        return expr
    
    # Analyze each term to extract factors
    term_factors = []
    for term in terms:
        factors = analyze_monomial_factors(term)
        term_factors.append((term, factors))
    
    # Find the best factorization that minimizes R1CS cost
    best_factorization = find_optimal_factorization(term_factors)
    
    if best_factorization:
        print(f"    R1CS optimization found better factorization: {expr} -> {best_factorization}")
        return best_factorization
    else:
        return expr


def analyze_monomial_factors(term):
    """
    Analyze a monomial term and extract all possible factors.
    Returns a dict of {factor_expr: remaining_expr} pairs.
    """
    factors = {}
    
    # Handle different term types
    if term.is_Symbol:
        factors[term] = Integer(1)
        factors[Integer(1)] = term
    elif term.is_Number:
        factors[Integer(1)] = term
    elif isinstance(term, Pow) and term.base.is_Symbol and term.exp.is_Integer:
        # Handle power terms like in0**2
        base = term.base
        exp = int(term.exp)
        symbol_powers = {base: exp}
        coeff = 1
        
        # Generate all possible factors for this power term
        symbols = [base]
        from itertools import combinations
        
        # For each symbol, we can include 0 to full_power of it in the factor
        for i in range(len(symbols) + 1):
            for symbol_combo in combinations(symbols, i):
                # For each symbol in the combination, try different power levels
                if symbol_combo:
                    # Try extracting different powers of each symbol
                    _generate_factor_combinations(symbol_combo, symbol_powers, coeff, factors)
                else:
                    # Just the coefficient (which is 1 for pure power terms)
                    factors[Integer(1)] = term
    elif isinstance(term, Mul):
        # Analyze multiplication terms
        coeff = 1
        symbol_powers = {}
        
        for arg in term.args:
            if arg.is_Number:
                coeff *= int(arg)
            elif arg.is_Symbol:
                symbol_powers[arg] = symbol_powers.get(arg, 0) + 1
            elif isinstance(arg, Pow) and arg.base.is_Symbol and arg.exp.is_Integer:
                base, exp = arg.base, int(arg.exp)
                symbol_powers[base] = symbol_powers.get(base, 0) + exp
        
        # Generate all possible factors
        symbols = list(symbol_powers.keys())
        
        # Generate all combinations of symbols and their partial powers
        from itertools import combinations
        
        # For each symbol, we can include 0 to full_power of it in the factor
        for i in range(len(symbols) + 1):
            for symbol_combo in combinations(symbols, i):
                # For each symbol in the combination, try different power levels
                if symbol_combo:
                    # Try extracting different powers of each symbol
                    _generate_factor_combinations(symbol_combo, symbol_powers, coeff, factors)
                elif coeff != 1:
                    # Just the coefficient
                    factors[Integer(coeff)] = Mul(*[sym**power for sym, power in symbol_powers.items()])
    
    return factors


from itertools import product

def _generate_factor_combinations(symbol_combo, symbol_powers, coeff, factors):
    """
    Generate factors that include the chosen symbols with ANY exponent from 1..power,
    and put the leftover exponents into the remaining part. Also carry the numeric coeff.
    """
    # Build the exponent ranges for each chosen symbol (1..full power)
    exp_ranges = [range(1, symbol_powers[sym] + 1) for sym in symbol_combo]

    for exps in product(*exp_ranges):
        factor_parts = []
        remaining_parts = []

        # numeric coefficient stays with the factor (so we can consider pulling consts)
        if coeff != 1:
            factor_parts.append(Integer(coeff))

        # For every symbol in the monomial, split between factor and remainder
        chosen = dict(zip(symbol_combo, exps))
        for sym, power in symbol_powers.items():
            if sym in chosen:
                k = chosen[sym]                # exponent in the factor (1..power)
                r = power - k                  # leftover exponent for the remainder
                # add factor piece
                factor_parts.append(sym if k == 1 else sym**k)
                # add remainder piece if any leftover
                if r > 0:
                    remaining_parts.append(sym if r == 1 else sym**r)
            else:
                # symbol not chosen goes entirely to remainder
                remaining_parts.append(sym if power == 1 else sym**power)

        # build expressions
        factor_expr = factor_parts[0] if len(factor_parts) == 1 else Mul(*factor_parts) if factor_parts else Integer(1)
        remaining_expr = remaining_parts[0] if len(remaining_parts) == 1 else Mul(*remaining_parts) if remaining_parts else Integer(1)

        # record
        factors[factor_expr] = remaining_expr

    # Also consider the pure-coefficient factor (pull only the const)
    if coeff != 1:
        remaining_parts = []
        for sym, power in symbol_powers.items():
            remaining_parts.append(sym if power == 1 else sym**power)
        remaining_expr = remaining_parts[0] if len(remaining_parts) == 1 else Mul(*remaining_parts) if remaining_parts else Integer(1)
        factors[Integer(coeff)] = remaining_expr



def find_optimal_factorization(term_factors):
    """
    Find the factorization that minimizes the estimated R1CS constraint count.
    """
    if len(term_factors) < 2:
        return None
    
    # Find common factors across terms
    common_factors = {}
    
    for term, factors in term_factors:
        for factor, remaining in factors.items():
            if factor == Integer(1):
                continue
            if factor not in common_factors:
                common_factors[factor] = []
            common_factors[factor].append((term, remaining))
    
    # Find factors that appear in multiple terms
    shared_factors = {f: terms for f, terms in common_factors.items() if len(terms) > 1}
    
    if not shared_factors:
        return None
    
    # Score each possible factorization
    best_factor = None
    best_score = float('inf')
    best_factorization = None
    
    for factor, term_pairs in shared_factors.items():
        # Estimate R1CS cost for this factorization
        score = estimate_r1cs_cost(factor, term_pairs, term_factors)
        
        if score < best_score:
            best_score = score
            best_factor = factor
            best_factorization = term_pairs
    
    if best_factor is None:
        return None
    
    # Build the factorized expression
    factored_terms = []
    remaining_terms = []
    used_terms = set()
    
    for term, remaining in best_factorization:
        factored_terms.append(remaining)
        used_terms.add(term)
    
    # Add terms that weren't factored
    for term, _ in term_factors:
        if term not in used_terms:
            remaining_terms.append(term)
    
    # Construct the result
    result_parts = []
    
    if factored_terms:
        if len(factored_terms) == 1:
            inner_expr = factored_terms[0]
        else:
            inner_expr = Add(*factored_terms)
        result_parts.append(best_factor * inner_expr)
    
    result_parts.extend(remaining_terms)
    
    if len(result_parts) == 1:
        return result_parts[0]
    else:
        return Add(*result_parts)


def estimate_r1cs_cost(factor, term_pairs, all_term_factors):
    """
    Estimate the R1CS constraint cost for a given factorization.
    Lower is better.
    """
    # Cost of computing the factor itself
    factor_cost = estimate_expr_cost(factor)
    
    # Cost of computing the inner sum
    inner_terms = [remaining for _, remaining in term_pairs]
    inner_cost = len(inner_terms) - 1  # additions are free, but structure matters
    
    # Cost of the final multiplication
    mult_cost = 1
    
    # Penalty for unused terms (they still need to be computed separately)
    used_terms = {term for term, _ in term_pairs}
    unused_terms = [term for term, _ in all_term_factors if term not in used_terms]
    unused_cost = sum(estimate_expr_cost(term) for term in unused_terms)
    
    return factor_cost + inner_cost + mult_cost + unused_cost


def estimate_expr_cost(expr):
    """
    Estimate the R1CS cost of computing an expression.
    """
    if expr.is_Symbol or expr.is_Number:
        return 0
    elif isinstance(expr, Mul):
        # Count the number of multiplications needed
        non_const_factors = [arg for arg in expr.args if not arg.is_Number]
        if len(non_const_factors) <= 1:
            return 0
        return len(non_const_factors) - 1
    elif isinstance(expr, Pow) and expr.exp.is_Integer:
        exp = int(expr.exp)
        if exp <= 1:
            return 0
        return exp - 1  # Approximate cost
    elif isinstance(expr, Add):
        return sum(estimate_expr_cost(arg) for arg in expr.args)
    else:
        return 1  # Unknown expression type


def extract_monomial_gcd(terms: List, p: int):
    """
    Extract the GCD of a list of polynomial terms.
    Returns (gcd_expr, [term/gcd for term in terms])
    """
    if not terms:
        return Integer(1), []
    
    if len(terms) == 1:
        return Integer(1), terms
    
    # Compute GCD of all terms
    try:
        gcd_expr = terms[0]
        for t in terms[1:]:
            gcd_expr = gcd(gcd_expr, t)
            if gcd_expr == 1:
                # No common factor
                return Integer(1), terms
        
        # Factor out the GCD
        reduced_terms = []
        for t in terms:
            quot = expand(t / gcd_expr)
            reduced_terms.append(quot)
        
        return gcd_expr, reduced_terms
    except Exception:
        return Integer(1), terms


def aggressive_partial_factor(expr, outputs: List[Symbol], p: int):
    """
    Aggressively factor by:
    1. Separating terms with/without outputs
    2. Extracting GCD from non-output terms
    3. Further factoring each group
    4. Reassembling the expression
    """
    e = expand(expr)
    
    if not isinstance(e, Add):
        try:
            return factor(e)
        except:
            return e
    
    outs_set = set(outputs)
    non_output_terms = []
    output_terms = []
    
    # Separate terms
    for t in e.args:
        if t.free_symbols & outs_set:
            output_terms.append(t)
        else:
            non_output_terms.append(t)
    
    # If everything has outputs or nothing has outputs, try symbolic factoring
    if not non_output_terms or not output_terms:
        try:
            return factor(e)
        except:
            return e
    
    # Extract GCD from non-output terms
    gcd_expr, reduced_non_output = extract_monomial_gcd(non_output_terms, p)
    
    # Create factored chunk
    if gcd_expr != 1:
        # Factor the sum of reduced terms
        inner_sum = Add(*reduced_non_output)
        try:
            inner_factored = factor(inner_sum)
            factored_chunk = expand(gcd_expr * inner_factored)
        except Exception:
            factored_chunk = expand(gcd_expr * inner_sum)
    else:
        # No GCD found, try factoring the whole chunk
        chunk = Add(*non_output_terms)
        try:
            factored_chunk = factor(chunk)
        except Exception:
            factored_chunk = chunk
    
    # Combine with output terms
    result = factored_chunk + Add(*output_terms)
    return result


def best_collect(expr, vars_pool: Iterable[Symbol], p: int):
    """
    Heuristic: try collecting by each symbol; pick the one that reduces
    the expression size the most (measured by string length and mul count).
    """
    e0 = expand(expr)
    
    def count_mul(e):
        """Count number of Mul nodes in expression tree"""
        if e.is_Atom:
            return 0
        if isinstance(e, Mul):
            return 1 + sum(count_mul(a) for a in e.args)
        if isinstance(e, Add):
            return sum(count_mul(a) for a in e.args)
        if isinstance(e, Pow):
            return count_mul(e.base)
        return 0

    def score_expr(e):
        """Score expression: prefer fewer Mul nodes and shorter string"""
        return (count_mul(e), len(str(e)))

    best_e = e0
    best_score = score_expr(e0)
    
    for v in vars_pool:
        if v not in e0.free_symbols:
            continue
        try:
            # Try collecting by this variable
            cand = collect(expand(e0), v)
            sc = score_expr(cand)
            if sc < best_score:
                best_e, best_score = cand, sc
        except Exception:
            continue
    
    return best_e


def factor_over_gfp(expr, p: int):
    """Symbolic factoring with fallback to GF(p) if needed"""
    try:
        # Try symbolic factoring first
        symbolic_fact = factor(expr)
        if symbolic_fact != expr:
            return symbolic_fact
        # If symbolic factoring didn't help, try over GF(p)
        return factor(expr, modulus=p)
    except Exception:
        return expand(expr)


def recognize_power_difference_patterns(expr, outputs: List[Symbol] = None):
    """
    Recognize patterns like a^n - b^n and suggest optimal R1CS factorizations.
    Returns:
    - For even powers: (a^n - b^n)(a^n + b^n)
    - For powers divisible by 3: (a^N - b^N) * ( a^N*(a^N + b^N) + b^N*b^N )
    - For other odd powers: returns the original expression to prevent suboptimal factoring
    
    The function also returns a boolean indicating whether this is a power difference
    that should skip other factorization attempts.
    """
    if not isinstance(expr, Add):
        return None, False

    output_terms = []
    power_terms = []

    if outputs is None:
        outputs = []
    output_set = set(outputs)

    for term in expr.args:
        if term.is_Symbol and term in output_set:
            output_terms.append(term)
        elif term.free_symbols & output_set:
            output_terms.append(term)
        else:
            power_terms.append(term)

    if len(power_terms) != 2:
        return None, False

    term_info = []
    for term in power_terms:
        if isinstance(term, Pow) and term.base.is_Symbol and term.exp.is_Integer:
            term_info.append((1, term.base, int(term.exp)))
        elif isinstance(term, Mul) and len(term.args) == 2:
            coeff_part = None
            pow_part = None
            for arg in term.args:
                if arg.is_Number:
                    coeff_part = int(arg)
                elif isinstance(arg, Pow) and arg.base.is_Symbol and arg.exp.is_Integer:
                    pow_part = (arg.base, int(arg.exp))
            if coeff_part is not None and pow_part is not None:
                term_info.append((coeff_part, pow_part[0], pow_part[1]))
        elif term.is_Symbol:
            term_info.append((1, term, 1))
        elif (isinstance(term, Mul) and len(term.args) == 2 and
              term.args[0].is_Number and term.args[1].is_Symbol):
            term_info.append((int(term.args[0]), term.args[1], 1))

    if len(term_info) != 2:
        return None, False

    coeff1, var1, exp1 = term_info[0]
    coeff2, var2, exp2 = term_info[1]

    # Normalize so the negative term is second
    if coeff1 < 0:
        coeff1, var1, exp1, coeff2, var2, exp2 = coeff2, var2, exp2, coeff1, var1, exp1

    # a^n - b^n with equal exponents
    if coeff1 == 1 and coeff2 == -1 and exp1 == exp2 and var1 != var2:
        exp = exp1
        a, b = var1, var2

        # Only factor in specific cases that reduce R1CS complexity:
        
        # 1) Even power: a^(2n) - b^(2n) -> (a^n - b^n)(a^n + b^n)
        if exp % 2 == 0:
            half_exp = exp // 2
            a_half = a**half_exp if half_exp > 1 else a
            b_half = b**half_exp if half_exp > 1 else b
            factored_power_part = (a_half - b_half) * (a_half + b_half)
            print(f"        Recognized even power difference: {a}^{exp} - {b}^{exp} "
                  f"-> ({a_half} - {b_half}) * ({a_half} + {b_half})")
        
        # 2) Powers divisible by 3: a^(3N) - b^(3N) -> (a^N - b^N) * ( a^N*(a^N + b^N) + b^N*b^N )
        elif exp % 3 == 0:
            N = exp // 3
            aN = a**N if N > 1 else a
            bN = b**N if N > 1 else b
            # Keep the structured second factor exactly as requested
            second_factor = aN * (aN + bN) + bN * bN
            factored_power_part = (aN - bN) * second_factor
            print(f"        Recognized power divisible by 3: {a}^{exp} - {b}^{exp} "
                  f"-> ({aN} - {bN}) * ({second_factor})")
        
        # 3) For other odd powers (not divisible by 2 or 3), skip factoring
        # as the standard factorization increases R1CS complexity
        else:
            print(f"        Skipping factorization for odd power {exp} not divisible by 3")
            factored_power_part = None

        # For odd powers not divisible by 2 or 3, return original expr and True to skip other factoring
        if exp % 2 != 0 and exp % 3 != 0:
            return expr, True
        
        if factored_power_part is not None:
            result = (Add(*output_terms) - factored_power_part) if output_terms else factored_power_part
            return result, True  # True indicates this is a power difference pattern

    return None, False  # Not a power difference pattern


def systematic_factorization_search(expr, outputs: List[Symbol] = None):
    """
    Systematically search for the best factorization by trying different combinations of terms.
    
    For an expression like in0**2 + in0*in1 + in0*in2 - in0 - in1**2 + in1*in2,
    this will try factoring:
    1. All 6 terms together
    2. All combinations of 5 terms (leaving 1 out)
    3. All combinations of 4 terms (leaving 2 out)
    4. etc.
    
    Returns the factorization with minimum estimated R1CS cost, or None if no improvement found.
    """
    from itertools import combinations
    
    if not isinstance(expr, Add):
        return None
    
    # Separate output terms from other terms
    if outputs is None:
        outputs = []
    output_set = set(outputs)
    
    output_terms = []
    factorable_terms = []
    
    for term in expr.args:
        if term.is_Symbol and term in output_set:
            output_terms.append(term)
        elif term.free_symbols & output_set:
            output_terms.append(term)
        else:
            factorable_terms.append(term)
    
    if len(factorable_terms) < 3:
        return None  # Need at least 3 terms to make factoring worthwhile
    
    best_factorization = None
    best_cost = estimate_expr_cost(expr)
    
    print(f"        Systematic factorization search on {len(factorable_terms)} terms...")
    
    # Try different numbers of terms to factor together, from all terms down to 3 terms
    for num_terms in range(len(factorable_terms), 2, -1):
        print(f"          Trying combinations of {num_terms} terms...")
        
        # Try all combinations of num_terms
        for term_combo in combinations(factorable_terms, num_terms):
            # Create sub-expression from selected terms
            sub_expr = Add(*term_combo)
            
            try:
                # Try to factor this sub-expression
                factored_sub = factor(sub_expr)
                
                # Skip if factoring didn't change anything
                if factored_sub == sub_expr:
                    continue
                
                # Calculate cost of the factored version
                factored_cost = estimate_expr_cost(factored_sub)
                
                # Add cost of remaining unfactored terms
                remaining_terms = [t for t in factorable_terms if t not in term_combo]
                remaining_cost = sum(estimate_expr_cost(t) for t in remaining_terms)
                
                total_cost = factored_cost + remaining_cost
                
                # If this is better, construct the full expression
                if total_cost < best_cost:
                    result_parts = [factored_sub] + remaining_terms + output_terms
                    candidate = Add(*result_parts) if len(result_parts) > 1 else result_parts[0]
                    
                    print(f"            Found better factorization: {sub_expr} -> {factored_sub}")
                    print(f"            Cost improvement: {best_cost} -> {total_cost}")
                    
                    best_factorization = candidate
                    best_cost = total_cost
                    
                    # If we found a really good factorization, we can stop early
                    if total_cost <= 1:
                        break
                        
            except Exception as e:
                # Factoring failed, continue with next combination
                continue
        
        # If we found a good factorization, no need to try smaller combinations
        if best_factorization is not None and best_cost <= 1:
            break
    
    return best_factorization


def factoring_stage(polys: List, outputs: List[Symbol], ios: List[Symbol], p: int) -> List:
    """
    Enhanced Step 2: Apply multiple factoring strategies
    1. Check for special power difference patterns (x^n - y^n)
    2. For non-special cases:
       a. Selective GCD factoring (separate output terms)
       b. Aggressive partial factoring
       c. Collect by I/O variables for better structure
       d. Full factoring over GF(p)
    
    Special handling for power differences:
    - Even powers: factor as (x^n - y^n)(x^n + y^n)
    - Powers divisible by 3: factor as (x^n - y^n)(x^n(x^n + y^n) + y^n*y^n)
    - Other odd powers: keep as x^n - y^n (no factoring)
    """
    print(f"  Applying factoring strategies to {len(polys)} polynomials...")
    result = []
    
    def expr_complexity(expr):
        """Measure expression complexity - prefer fewer terms and lower string length"""
        if isinstance(expr, Add):
            return (len(expr.args), len(str(expr)))
        return (1, len(str(expr)))
    
    for idx, e in enumerate(polys):
        print(f"    Factoring polynomial {idx}: {e}")
        
        # First check if this is a power difference pattern
        power_pattern_result, is_power_diff = recognize_power_difference_patterns(e, outputs)
        if is_power_diff:
            # For power differences, we only use the pattern-specific factorization
            # or keep as is for odd powers not divisible by 2 or 3
            print(f"      Power pattern recognition: {e} -> {power_pattern_result}")
            best_expr = power_pattern_result
        else:
            # For non-power-difference expressions, try all strategies
            best_expr = e
            best_score = expr_complexity(e)
            
            # Strategy 1: Systematic factorization search
            systematic_result = systematic_factorization_search(best_expr, outputs)
            if systematic_result is not None:
                print(f"      Systematic factorization: {best_expr} -> {systematic_result}")
                score_systematic = expr_complexity(systematic_result)
                if score_systematic < best_score:
                    best_expr, best_score = systematic_result, score_systematic
            
            # Strategy 2: Selective GCD factor
            e2 = selective_gcd_factor(best_expr, outputs, p)
            print(f"      After selective GCD: {e2}")
            score2 = expr_complexity(e2)
            if score2 < best_score:
                best_expr, best_score = e2, score2
            
            # Strategy 3: Aggressive partial factoring
            e3 = aggressive_partial_factor(best_expr, outputs, p)
            print(f"      After aggressive factoring: {e3}")
            score3 = expr_complexity(e3)
            if score3 < best_score:
                best_expr, best_score = e3, score3
            
            # Strategy 4: Try collecting by I/O variables
            e4 = best_collect(best_expr, ios, p)
            print(f"      After collection: {e4}")
            score4 = expr_complexity(e4)
            if score4 < best_score:
                best_expr, best_score = e4, score4
            
            # Strategy 5: Final full factoring pass
            e5 = factor_over_gfp(best_expr, p)
            print(f"      After final factoring: {e5}")
            score5 = expr_complexity(e5)
            if score5 < best_score:
                best_expr = e5
        
        print(f"      Final result: {best_expr}")
        result.append(best_expr)
    
    return result


# -----------------------------
# Linear form representation for R1CS
# -----------------------------

@dataclass
class LinearForm:
    # lin(x) = const + sum coeffs[sym]*sym
    coeffs: Dict[Symbol, int]
    const: int

    @staticmethod
    def zero() -> "LinearForm":
        return LinearForm({}, 0)

    @staticmethod
    def one(p: int) -> "LinearForm":
        return LinearForm({}, 1 % p)

    @staticmethod
    def var(sym: Symbol) -> "LinearForm":
        return LinearForm({sym: 1}, 0)

    @staticmethod
    def const_(c: int, p: int) -> "LinearForm":
        return LinearForm({}, c % p)

    def add(self, other: "LinearForm", p: int) -> "LinearForm":
        out = dict(self.coeffs)
        for s, c in other.coeffs.items():
            out[s] = (out.get(s, 0) + c) % p
        return LinearForm(out, (self.const + other.const) % p)

    def scale(self, k: int, p: int) -> "LinearForm":
        k = k % p
        return LinearForm({s: (c * k) % p for s, c in self.coeffs.items()},
                          (self.const * k) % p)


# -----------------------------
# Step 3: Compilation to R1CS a*b + c = 0
# -----------------------------

@dataclass
class R1CSConstraint:
    A: LinearForm
    B: LinearForm
    C: LinearForm  # A*B + C == 0


@dataclass
class R1CSBuilder:
    p: int
    outs: List[Symbol]
    ins: List[Symbol]
    tmp_counter: int = 0
    # Cache for common subexpressions: (expr_key) -> tmp_symbol
    expr_cache: Optional[Dict[str, Symbol]] = None
    # Track desired powers for optimization
    desired_powers: Optional[Dict[Symbol, Set[int]]] = None

    def __post_init__(self):
        if self.expr_cache is None:
            self.expr_cache = {}
    def _cache_power(self, base: Symbol, exp: int, sym: Symbol):
        if exp >= 2:
            self.expr_cache[f"{base}**{exp}"] = sym
    def _pow_binary_on_lf(self, lf_base: LinearForm, k: int, cons_acc: List[R1CSConstraint]):
        result = LinearForm.const_(1, self.p)
        power = lf_base
        exp_val = k
        made_tmp_for_result = None

        while exp_val > 0:
            if exp_val % 2 == 1:
                if result.const == 1 and not result.coeffs:
                    result = power
                else:
                    t, cons = self.emit_mul(result, power)
                    if cons is not None:
                        cons_acc.append(cons)
                    result = LinearForm.var(t)
                    made_tmp_for_result = t
            exp_val //= 2
            if exp_val > 0:
                t, cons = self.emit_mul(power, power)
                if cons is not None:
                    cons_acc.append(cons)
                power = LinearForm.var(t)

        # If base is a single symbol and result is a tmp, cache a**k
        if len(lf_base.coeffs) == 1 and lf_base.const == 0 and made_tmp_for_result is not None:
            base_sym = list(lf_base.coeffs.keys())[0]
            self._cache_power(base_sym, k, made_tmp_for_result)
        return result, cons_acc
    def new_tmp(self) -> Symbol:
        t = symbols(f"tmp{self.tmp_counter}")
        self.tmp_counter += 1
        return t
    
    def check_cache_for_expr(self, expr) -> Optional[Symbol]:
        """
        Check if an expression is already cached and return the cached variable.
        This handles common patterns like var**2, var1*var2, etc.
        """

        if isinstance(expr, Pow) and expr.base.is_Symbol and expr.exp.is_Integer:
            key = f"{expr.base}**{int(expr.exp)}"
            if key in self.expr_cache:
                return self.expr_cache[key]
        elif isinstance(expr, Mul) and len(expr.args) == 2:
            # Check for var1*var2 pattern
            args = expr.args
            if all(arg.is_Symbol for arg in args):
                if args[0] == args[1]:
                    # var**2 pattern
                    cache_key = f"{args[0]}**2"
                else:
                    # var1*var2 pattern
                    vars_sorted = sorted([str(args[0]), str(args[1])])
                    cache_key = f"{vars_sorted[0]}*{vars_sorted[1]}"
                return self.expr_cache.get(cache_key)
        return None

    def emit_mul(self, left: LinearForm, right: LinearForm) -> Tuple[Symbol, R1CSConstraint]:
        """
        Create t = left * right  ->  left*right + (-t) == 0
        Check cache first to avoid recomputing the same expression.
        """
        # Create a key for caching - cache simple multiplications and powers
        cache_key = None
        if (len(left.coeffs) == 1 and left.const == 0 and 
            len(right.coeffs) == 1 and right.const == 0):
            # Both are single variables
            left_var = list(left.coeffs.keys())[0]
            right_var = list(right.coeffs.keys())[0]
            left_coeff = list(left.coeffs.values())[0]
            right_coeff = list(right.coeffs.values())[0]
            if left_coeff == 1 and right_coeff == 1:
                # Simple var * var multiplication - cache it
                if left_var == right_var:
                    # This is a power: var^2
                    cache_key = f"{left_var}**2"
                else:
                    # Different variables
                    vars_sorted = sorted([str(left_var), str(right_var)])
                    cache_key = f"{vars_sorted[0]}*{vars_sorted[1]}"
                
                if cache_key in self.expr_cache:
                    return self.expr_cache[cache_key], None  # Return cached, no new constraint
        
        t = self.new_tmp()
        cons = R1CSConstraint(
            A=left,
            B=right,
            C=LinearForm.var(t).scale(-1, self.p)
        )
        
        # Cache the result if we have a cache key
        if cache_key:
            self.expr_cache[cache_key] = t
            
        return t, cons

    def emit_lin_eq(self, lhs: LinearForm, rhs: LinearForm) -> R1CSConstraint:
        """
        Encode lhs == rhs  ->  (lhs - rhs) * 1 + 0 == 0
        """
        diff = lhs.add(rhs.scale(-1, self.p), self.p)
        return R1CSConstraint(A=diff, B=LinearForm.one(self.p), C=LinearForm.const_(0, self.p))

    def try_as_linear(self, expr) -> Tuple[Optional[LinearForm], List[R1CSConstraint]]:
        """
        Try to convert expression to LinearForm without introducing tmp variables.
        Returns (LinearForm, constraints) if possible, (None, []) if not linear.
        """
        e = expr
        if e.is_Number:
            return LinearForm.const_(int(e), self.p), []
        if e.is_Symbol:
            return LinearForm.var(e), []
        
        if isinstance(e, Add):
            # Sum of linear terms is linear
            acc = LinearForm.zero()
            all_cons: List[R1CSConstraint] = []
            for term in e.args:
                term_lf, term_cons = self.try_as_linear(term)
                if term_lf is None:
                    return None, []  # Not linear
                acc = acc.add(term_lf, self.p)
                all_cons.extend(term_cons)
            return acc, all_cons
            
        if isinstance(e, Mul):
            # Check if it's just a constant times a symbol
            num = 1
            sym = None
            for arg in e.args:
                if arg.is_Number:
                    num = mulm(num, int(arg), self.p)
                elif arg.is_Symbol:
                    if sym is None:
                        sym = arg
                    else:
                        return None, []  # Multiple symbols = not linear
                else:
                    return None, []  # Complex factor = not linear
            
            if sym is not None:
                lf = LinearForm.var(sym).scale(num, self.p)
                return lf, []
            else:
                # Just a constant
                return LinearForm.const_(num, self.p), []
        
        # Not a linear expression
        return None, []

    # ---- expression lowering ----

    def lower_expr(self, expr) -> Tuple[LinearForm, List[R1CSConstraint]]:
        """
        Lower an arbitrary polynomial 'expr' to a LinearForm by inserting tmp vars and
        constraints so that each product is handled by a single a*b + (-tmp) == 0.
        Returns (LinearForm representing expr, list of emitted constraints).
        """
        e = expr  # Don't expand immediately - preserve factored structure
        if e.is_Number:
            return LinearForm.const_(int(e), self.p), []
        if e.is_Symbol:
            return LinearForm.var(e), []
        
        # Check cache first for common expressions
        cached_var = self.check_cache_for_expr(e)
        if cached_var is not None:
            return LinearForm.var(cached_var), []

        if isinstance(e, Add):
            acc = LinearForm.zero()
            all_cons: List[R1CSConstraint] = []
            for term in e.args:
                lf, cons = self.lower_expr(term)
                acc = acc.add(lf, self.p)
                all_cons.extend(cons)
            return acc, all_cons

        if isinstance(e, Mul):
            # Separate numeric coefficient and non-numeric factors
            num = 1
            factors = []
            for arg in e.args:
                if arg.is_Number:
                    num = mulm(num, int(arg), self.p)
                else:
                    factors.append(arg)

            cons_acc: List[R1CSConstraint] = []
            if not factors:
                return LinearForm.const_(num, self.p), []

            # Special optimization: detect patterns like var^n * (linear_expr)
            # where linear_expr is a sum/difference that can be handled directly
            if len(factors) == 2:
                left, right = factors
                left_lf, left_cons = self.try_as_linear(left)
                right_lf, right_cons = self.try_as_linear(right)
                
                # If one factor is linear, we can multiply directly
                if left_lf is not None and right_lf is None:
                    # left is linear, right needs lowering
                    right_lf, right_cons = self.lower_expr(right)
                    cons_acc.extend(right_cons)
                    t, cons = self.emit_mul(left_lf, right_lf)
                    if cons is not None:
                        cons_acc.append(cons)
                    result_lf = LinearForm.var(t)
                elif right_lf is not None and left_lf is None:
                    # right is linear, left needs lowering  
                    left_lf, left_cons = self.lower_expr(left)
                    cons_acc.extend(left_cons)
                    t, cons = self.emit_mul(left_lf, right_lf)
                    if cons is not None:
                        cons_acc.append(cons)
                    result_lf = LinearForm.var(t)
                elif left_lf is not None and right_lf is not None:
                    # Both are linear - direct multiplication
                    t, cons = self.emit_mul(left_lf, right_lf)
                    if cons is not None:
                        cons_acc.append(cons)
                    result_lf = LinearForm.var(t)
                else:
                    # Neither is linear - fall back to sequential processing
                    result_lf = None
                
                if result_lf is not None:
                    # Apply numeric coefficient
                    if num % self.p != 1:
                        if num % self.p == 0:
                            return LinearForm.const_(0, self.p), cons_acc
                        result_lf = result_lf.scale(num, self.p)
                    return result_lf, cons_acc

            # Fallback: sequential multiplication
            lf_left, cons_left = self.lower_expr(factors[0])
            cons_acc.extend(cons_left)

            current_lf = lf_left
            for f in factors[1:]:
                lf_right, cons_right = self.lower_expr(f)
                cons_acc.extend(cons_right)
                t, cons = self.emit_mul(current_lf, lf_right)
                if cons is not None:
                    cons_acc.append(cons)
                current_lf = LinearForm.var(t)

            # Apply numeric coefficient
            if num % self.p != 1:
                if num % self.p == 0:
                    return LinearForm.const_(0, self.p), cons_acc
                current_lf = current_lf.scale(num, self.p)
            
            return current_lf, cons_acc

        if isinstance(e, Pow):
            base, exp = e.base, e.exp
            if not exp.is_Integer or int(exp) < 0:
                raise ValueError("Negative or non-integer power encountered.")
            k = int(exp)
            if k == 0:
                return LinearForm.const_(1, self.p), []
            if k == 1:
                return self.lower_expr(base)

            if not base.is_Symbol:
                # Fallback: lower base, then binary exponentiation on the tmp
                lf_base, cons_b = self.lower_expr(base)
                return self._pow_binary_on_lf(lf_base, k, cons_b)

            # ---- base is a Symbol: use hint-aware chain ----
            cons_acc: List[R1CSConstraint] = []
            # 1) If exactly cached, reuse
            cached = self.expr_cache.get(f"{base}**{k}")
            if cached is not None:
                return LinearForm.var(cached), []

            # 2) If k even and k//2 is desired elsewhere, build that first, then square
            desired = (self.desired_powers or {}).get(base, set())
            if k % 2 == 0 and (k // 2) in desired:
                lf_half, cons_half = self.lower_expr(Pow(base, k // 2))
                cons_acc.extend(cons_half)
                t, cons = self.emit_mul(lf_half, lf_half)
                if cons is not None:
                    cons_acc.append(cons)
                # cache a**k
                if len(lf_half.coeffs) == 1 and lf_half.const == 0:
                    half_sym = list(lf_half.coeffs.keys())[0]
                    # if half_sym is a tmp symbol, great; cache both
                    self._cache_power(base, k // 2, half_sym)
                # resulting tmp symbol of squaring:
                res_sym = list(cons.C.coeffs.keys())[0] if cons is not None else self.expr_cache.get(f"{base}**{k}")
                if res_sym is not None:
                    self._cache_power(base, k, res_sym)
                return LinearForm.var(res_sym), cons_acc

            # 3) Otherwise, standard binary exponentiation
            lf_base = LinearForm.var(base)
            return self._pow_binary_on_lf(lf_base, k, [])

    def lower_poly_eq_zero(self, expr) -> List[R1CSConstraint]:
        """
        Lower a single polynomial equation expr == 0 to R1CS, without emitting a
        trailing linear constraint. We try to absorb the final linear tie-off into
        the last multiplication constraint whenever possible.

        Special cases (kept):
        - out - complex_expr = 0  -> assign result directly to 'out'
        - (a+b)*(c+d) = linear   -> one constraint with linear on C
        """
        output_vars = set(self.outs)

        # Case 1: out - product_expr = 0  -> write directly into out
        if isinstance(expr, Add) and len(expr.args) == 2:
            pos_terms = []
            neg_terms = []
            for term in expr.args:
                if isinstance(term, Mul) and len(term.args) >= 1 and term.args[0] == -1:
                    if len(term.args) == 2:
                        neg_terms.append(term.args[1])
                    else:
                        neg_terms.append(Mul(*term.args[1:]))
                else:
                    pos_terms.append(term)
            if (len(pos_terms) == 1 and pos_terms[0] in output_vars and len(neg_terms) == 1):
                out_var = pos_terms[0]
                complex_expr = neg_terms[0]
                return self.lower_expr_to_output(complex_expr, out_var)

        # Case 2: factored_product = linear  -> put linear on C directly
        if isinstance(expr, Add):
            linear_terms = []
            factored_products = []
            for term in expr.args:
                if isinstance(term, Mul) and len(term.args) >= 2:
                    has_factored_structure = any(isinstance(arg, Add) for arg in term.args if not arg.is_Number)
                    if has_factored_structure:
                        factored_products.append(term)
                    else:
                        linear_terms.append(term)
                else:
                    linear_terms.append(term)
            if len(factored_products) == 1 and len(linear_terms) > 0:
                factored_product = factored_products[0]
                linear_sum = Add(*linear_terms) if len(linear_terms) > 1 else linear_terms[0]
                return self.lower_factored_product_assignment(factored_product, -linear_sum)

        # Fallback: lower arbitrarily, then absorb linear into the last multiply if we can
        lf, cons = self.lower_expr(expr)

        # Try to absorb linear tie-off into the last multiplication constraint.
        # We look for: last constraint is A*B + (-t_last) == 0, and lf = (± t_last) + L.
        if cons:
            last = cons[-1]
            # Detect last as a standard multiplication with C = -t_last
            if (last.C.const % self.p == 0 and
                len(last.C.coeffs) == 1 and
                list(last.C.coeffs.values())[0] % self.p == (self.p - 1)):
                t_last = list(last.C.coeffs.keys())[0]

                # Coefficient of t_last inside lf (mod p)
                coeff_mod = lf.coeffs.get(t_last, 0) % self.p
                # Interpret as signed: +1 or -1
                if coeff_mod == 1 or coeff_mod == (self.p - 1):
                    # Build L = lf - (±1)*t_last  (i.e., remove the t_last term)
                    # Copy lf without t_last
                    new_coeffs = dict(lf.coeffs)
                    new_coeffs.pop(t_last, None)
                    L = LinearForm(new_coeffs, lf.const)

                    # If lf = t_last + L, then A*B + L == 0
                    # If lf = -t_last + L, then A*B - L == 0  --> C = (-L)
                    if coeff_mod == 1:
                        last.C = L
                    else:
                        last.C = L.scale(-1, self.p)

                    # We successfully absorbed the linear tie-off; return cons (no extra linear eq).
                    return cons

        # Could not absorb; as a last resort, we *still* avoid a separate linear constraint by
        # multiplying by 1:  (lf) * 1 + 0 == 0  (keeps the "A*B + C = 0" shape with B linear)
        # If you want to forbid *any* purely-linear final constraint, comment the next two lines
        # and instead raise or try further restructuring.
        # return [self.emit_lin_eq(lf, LinearForm.const_(0, self.p))]

        # Strict: never emit a linear tail. Try to anchor lf against an existing multiply.
        # If no multiply exists, force one by multiplying with a constant-1 tmp:
        # Create t = 1*1 (no-op) then set A=lf, B=1, C=0 anchored on that constraint.
        # Simpler: inject a dummy multiplication by 1 via new tmp gate:
        one_lf = LinearForm.one(self.p)
        t, mul_cons = self.emit_mul(one_lf, one_lf)  # t = 1*1
        # Now rewrite that gate to encode lf*1 + 0 == 0
        mul_cons.A = lf
        mul_cons.B = one_lf
        mul_cons.C = LinearForm.const_(0, self.p)
        return cons + [mul_cons]
    
    def lower_expr_to_output(self, expr, output_var: Symbol) -> List[R1CSConstraint]:
        """
        Lower expression but assign the final result directly to output_var instead of a tmp.
        """
        if isinstance(expr, Mul) and len(expr.args) >= 2:
            # For multiplication, process all but the last factor normally,
            # then use output_var for the final result
            num = 1
            factors = []
            for arg in expr.args:
                if arg.is_Number:
                    num = mulm(num, int(arg), self.p)
                else:
                    factors.append(arg)
            
            if len(factors) >= 2:
                cons_acc: List[R1CSConstraint] = []
                
                # Process first factor
                lf_left, cons_left = self.lower_expr(factors[0])
                cons_acc.extend(cons_left)
                current_lf = lf_left
                
                # Process middle factors (if any)
                for f in factors[1:-1]:
                    lf_right, cons_right = self.lower_expr(f)
                    cons_acc.extend(cons_right)
                    t, cons = self.emit_mul(current_lf, lf_right)
                    if cons is not None:
                        cons_acc.append(cons)
                    current_lf = LinearForm.var(t)
                
                # Final factor - assign result to output_var
                final_factor = factors[-1]
                lf_final, cons_final = self.lower_expr(final_factor)
                cons_acc.extend(cons_final)
                
                # Apply numeric coefficient to final factor if needed
                if num % self.p != 1:
                    lf_final = lf_final.scale(num, self.p)
                
                # Create final constraint: current_lf * lf_final + (-output_var) = 0
                final_cons = R1CSConstraint(
                    A=current_lf,
                    B=lf_final,
                    C=LinearForm.var(output_var).scale(-1, self.p)
                )
                cons_acc.append(final_cons)
                return cons_acc
        
        # Fallback
        lf, cons = self.lower_expr(expr)
        final_cons = R1CSConstraint(
            A=lf.scale(-1, self.p),
            B=LinearForm.one(self.p),
            C=LinearForm.var(output_var).scale(-1, self.p)
        )
        cons.append(final_cons)
        return cons
    
    def lower_factored_product_assignment(self, factored_product, rhs_expr) -> List[R1CSConstraint]:
        """
        Handle assignment: factored_product = rhs_expr
        This avoids creating a temporary variable for factored_product by directly
        assigning the RHS to the product result.
        
        For example: (a+b)*(c+d) = x + y becomes:
        - Constraint: (a+b) * (c+d) + (-x-y) = 0
        """
        # Lower the RHS expression to get its linear form
        rhs_lf, rhs_cons = self.lower_expr(rhs_expr)
        
        # Process the factored product
        if isinstance(factored_product, Mul):
            # Separate numeric coefficient and factors
            num = 1
            factors = []
            for arg in factored_product.args:
                if arg.is_Number:
                    num = mulm(num, int(arg), self.p)
                else:
                    factors.append(arg)
            
            if len(factors) == 2:
                # Simple case: two factors
                left_lf, left_cons = self.lower_expr(factors[0])
                right_lf, right_cons = self.lower_expr(factors[1])
                
                # Apply numeric coefficient to the RHS
                if num % self.p != 1:
                    rhs_lf = rhs_lf.scale(num, self.p)
                
                # Create the constraint: left * right + (-rhs) = 0
                # Which is equivalent to: left * right = rhs
                constraint = R1CSConstraint(
                    A=left_lf,
                    B=right_lf,
                    C=rhs_lf.scale(-1, self.p)
                )
                
                all_cons = rhs_cons + left_cons + right_cons + [constraint]
                return all_cons
            else:
                # Multiple factors - fall back to sequential processing
                # This is more complex but handles cases like a*b*c = rhs
                cons_acc = list(rhs_cons)
                
                # Process first factor
                current_lf, cons_left = self.lower_expr(factors[0])
                cons_acc.extend(cons_left)
                
                # Process middle factors
                for f in factors[1:-1]:
                    lf_right, cons_right = self.lower_expr(f)
                    cons_acc.extend(cons_right)
                    t, cons = self.emit_mul(current_lf, lf_right)
                    if cons is not None:
                        cons_acc.append(cons)
                    current_lf = LinearForm.var(t)
                
                # Final factor - assign result to RHS
                final_factor = factors[-1]
                lf_final, cons_final = self.lower_expr(final_factor)
                cons_acc.extend(cons_final)
                
                # Apply numeric coefficient
                if num % self.p != 1:
                    rhs_lf = rhs_lf.scale(num, self.p)
                
                # Create final constraint: current_lf * lf_final + (-rhs) = 0
                final_cons = R1CSConstraint(
                    A=current_lf,
                    B=lf_final,
                    C=rhs_lf.scale(-1, self.p)
                )
                cons_acc.append(final_cons)
                return cons_acc
        
        # Fallback - shouldn't happen for well-formed factored products
        return self.lower_expr_to_output(factored_product, Symbol("fallback_tmp"))


# -----------------------------
# Helper functions for printing
# -----------------------------

def format_coefficient(c: int, p: int = None) -> int:
    """Format coefficient, showing negative equivalent for large numbers close to p"""
    if p is None:
        return c
    
    # Normalize coefficient to be in range [0, p)
    c_mod = c % p
    
    # If the coefficient is in the upper half of the field, show it as negative
    # This makes large numbers like p-1, p-2, etc. display as -1, -2, etc.
    if c_mod > p // 2:
        return c_mod - p
    
    return c_mod

def format_linear_form(lf: LinearForm, p: int = None) -> str:
    """Format a LinearForm for readable printing"""
    terms = []
    
    def format_coeff(c):
        """Format coefficient, showing -1 instead of p-1"""
        return format_coefficient(c, p)
    
    # Add constant term if non-zero
    if lf.const != 0:
        const_formatted = format_coeff(lf.const)
        terms.append(str(const_formatted))
    
    # Add variable terms
    for sym, coeff in lf.coeffs.items():
        if coeff == 0:
            continue
        coeff_formatted = format_coeff(coeff)
        if coeff_formatted == 1:
            terms.append(str(sym))
        elif coeff_formatted == -1:
            terms.append(f"-{sym}")
        else:
            terms.append(f"{coeff_formatted}*{sym}")
    
    if not terms:
        return "0"
    
    result = terms[0]
    for term in terms[1:]:
        if term.startswith('-'):
            result += f" {term}"
        else:
            result += f" + {term}"
    
    return result


# -----------------------------
# JSON writer (R1CS triplets)
# -----------------------------

def collect_symbols_from_constraints(cons_list: List[R1CSConstraint]) -> Set[Symbol]:
    syms: Set[Symbol] = set()
    for c in cons_list:
        syms.update(c.A.coeffs.keys())
        syms.update(c.B.coeffs.keys())
        syms.update(c.C.coeffs.keys())
    return syms


def build_index_map(meta: CircuitMeta, outs: List[Symbol], ins: List[Symbol],
                    tmp_syms: List[Symbol]) -> Tuple[Dict[Symbol, int], List]:
    """
    Index layout:
      0: constant 1
      1..n_outputs: outs
      next n_inputs: ins
      remaining: tmp_syms (sorted by tmp index order)
    """
    idx2sym: List = [Integer(1)] + outs + ins + tmp_syms
    sym2idx: Dict[Symbol, int] = {}
    for idx, s in enumerate(idx2sym):
        if isinstance(s, Integer):
            continue
        sym2idx[s] = idx
    return sym2idx, idx2sym


def lf_to_sparse_dict(lf: LinearForm, sym2idx: Dict[Symbol, int], p: int) -> Dict[str, str]:
    d: Dict[str, str] = {}
    # constant -> index 0
    if lf.const % p != 0:
        d["0"] = str(modp(lf.const, p))
    for s, c in lf.coeffs.items():
        if c % p == 0:
            continue
        d[str(sym2idx[s])] = str(modp(c, p))
    return d


def write_r1cs_json(output_path: str, meta: CircuitMeta,
                    outs: List[Symbol], ins: List[Symbol],
                    constraints: List[R1CSConstraint], p: int):
    # Gather tmp symbols in order of appearance
    seen_tmp: List[Symbol] = []
    seen_set: Set[Symbol] = set()

    def track(lf: LinearForm):
        for s in lf.coeffs.keys():
            if s not in outs and s not in ins and not s.name.startswith("out") and not s.name.startswith("in"):
                # tmp by our naming (tmpK)
                if s not in seen_set:
                    seen_set.add(s)
                    seen_tmp.append(s)

    for c in constraints:
        track(c.A); track(c.B); track(c.C)

    sym2idx, _ = build_index_map(meta, outs, ins, seen_tmp)

    out_constraints = []
    for c in constraints:
        A = lf_to_sparse_dict(c.A, sym2idx, p)
        B = lf_to_sparse_dict(c.B, sym2idx, p)
        C = lf_to_sparse_dict(c.C, sym2idx, p)
        out_constraints.append([A, B, C])

    n_inputs = meta.n_pub + meta.n_prv
    n_vars = 1 + meta.n_outputs + n_inputs + len(seen_tmp)
    
    # Create map array (sequential indices)
    var_map = list(range(n_vars))
    
    payload = {
        "n8": 32,
        "prime": str(meta.prime),
        "nVars": n_vars,
        "nOutputs": meta.n_outputs,
        "nPubInputs": meta.n_pub,
        "nPrvInputs": meta.n_prv,
        "nLabels": n_vars,
        "nConstraints": len(out_constraints),
        "useCustomGates": False,
        "constraints": out_constraints,
        "map": var_map,
        "customGates": [],
        "customGatesUses": []
    }
    with open(output_path, "w") as f:
        json.dump(payload, f, indent=1)


# -----------------------------
# Orchestration
# -----------------------------

def optimize_to_r1cs(input_path: str, output_path: str,
                     enable_groebner: bool = True) -> None:
    print(f"Loading circuit from {input_path}")
    # Load
    meta, raw_constraints = load_circuit(input_path)
    one, outs, ins, internals = build_symbols(meta)
    ios = outs + ins

    print(f"Circuit metadata: {meta.n_vars} vars, {meta.n_outputs} outputs, {meta.n_pub} pub inputs, {meta.n_prv} prv inputs")
    print(f"Symbols: {len(internals)} internal vars, {len(ios)} I/O vars")

    # Convert original R1CS triplets -> polynomial list (A*B - C == 0)
    polys, idx2sym = constraints_to_polys(raw_constraints, meta, (one, outs, ins, internals))
    print(f"Converted {len(raw_constraints)} R1CS constraints to {len(polys)} polynomial equations")

    # ---- Step 1: eliminate internals iteratively, then Groebner eliminate ----
    print("\n=== STEP 1: Eliminating internal variables ===")
    print(f"Starting with {len(polys)} polynomial equations")
    print(f"Trying to eliminate {len(internals)} internal variables: {[str(v) for v in internals]}")
    
    polys_before = len(polys)
    polys = eliminate_internals_iterative(polys, internals, meta.prime)
    polys_after_iterative = len(polys)
    print(f"After iterative elimination: {polys_after_iterative} equations (eliminated {polys_before - polys_after_iterative})")
    
    if enable_groebner:
        print("Applying Groebner basis elimination...")
        polys = groebner_eliminate(polys, internals, ios, meta.prime, enable=enable_groebner)
        polys_after_groebner = len(polys)
        print(f"After Groebner elimination: {polys_after_groebner} equations (eliminated {polys_after_iterative - polys_after_groebner})")
    else:
        print("Groebner elimination disabled")

    # If elimination removed all equations (degenerate), keep 0 == 0 to avoid empty R1CS
    if not polys:
        print("Warning: All equations eliminated, adding trivial 0 == 0 constraint")
        polys = [Integer(0)]

    print(f"Step 1 result: {len(polys)} polynomial equations remaining")
    for i, poly in enumerate(polys):
        print(f"  P{i}: {poly}")

    # ---- Step 2: factor / collect over GF(p) ----
    print("\n=== STEP 2: Factoring and collecting over GF(p) ===")
    print(f"Factoring {len(polys)} polynomial equations...")
    polys = factoring_stage(polys, outs, ios, meta.prime)
    print(f"Step 2 result: {len(polys)} factored polynomial equations")
    for i, poly in enumerate(polys):
        print(f"  P{i}: {poly}")

    # ---- Step 3: lower to R1CS a*b + c == 0 ----
    print("\n=== STEP 3: Compiling to R1CS format ===")
    print(f"Converting {len(polys)} polynomial equations to R1CS constraints...")
    desired_powers = collect_desired_powers(polys)
    builder = R1CSBuilder(p=meta.prime, outs=outs, ins=ins, expr_cache={}, desired_powers=desired_powers)
    # builder = R1CSBuilder(p=meta.prime, outs=outs, ins=ins)
    all_constraints: List[R1CSConstraint] = []
    for i, e in enumerate(polys):
        print(f"  Converting polynomial {i}: {e}")
        constraints = builder.lower_poly_eq_zero(e)
        print(f"    -> Generated {len(constraints)} R1CS constraints:")
        for j, c in enumerate(constraints):
            print(f"      Constraint {j}: A*B + C = 0")
            print(f"        A = {format_linear_form(c.A, meta.prime)}")
            print(f"        B = {format_linear_form(c.B, meta.prime)}")
            print(f"        C = {format_linear_form(c.C, meta.prime)}")
            # Add interpretation for multiplication constraints
            if len(c.C.coeffs) == 1 and list(c.C.coeffs.values())[0] == meta.prime - 1:
                tmp_var = list(c.C.coeffs.keys())[0]
                a_form = format_linear_form(c.A, meta.prime)
                b_form = format_linear_form(c.B, meta.prime)
                print(f"        -> {tmp_var} = ({a_form}) * ({b_form})")
            elif c.B.coeffs == {symbols('1'): 1} or (len(c.B.coeffs) == 0 and c.B.const == 1):
                # Linear constraint
                a_form = format_linear_form(c.A, meta.prime)
                print(f"        -> Linear constraint: {a_form} = 0")
        all_constraints.extend(constraints)

    print(f"Step 3 result: {len(all_constraints)} total R1CS constraints generated")
    print(f"Temporary variables used: {builder.tmp_counter}")

    # Emit final R1CS JSON
    print(f"\nWriting optimized R1CS to {output_path}")
    write_r1cs_json(output_path, meta, outs, ins, all_constraints, meta.prime)
    print("Optimization complete!")


# -----------------------------
# CLI
# -----------------------------

def main():
    ap = argparse.ArgumentParser(description="Optimize circuit and re-emit R1CS in a*b + c = 0 form.")
    ap.add_argument("input", help="Input R1CS-like JSON (with constraints as A,B,C sparse dicts).")
    ap.add_argument("-o", "--output", default="final_r1cs.json", help="Output R1CS JSON path.")
    ap.add_argument("--no-groebner", action="store_true", help="Disable Groebner elimination stage.")
    args = ap.parse_args()

    optimize_to_r1cs(args.input, args.output, enable_groebner=not args.no_groebner)


if __name__ == "__main__":
    main()
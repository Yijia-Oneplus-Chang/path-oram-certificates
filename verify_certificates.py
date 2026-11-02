#!/usr/bin/env python3
"""Verify the manuscript's finite arithmetic premises and saved tables.

Only the Python standard library is required. All proof checks use exact
integers or fractions and remain active under ``python3 -O``.
Run ``python3 verify_certificates.py`` for read-only verification.
Use ``--write-output PATH`` to generate a table at a new, explicit path.
"""
import argparse
from fractions import Fraction as F
from pathlib import Path
import json
import sys


class CertificateError(Exception):
    """A finite premise or a saved certificate entry failed verification."""


def require(condition, message):
    if not condition:
        raise CertificateError(message)


def ceil_grid(x, D):
    return F((x.numerator * D + x.denominator - 1) // x.denominator, D)


def z4_certificate():
    D, H, u = 2**80, 12, F(7, 5)
    a, r = 1 / (u**4 * (2-u)), F(23, 20)
    p = [F(1)] + [F(0)] * (4*H+3)
    M = F(1)
    rows = []
    for h in range(H+1):
        K = 4*(H-h)+3
        require(len(p) == K+1, f'Z=4, h={h}: incorrect coefficient count')
        b = [sum((p[i]*p[j-i] for i in range(j+1)), F(0)) for j in range(K+1)]
        c = sum((b[j] * sum((F(1, 2**(g+1))*(1-u**(g+j-4))
                    for g in range(4-j)), F(0)) for j in range(4)), F(0))
        require(M < F(111,100), f'Z=4, h={h}: moment upper bound is not below 111/100')
        rows.append({'h':h, 'M_upper':str(M), 'B_upper':str(c),
                     'probability_upper':[str(x) for x in p]})
        if h == H:
            break
        nxt = [ceil_grid(sum((b[j]*(1-F(1,2**(5-j))) for j in range(5)), F(0)),D)]
        for k in range(1,K-3):
            nxt.append(ceil_grid(sum((b[j]*F(1,2**(k+5-j))
                                      for j in range(k+5)), F(0)), D))
        M = ceil_grid(a*M*M+c, D)
        p = nxt
    require(M < F(111,100) < r, 'Z=4: terminal moment bound does not imply M_12 < r')
    require(c < F(57621,100000), 'Z=4: B_12 upper bound is not below 57621/100000')
    margin = r-a*r*r-F(57621,100000)
    require(margin == F(19687,720300000) and margin > 0,
            'Z=4: incorrect or nonpositive all-height closure margin')
    return {'grid':D, 'height':H, 'u':str(u), 'r':str(r),
            'a':str(a), 'margin':str(margin), 'rows':rows}


def z3_certificate():
    D, K, H = 2**48, 32, 14
    masses = [D]+[0]*K
    rows = [masses[:]]
    # G is clipped at 35; its mass at 35 includes the whole tail.
    g = [2**(34-j) for j in range(35)]+[1]
    for h in range(H):
        conv = [sum(masses[i]*masses[j-i] for i in range(max(0,j-K),min(K,j)+1))
                for j in range(2*K+1)]
        numer = [0]*(K+1)
        for j, w in enumerate(conv):
            for k, gmass in enumerate(g):
                out = min(K,max(0,j+k-3))
                numer[out] += w*gmass
        nxt = [0]+[numer[k]//(D*2**35) for k in range(1,K+1)]
        nxt[0] = D-sum(nxt[1:])
        # Recompute with geometric survival probabilities, rather than masses.
        tails = [sum(w*2**(35-max(0,k+3-j)) for j,w in enumerate(conv))
                 for k in range(1,K+1)]+[0]
        check = [0]+[(tails[k-1]-tails[k])//(D*2**35) for k in range(1,K+1)]
        check[0] = D-sum(check[1:])
        require(nxt == check, f'Z=3, h={h+1}: mass and tail formulas disagree')
        require(sum(nxt) == D and min(nxt) >= 0,
                f'Z=3, h={h+1}: masses are negative or do not sum to the grid denominator')
        masses = nxt
        rows.append(masses[:])
    mean = F(sum(k*masses[k] for k in range(K+1)),D)
    require(mean == F(775231701030647,140737488355328),
            'Z=3: terminal mean differs from the fraction stated in the manuscript')
    require(mean > F(11,2), 'Z=3: terminal mean is not above 11/2')
    loss_time, loss_coupling = F(1,8192), F(3,64)
    require(F(11,2)-loss_time-loss_coupling == F(44671,8192) > 5,
            'Z=3: subtracting the two error bounds does not leave the stated surplus')
    return {'grid':D, 'cap':K, 'height':H, 'mean_lower':str(mean), 'rows':rows,
            'time_loss_upper':str(loss_time), 'coupling_loss_upper':str(loss_coupling)}


def compare_tables(stored, recomputed, path='certificate'):
    """Report the first unequal entry, including its field or row index."""
    require(type(stored) is type(recomputed),
            f'{path}: stored type {type(stored).__name__}, '
            f'recomputed type {type(recomputed).__name__}')
    if isinstance(recomputed, dict):
        missing = sorted(set(recomputed) - set(stored))
        extra = sorted(set(stored) - set(recomputed))
        require(not missing and not extra,
                f'{path}: missing fields {missing}, unexpected fields {extra}')
        for key in recomputed:
            compare_tables(stored[key], recomputed[key], f'{path}.{key}')
    elif isinstance(recomputed, list):
        require(len(stored) == len(recomputed),
                f'{path}: stored length {len(stored)}, recomputed length {len(recomputed)}')
        for index, (old, new) in enumerate(zip(stored, recomputed)):
            compare_tables(old, new, f'{path}[{index}]')
    else:
        require(stored == recomputed,
                f'{path}: stored {stored!r}, recomputed {recomputed!r}')


def unique_fields(pairs):
    """Reject duplicate JSON fields rather than silently ignoring an entry."""
    result = {}
    for key, value in pairs:
        require(key not in result, f'duplicate JSON field {key!r}')
        result[key] = value
    return result


def print_results(data):
    print('PASS: Z=4 finite conditions (Lemma 4.21; Appendix B, Algorithm 1).')
    print('M_12 upper =', data['z4']['rows'][-1]['M_upper'], '< 111/100 < 23/20')
    print('B_12 upper =', data['z4']['rows'][-1]['B_upper'], '< 57621/100000')
    print('Closure margin =', data['z4']['margin'], '> 0')
    print('PASS: Z=3 all 15 rows agree under the mass and tail formulas.')
    print('PASS: Z=3 finite local surplus (Lemma 5.5; Appendix B, Algorithm 2).')
    print('Lower-distribution mean =', data['z3']['mean_lower'], '> 11/2')
    print('PASS: lower-bound error margins (Appendix C).')
    print('11/2 - 1/8192 - 3/64 = 44671/8192 > 5')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--certificate', type=Path,
                      help='verify this saved JSON table instead of the bundled table')
    mode.add_argument('--write-output', type=Path,
                      help='generate a JSON table at a new path (refuse to overwrite any file)')
    args = parser.parse_args()
    try:
        if args.write_output is None:
            source = args.certificate or Path(__file__).with_name('certificates.json')
            stored = json.loads(source.read_text(encoding='utf-8'),
                                object_pairs_hook=unique_fields)
        data = {'z4': z4_certificate(), 'z3': z3_certificate()}
        if args.write_output is None:
            compare_tables(stored, data)
            print('PASS: saved certificate matches recomputation in every field.')
        else:
            # Exclusive creation also prevents overwriting a symlink target.
            with args.write_output.open('x', encoding='utf-8', newline='\n') as output:
                output.write(json.dumps(data, indent=2)+'\n')
            print('WROTE:', args.write_output)
        print_results(data)
    except (CertificateError, OSError, UnicodeError, json.JSONDecodeError) as error:
        print(f'FAIL: {error}', file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())

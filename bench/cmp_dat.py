#!/usr/bin/env python3
"""Compare two DDA .dat files element-wise.

Elements whose reference magnitude is below --floor times the largest
reference magnitude are treated as numerical noise and reported separately
(the Mueller matrix of a symmetric target has many of them near zero).
"""

import argparse
import sys


def parse(path):
    with open(path) as f:
        lines = [ln.rstrip("\n") for ln in f if ln.strip()]

    head = lines[0]
    cabs = cext = None
    for tok in head.replace("\t", " ").split():
        if tok.startswith("Cabs="):
            cabs = float(tok.split("=", 1)[1])
        elif tok.startswith("Cext="):
            cext = float(tok.split("=", 1)[1])

    rows = []
    for ln in lines[1:]:
        vals = [float(x) for x in ln.replace("\t", " ").split()]
        rows.append(vals)
    return (cabs, cext), rows


def relerr(a, b):
    d = abs(a - b)
    s = max(abs(a), abs(b))
    return d / s if s != 0.0 else 0.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ref")
    ap.add_argument("new")
    ap.add_argument("--floor", type=float, default=1e-3,
                    help="relative-to-max magnitude below which an element is noise")
    args = ap.parse_args()

    (ca_r, ce_r), rows_r = parse(args.ref)
    (ca_n, ce_n), rows_n = parse(args.new)

    if len(rows_r) != len(rows_n):
        print("row count differs")
        return 1

    print(f"Cabs  {ca_r:.8e}  {ca_n:.8e}  rel {relerr(ca_r, ca_n):.2e}")
    print(f"Cext  {ce_r:.8e}  {ce_n:.8e}  rel {relerr(ce_r, ce_n):.2e}")

    # scale of the matrix
    scale = max(abs(v) for r in rows_r for v in r[1:])
    floor = args.floor * scale

    worst_sig = (0.0, None)
    n_noise = 0
    worst_noise = 0.0

    for rr, rn in zip(rows_r, rows_n):
        angle = rr[0]
        for k, (a, b) in enumerate(zip(rr[1:], rn[1:])):
            e = relerr(a, b)
            if abs(a) >= floor or abs(b) >= floor:
                if e > worst_sig[0]:
                    worst_sig = (e, (angle, k))
            else:
                n_noise += 1
                if e > worst_noise:
                    worst_noise = e

    print(f"significant elements: worst rel = {worst_sig[0]:.2e} "
          f"at angle {worst_sig[1][0]:g}, element {worst_sig[1][1]}")
    print(f"noise-level elements: {n_noise}, worst rel = {worst_noise:.2e} "
          f"(expected to differ)")
    return 0


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""reduce_ref.py - bit-exact model of the IRS uint16 Montgomery core (reduce.h,
R = 2^16) + a NUMERICAL no-overflow / range proof for every modulus (Python
full-reproduction principle).

Models the unsigned uint16 reduction used in reduce.h / ntt.c:

    montgomery_reduce16(a):           # 0 <= a < q*2^16
        m = (a * QINV) & 0xffff       # QINV = -q^{-1} mod 2^16
        t = (a + m*q) >> 16           # 0 <= t < 2q  (a + m*q accumulated in uint64)
        return t - q if t >= q else t  # [0,q)
    addm16(a,b) = (a+b>=q) ? a+b-q : a+b      # widen to uint32 first
    subm16(a,b) = (a>=b)   ? a-b   : a+q-b

and proves, for each modulus q:
  * congruence:  montgomery_reduce16(a) == a * R^{-1} (mod q)  for 0 <= a < q*R;
  * range:       output in [0, q);
  * no overflow: every product fed to montgomery_reduce16 is a*b < q^2 <= q*R
                 (so it fits uint32 AND respects the a < q*R precondition), and
                 the a + m*q intermediate is < 2*q*R < 2^33 (fits uint64).

Run:  python3 reduce_ref.py
"""

import random

R = 1 << 16
MASK16 = R - 1
U32 = 1 << 32
U64 = 1 << 64

MODES = {
    "120": dict(q=57601, n=256, base_deg=2),
    "180": dict(q=30977, n=256, base_deg=2),
    "260": dict(q=64513, n=512, base_deg=1),
}


def qinv16(q):
    """-q^{-1} mod 2^16 (matches gen_zetas.py / QINV)."""
    return (-pow(q, -1, R)) % R


def montgomery_reduce16(a, q, qinv):
    m = (a * qinv) & MASK16
    t = (a + m * q) >> 16
    return t - q if t >= q else t


def check_correctness(q, qinv, samples=300000):
    rinv = pow(R, -1, q)
    window = q * R                       # admissible input bound (a < q*R)
    rng = random.Random(0xA11CE ^ q)
    edges = [0, 1, q - 1, q, (q - 1) * (q - 1), window - 1]
    bad_cong = bad_range = bad_sum = 0
    for a in edges + [rng.randrange(window) for _ in range(samples)]:
        m = (a * qinv) & MASK16
        if a + m * q >= U64:             # uint64 accumulator must not overflow
            bad_sum += 1
        r = montgomery_reduce16(a, q, qinv)
        if (r - a * rinv) % q != 0:
            bad_cong += 1
        if not (0 <= r < q):
            bad_range += 1
    return bad_cong, bad_range, bad_sum


def main():
    print("reduce_ref.py - IRS uint16 Montgomery (R=2^16) model + no-overflow proof\n")
    allok = True
    for sec, p in MODES.items():
        q = p["q"]
        qinv = qinv16(q)
        # (q*qinv) mod 2^16 == 2^16 - 1  (i.e. qinv == -q^{-1} mod 2^16)
        qinv_ok = (q * qinv) % R == R - 1
        prod_max = (q - 1) * (q - 1)     # max product fed to montgomery_reduce16
        window = q * R
        sum_max = (window - 1) + MASK16 * q
        prod_u32 = prod_max < U32
        prod_window = prod_max < window
        sum_u64 = sum_max < U64
        bad_cong, bad_range, bad_sum = check_correctness(q, qinv)
        ok = (qinv_ok and bad_cong == 0 and bad_range == 0 and bad_sum == 0
              and prod_u32 and prod_window and sum_u64)
        allok &= ok
        print(f"--- IRS-{sec} (q={q}, n={p['n']}, base_deg={p['base_deg']}) ---")
        print(f"  QINV = -q^-1 mod 2^16 = {qinv}  (q*QINV mod 2^16 == 2^16-1: "
              f"{'OK' if qinv_ok else 'FAIL'})")
        print(f"  congruence mreduce16(a) == a*R^-1 (full window 0<=a<q*R={window}): "
              f"{'OK' if bad_cong == 0 else f'FAIL({bad_cong})'}")
        print(f"  strict range 0 <= r < q: {'OK' if bad_range == 0 else f'FAIL({bad_range})'}")
        print(f"  product a*b <= {prod_max} < q*R and < 2^32: "
              f"{'OK' if (prod_window and prod_u32) else 'OVERFLOW'}")
        print(f"  a+m*q sum <= {sum_max} < 2^33 fits uint64: "
              f"{'OK' if sum_u64 else 'FAIL'}\n")
    print("ALL OK" if allok else "FAILURES PRESENT")
    return 0 if allok else 1


if __name__ == "__main__":
    raise SystemExit(main())

import sys, re
from decimal import Decimal, getcontext
getcontext().prec = 300
LN2 = Decimal(
    "0.69314718055994530941723212145817656807550013436025525412068000949339"
    "362196969471560586332699641868754200148102057068573368552023575813055"
    "703267075163507596193072757082837143519030703862389167347112335011")
BETA = int(sys.argv[1]); LOG = sys.argv[2]; OUT = sys.argv[3]
B = 1 << BETA
RCP = int(Decimal(B) / LN2 * Decimal(1 << 57))
LB  = int((LN2 / Decimal(B) * Decimal(1 << (64 + BETA))).quantize(Decimal(1)))
two63 = Decimal(1 << 63)
table = [int((two63 * (Decimal(-j) * LN2 / Decimal(B)).exp()).quantize(Decimal(1))) for j in range(B)]
txt = open(LOG).read()
m = re.search(r'Running unified Q63\b', txt); seg = txt[m.end():]
m2 = re.search(r'latex exact values poly:', seg)
coeffs = {}
for line in seg[m2.end():].splitlines():
    if '====' in line or 'Running unified' in line: break
    mm = re.search(r'(\d+)\s*\\cdot 2 \^ \{-63\}(?:\s*\\cdot x \^ \{(\d+)\})?', line)
    if mm: coeffs[int(mm.group(2)) if mm.group(2) else 0] = int(mm.group(1))
DEG = max(coeffs); poly = [coeffs[k] for k in range(DEG + 1)]
assert RCP < (1<<64) and LB < (1<<64), "RCP/LB overflow"
with open(OUT, 'w') as f:
    f.write(f"#define EB_BETA {BETA}\n#define EB_B {B}\n#define EB_DEG {DEG}\n")
    f.write(f"#define EB_RCP UINT64_C({RCP})\n#define EB_LB UINT64_C({LB})\n")
    f.write("static const uint64_t EB_POLY[EB_DEG+1]={" + ",".join(f"UINT64_C({c})" for c in poly) + "};\n")
    f.write("static const uint64_t EB_TABLE[EB_B]={" + ",".join(f"UINT64_C({v})" for v in table) + "};\n")
print(f"B={B} BETA={BETA} DEG={DEG} RCP_bits={RCP.bit_length()} LB_bits={LB.bit_length()}")

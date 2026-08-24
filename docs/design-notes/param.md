# IRS Parameter Notes

The public implementation supports Lithium-120, Lithium-180, and Lithium-260. The canonical machine-readable parameters are in `ref/params.h`; generated threshold constants are maintained by `tools/gen_bounds.py` and regenerated through `make tables`.

| Mode | n | q | ell | m | kappa | r | NTT shape |
|---|---:|---:|---:|---:|---:|---:|---|
| Lithium-120 | 256 | 57601 | 2 | 2 | 30 | 182 | incomplete, base degree 2 |
| Lithium-180 | 256 | 30977 | 3 | 3 | 58 | 204 | incomplete, base degree 2 |
| Lithium-260 | 512 | 64513 | 2 | 2 | 58 | 257 | complete, base degree 1 |

The target profile is the no-abort parameterization used by this implementation, not the paper's example parameter table. Signature reserves and rANS tables are derived from the current implementation PMFs in `tools/rans_model.py`, `tools/gen_rans_tables.py`, and `tools/SigSize.py`.

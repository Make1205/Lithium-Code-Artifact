# Rejection Truncation Notes

The iterative rejection sampler uses the exp-recurrence three-region method with `LITHIUM_REJ_T = 33`. The implementation validates the conditional-output relative-error target with `tools/irs_ref.py` and keeps keygen/sign/verify norm thresholds in generated `@@AUTOGEN@@` regions.

The keygen norm window `B_k' <= ||s|| < B_k` supplies the lower bound used by the rejection analysis. The committed thresholds in `ref/params.h` are the source consumed by the C implementation.

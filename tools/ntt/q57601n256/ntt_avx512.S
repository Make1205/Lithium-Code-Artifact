/* Auto-generated AVX-512BW full-ZMM NTT.  Edit gen_ntt.py. */
#include "../config.h"
.text
.global LITHIUM_NAMESPACE(ntt_avx512)
.global LITHIUM_NAMESPACE(invntt_tomont_avx512)
.global LITHIUM_NAMESPACE(basemul_avx512)

.macro load_consts512 qptr
  vmovdqa64 0(\qptr),%zmm16
  vmovdqa64 64(\qptr),%zmm17
.endm
.macro load_tw512 base,disp
  vmovdqa64 \disp(\base),%zmm18
  vmovdqa64 \disp+64(\base),%zmm19
.endm
.macro montmul512 out,in
  vpmullw  %zmm18,\in,%zmm20
  vpmulhuw %zmm19,\in,\out
  vpmulhuw %zmm16,%zmm20,%zmm20
  vpcmpuw  $1,%zmm20,\out,%k1
  vpsubw   %zmm20,\out,\out
  vpsubw   %zmm17,\out,\out{%k1}
.endm
.macro montmul512_pair out,in,zl,zh
  vpmullw  \zl,\in,%zmm20
  vpmulhuw \zh,\in,\out
  vpmulhuw %zmm16,%zmm20,%zmm20
  vpcmpuw  $1,%zmm20,\out,%k1
  vpsubw   %zmm20,\out,\out
  vpsubw   %zmm17,\out,\out{%k1}
.endm
.macro addmod512 dst,a,t
  vpsubw  \t,%zmm16,%zmm21
  vpcmpuw $1,%zmm21,\a,%k1
  vpsubw  %zmm21,\a,\dst
  vpsubw  %zmm17,\dst,\dst{%k1}
.endm
.macro submod512 dst,a,t
  vpcmpuw $1,\t,\a,%k1
  vpsubw  \t,\a,\dst
  vpsubw  %zmm17,\dst,\dst{%k1}
.endm
.macro ctbf512 a,b
  montmul512 %zmm22,\b
  submod512  \b,\a,%zmm22
  addmod512  \a,\a,%zmm22
.endm
.macro gsbf512 a,b
  submod512  %zmm22,\a,\b
  addmod512  \a,\a,\b
  montmul512 \b,%zmm22
.endm
.macro ctbf512_lazy a,b
  montmul512 %zmm22,\b
  vpsubw  %zmm22,%zmm16,%zmm21
  vpaddw  %zmm21,\a,\b
  vpaddw  %zmm22,\a,\a
.endm
.macro ctbf512_csub a,b
  montmul512 %zmm22,\b
  csub512 \a
  vpsubw  %zmm22,%zmm16,%zmm21
  vpaddw  %zmm21,\a,\b
  csub512 \b
  vpaddw  %zmm22,\a,\a
  csub512 \a
.endm
.macro gsbf512_lazy a,b
  vpsubw  \b,%zmm16,%zmm22
  vpaddw  \a,%zmm22,%zmm22
  vpaddw  \b,\a,\a
  montmul512 \b,%zmm22
.endm
.macro gsbf512_csub a,b
  csub512 \a
  csub512 \b
  vpsubw  \b,%zmm16,%zmm22
  vpaddw  \a,%zmm22,%zmm22
  vpaddw  \b,\a,\a
  csub512 \a
  montmul512 \b,%zmm22
.endm

.macro dmul512 out,x,y
  vpmullw  \x,\y,%zmm20
  vpmullw  %zmm28,%zmm20,%zmm20
  vpmulhuw \x,\y,\out
  vpmulhuw %zmm16,%zmm20,%zmm20
  vpcmpuw  $1,%zmm20,\out,%k1
  vpsubw   %zmm20,\out,\out
  vpsubw   %zmm17,\out,\out{%k1}
.endm
.macro shuf16_512 r0,r1
  vmovdqa64 \r0,%zmm31
  vmovdqa64 LITHIUM_NAMESPACE(shuf16_c)(%rip),%zmm30
  vpermt2w \r1,%zmm30,%zmm31
  vmovdqa64 LITHIUM_NAMESPACE(shuf16_d)(%rip),%zmm30
  vpermt2w \r1,%zmm30,\r0
  vmovdqa64 \r0,\r1
  vmovdqa64 %zmm31,\r0
.endm
.macro shuf8_512 r0,r1
  vmovdqa64 \r0,%zmm31
  vmovdqa64 LITHIUM_NAMESPACE(shuf8_c)(%rip),%zmm30
  vpermt2w \r1,%zmm30,%zmm31
  vmovdqa64 LITHIUM_NAMESPACE(shuf8_d)(%rip),%zmm30
  vpermt2w \r1,%zmm30,\r0
  vmovdqa64 \r0,\r1
  vmovdqa64 %zmm31,\r0
.endm
.macro shuf4_512 r0,r1
  vmovdqa64 \r0,%zmm31
  vmovdqa64 LITHIUM_NAMESPACE(shuf4_c)(%rip),%zmm30
  vpermt2w \r1,%zmm30,%zmm31
  vmovdqa64 LITHIUM_NAMESPACE(shuf4_d)(%rip),%zmm30
  vpermt2w \r1,%zmm30,\r0
  vmovdqa64 \r0,\r1
  vmovdqa64 %zmm31,\r0
.endm
.macro shuf2_512 r0,r1
  vmovdqa64 \r0,%zmm31
  vmovdqa64 LITHIUM_NAMESPACE(shuf2_c)(%rip),%zmm30
  vpermt2w \r1,%zmm30,%zmm31
  vmovdqa64 LITHIUM_NAMESPACE(shuf2_d)(%rip),%zmm30
  vpermt2w \r1,%zmm30,\r0
  vmovdqa64 \r0,\r1
  vmovdqa64 %zmm31,\r0
.endm
.macro shuf1_512 r0,r1
  vmovdqa64 \r0,%zmm31
  vmovdqa64 LITHIUM_NAMESPACE(shuf1_c)(%rip),%zmm30
  vpermt2w \r1,%zmm30,%zmm31
  vmovdqa64 LITHIUM_NAMESPACE(shuf1_d)(%rip),%zmm30
  vpermt2w \r1,%zmm30,\r0
  vmovdqa64 \r0,\r1
  vmovdqa64 %zmm31,\r0
.endm


LITHIUM_NAMESPACE(ntt_avx512):
  load_consts512 %rsi
  vmovdqu16 0(%rdi),%zmm0
  vmovdqu16 64(%rdi),%zmm1
  vmovdqu16 128(%rdi),%zmm2
  vmovdqu16 192(%rdi),%zmm3
  vmovdqu16 256(%rdi),%zmm4
  vmovdqu16 320(%rdi),%zmm5
  vmovdqu16 384(%rdi),%zmm6
  vmovdqu16 448(%rdi),%zmm7
  vmovdqa64 0(%rdx),%zmm18
  vmovdqa64 64(%rdx),%zmm19
  ctbf512 %zmm0,%zmm4
  vmovdqa64 128(%rdx),%zmm18
  vmovdqa64 192(%rdx),%zmm19
  ctbf512 %zmm1,%zmm5
  vmovdqa64 256(%rdx),%zmm18
  vmovdqa64 320(%rdx),%zmm19
  ctbf512 %zmm2,%zmm6
  vmovdqa64 384(%rdx),%zmm18
  vmovdqa64 448(%rdx),%zmm19
  ctbf512 %zmm3,%zmm7
  vmovdqa64 512(%rdx),%zmm18
  vmovdqa64 576(%rdx),%zmm19
  ctbf512 %zmm0,%zmm2
  vmovdqa64 640(%rdx),%zmm18
  vmovdqa64 704(%rdx),%zmm19
  ctbf512 %zmm1,%zmm3
  vmovdqa64 768(%rdx),%zmm18
  vmovdqa64 832(%rdx),%zmm19
  ctbf512 %zmm4,%zmm6
  vmovdqa64 896(%rdx),%zmm18
  vmovdqa64 960(%rdx),%zmm19
  ctbf512 %zmm5,%zmm7
  vmovdqa64 1024(%rdx),%zmm18
  vmovdqa64 1088(%rdx),%zmm19
  ctbf512 %zmm0,%zmm1
  vmovdqa64 1152(%rdx),%zmm18
  vmovdqa64 1216(%rdx),%zmm19
  ctbf512 %zmm2,%zmm3
  vmovdqa64 1280(%rdx),%zmm18
  vmovdqa64 1344(%rdx),%zmm19
  ctbf512 %zmm4,%zmm5
  vmovdqa64 1408(%rdx),%zmm18
  vmovdqa64 1472(%rdx),%zmm19
  ctbf512 %zmm6,%zmm7
  shuf16_512 %zmm0,%zmm1
  vmovdqa64 1536(%rdx),%zmm18
  vmovdqa64 1600(%rdx),%zmm19
  ctbf512 %zmm0,%zmm1
  shuf16_512 %zmm2,%zmm3
  vmovdqa64 1664(%rdx),%zmm18
  vmovdqa64 1728(%rdx),%zmm19
  ctbf512 %zmm2,%zmm3
  shuf16_512 %zmm4,%zmm5
  vmovdqa64 1792(%rdx),%zmm18
  vmovdqa64 1856(%rdx),%zmm19
  ctbf512 %zmm4,%zmm5
  shuf16_512 %zmm6,%zmm7
  vmovdqa64 1920(%rdx),%zmm18
  vmovdqa64 1984(%rdx),%zmm19
  ctbf512 %zmm6,%zmm7
  shuf8_512 %zmm0,%zmm1
  vmovdqa64 2048(%rdx),%zmm18
  vmovdqa64 2112(%rdx),%zmm19
  ctbf512 %zmm0,%zmm1
  shuf8_512 %zmm2,%zmm3
  vmovdqa64 2176(%rdx),%zmm18
  vmovdqa64 2240(%rdx),%zmm19
  ctbf512 %zmm2,%zmm3
  shuf8_512 %zmm4,%zmm5
  vmovdqa64 2304(%rdx),%zmm18
  vmovdqa64 2368(%rdx),%zmm19
  ctbf512 %zmm4,%zmm5
  shuf8_512 %zmm6,%zmm7
  vmovdqa64 2432(%rdx),%zmm18
  vmovdqa64 2496(%rdx),%zmm19
  ctbf512 %zmm6,%zmm7
  shuf4_512 %zmm0,%zmm1
  vmovdqa64 2560(%rdx),%zmm18
  vmovdqa64 2624(%rdx),%zmm19
  ctbf512 %zmm0,%zmm1
  shuf4_512 %zmm2,%zmm3
  vmovdqa64 2688(%rdx),%zmm18
  vmovdqa64 2752(%rdx),%zmm19
  ctbf512 %zmm2,%zmm3
  shuf4_512 %zmm4,%zmm5
  vmovdqa64 2816(%rdx),%zmm18
  vmovdqa64 2880(%rdx),%zmm19
  ctbf512 %zmm4,%zmm5
  shuf4_512 %zmm6,%zmm7
  vmovdqa64 2944(%rdx),%zmm18
  vmovdqa64 3008(%rdx),%zmm19
  ctbf512 %zmm6,%zmm7
  shuf2_512 %zmm0,%zmm1
  vmovdqa64 3072(%rdx),%zmm18
  vmovdqa64 3136(%rdx),%zmm19
  ctbf512 %zmm0,%zmm1
  shuf2_512 %zmm2,%zmm3
  vmovdqa64 3200(%rdx),%zmm18
  vmovdqa64 3264(%rdx),%zmm19
  ctbf512 %zmm2,%zmm3
  shuf2_512 %zmm4,%zmm5
  vmovdqa64 3328(%rdx),%zmm18
  vmovdqa64 3392(%rdx),%zmm19
  ctbf512 %zmm4,%zmm5
  shuf2_512 %zmm6,%zmm7
  vmovdqa64 3456(%rdx),%zmm18
  vmovdqa64 3520(%rdx),%zmm19
  ctbf512 %zmm6,%zmm7
  vmovdqu16 %zmm0,0(%rdi)
  vmovdqu16 %zmm1,64(%rdi)
  vmovdqu16 %zmm2,128(%rdi)
  vmovdqu16 %zmm3,192(%rdi)
  vmovdqu16 %zmm4,256(%rdi)
  vmovdqu16 %zmm5,320(%rdi)
  vmovdqu16 %zmm6,384(%rdi)
  vmovdqu16 %zmm7,448(%rdi)
  vzeroupper
  ret

LITHIUM_NAMESPACE(invntt_tomont_avx512):
  load_consts512 %rsi
  vmovdqu16 0(%rdi),%zmm0
  vmovdqu16 64(%rdi),%zmm1
  vmovdqu16 128(%rdi),%zmm2
  vmovdqu16 192(%rdi),%zmm3
  vmovdqu16 256(%rdi),%zmm4
  vmovdqu16 320(%rdi),%zmm5
  vmovdqu16 384(%rdi),%zmm6
  vmovdqu16 448(%rdi),%zmm7
  vmovdqa64 0(%rdx),%zmm18
  vmovdqa64 64(%rdx),%zmm19
  gsbf512 %zmm6,%zmm7
  shuf2_512 %zmm6,%zmm7
  vmovdqa64 128(%rdx),%zmm18
  vmovdqa64 192(%rdx),%zmm19
  gsbf512 %zmm4,%zmm5
  shuf2_512 %zmm4,%zmm5
  vmovdqa64 256(%rdx),%zmm18
  vmovdqa64 320(%rdx),%zmm19
  gsbf512 %zmm2,%zmm3
  shuf2_512 %zmm2,%zmm3
  vmovdqa64 384(%rdx),%zmm18
  vmovdqa64 448(%rdx),%zmm19
  gsbf512 %zmm0,%zmm1
  shuf2_512 %zmm0,%zmm1
  vmovdqa64 512(%rdx),%zmm18
  vmovdqa64 576(%rdx),%zmm19
  gsbf512 %zmm6,%zmm7
  shuf4_512 %zmm6,%zmm7
  vmovdqa64 640(%rdx),%zmm18
  vmovdqa64 704(%rdx),%zmm19
  gsbf512 %zmm4,%zmm5
  shuf4_512 %zmm4,%zmm5
  vmovdqa64 768(%rdx),%zmm18
  vmovdqa64 832(%rdx),%zmm19
  gsbf512 %zmm2,%zmm3
  shuf4_512 %zmm2,%zmm3
  vmovdqa64 896(%rdx),%zmm18
  vmovdqa64 960(%rdx),%zmm19
  gsbf512 %zmm0,%zmm1
  shuf4_512 %zmm0,%zmm1
  vmovdqa64 1024(%rdx),%zmm18
  vmovdqa64 1088(%rdx),%zmm19
  gsbf512 %zmm6,%zmm7
  shuf8_512 %zmm6,%zmm7
  vmovdqa64 1152(%rdx),%zmm18
  vmovdqa64 1216(%rdx),%zmm19
  gsbf512 %zmm4,%zmm5
  shuf8_512 %zmm4,%zmm5
  vmovdqa64 1280(%rdx),%zmm18
  vmovdqa64 1344(%rdx),%zmm19
  gsbf512 %zmm2,%zmm3
  shuf8_512 %zmm2,%zmm3
  vmovdqa64 1408(%rdx),%zmm18
  vmovdqa64 1472(%rdx),%zmm19
  gsbf512 %zmm0,%zmm1
  shuf8_512 %zmm0,%zmm1
  vmovdqa64 1536(%rdx),%zmm18
  vmovdqa64 1600(%rdx),%zmm19
  gsbf512 %zmm6,%zmm7
  shuf16_512 %zmm6,%zmm7
  vmovdqa64 1664(%rdx),%zmm18
  vmovdqa64 1728(%rdx),%zmm19
  gsbf512 %zmm4,%zmm5
  shuf16_512 %zmm4,%zmm5
  vmovdqa64 1792(%rdx),%zmm18
  vmovdqa64 1856(%rdx),%zmm19
  gsbf512 %zmm2,%zmm3
  shuf16_512 %zmm2,%zmm3
  vmovdqa64 1920(%rdx),%zmm18
  vmovdqa64 1984(%rdx),%zmm19
  gsbf512 %zmm0,%zmm1
  shuf16_512 %zmm0,%zmm1
  vmovdqa64 2048(%rdx),%zmm18
  vmovdqa64 2112(%rdx),%zmm19
  gsbf512 %zmm6,%zmm7
  vmovdqa64 2176(%rdx),%zmm18
  vmovdqa64 2240(%rdx),%zmm19
  gsbf512 %zmm4,%zmm5
  vmovdqa64 2304(%rdx),%zmm18
  vmovdqa64 2368(%rdx),%zmm19
  gsbf512 %zmm2,%zmm3
  vmovdqa64 2432(%rdx),%zmm18
  vmovdqa64 2496(%rdx),%zmm19
  gsbf512 %zmm0,%zmm1
  vmovdqa64 2560(%rdx),%zmm18
  vmovdqa64 2624(%rdx),%zmm19
  gsbf512 %zmm5,%zmm7
  vmovdqa64 2688(%rdx),%zmm18
  vmovdqa64 2752(%rdx),%zmm19
  gsbf512 %zmm4,%zmm6
  vmovdqa64 2816(%rdx),%zmm18
  vmovdqa64 2880(%rdx),%zmm19
  gsbf512 %zmm1,%zmm3
  vmovdqa64 2944(%rdx),%zmm18
  vmovdqa64 3008(%rdx),%zmm19
  gsbf512 %zmm0,%zmm2
  vmovdqa64 3072(%rdx),%zmm18
  vmovdqa64 3136(%rdx),%zmm19
  gsbf512 %zmm3,%zmm7
  vmovdqa64 3200(%rdx),%zmm18
  vmovdqa64 3264(%rdx),%zmm19
  gsbf512 %zmm2,%zmm6
  vmovdqa64 3328(%rdx),%zmm18
  vmovdqa64 3392(%rdx),%zmm19
  gsbf512 %zmm1,%zmm5
  vmovdqa64 3456(%rdx),%zmm18
  vmovdqa64 3520(%rdx),%zmm19
  gsbf512 %zmm0,%zmm4
  vmovdqa64 0(%rcx),%zmm18
  vmovdqa64 64(%rcx),%zmm19
  montmul512 %zmm0,%zmm0
  montmul512 %zmm1,%zmm1
  montmul512 %zmm2,%zmm2
  montmul512 %zmm3,%zmm3
  montmul512 %zmm4,%zmm4
  montmul512 %zmm5,%zmm5
  montmul512 %zmm6,%zmm6
  montmul512 %zmm7,%zmm7
  vmovdqu16 %zmm0,0(%rdi)
  vmovdqu16 %zmm1,64(%rdi)
  vmovdqu16 %zmm2,128(%rdi)
  vmovdqu16 %zmm3,192(%rdi)
  vmovdqu16 %zmm4,256(%rdi)
  vmovdqu16 %zmm5,320(%rdi)
  vmovdqu16 %zmm6,384(%rdi)
  vmovdqu16 %zmm7,448(%rdi)
  vzeroupper
  ret

LITHIUM_NAMESPACE(basemul_avx512):
  load_consts512 %rsi
  vmovdqa64 LITHIUM_NAMESPACE(ntt512_qinv)(%rip),%zmm28
  mov $0xaaaaaaaa,%eax
  kmovd %eax,%k3
  vmovdqu16 0(%rdx),%zmm0
  vmovdqu16 0(%rcx),%zmm8
  vmovdqu16 64(%rdx),%zmm1
  vmovdqu16 64(%rcx),%zmm9
  vmovdqu16 128(%rdx),%zmm2
  vmovdqu16 128(%rcx),%zmm10
  vmovdqu16 192(%rdx),%zmm3
  vmovdqu16 192(%rcx),%zmm11
  vmovdqu16 256(%rdx),%zmm4
  vmovdqu16 256(%rcx),%zmm12
  vmovdqu16 320(%rdx),%zmm5
  vmovdqu16 320(%rcx),%zmm13
  vmovdqu16 384(%rdx),%zmm6
  vmovdqu16 384(%rcx),%zmm14
  vmovdqu16 448(%rdx),%zmm7
  vmovdqu16 448(%rcx),%zmm15
  vprold $16,%zmm8,%zmm21
  dmul512 %zmm22,%zmm0,%zmm8
  dmul512 %zmm23,%zmm0,%zmm21
  vmovdqa64 0(%r8),%zmm18
  vmovdqa64 64(%r8),%zmm19
  montmul512 %zmm24,%zmm22
  vprold $16,%zmm24,%zmm25
  addmod512 %zmm24,%zmm22,%zmm25
  vprold $16,%zmm23,%zmm25
  addmod512 %zmm23,%zmm23,%zmm25
  vmovdqa64 %zmm24,%zmm26
  vmovdqu16 %zmm23,%zmm26{%k3}
  vmovdqu16 %zmm26,0(%rdi)
  vprold $16,%zmm9,%zmm21
  dmul512 %zmm22,%zmm1,%zmm9
  dmul512 %zmm23,%zmm1,%zmm21
  vmovdqa64 128(%r8),%zmm18
  vmovdqa64 192(%r8),%zmm19
  montmul512 %zmm24,%zmm22
  vprold $16,%zmm24,%zmm25
  addmod512 %zmm24,%zmm22,%zmm25
  vprold $16,%zmm23,%zmm25
  addmod512 %zmm23,%zmm23,%zmm25
  vmovdqa64 %zmm24,%zmm26
  vmovdqu16 %zmm23,%zmm26{%k3}
  vmovdqu16 %zmm26,64(%rdi)
  vprold $16,%zmm10,%zmm21
  dmul512 %zmm22,%zmm2,%zmm10
  dmul512 %zmm23,%zmm2,%zmm21
  vmovdqa64 256(%r8),%zmm18
  vmovdqa64 320(%r8),%zmm19
  montmul512 %zmm24,%zmm22
  vprold $16,%zmm24,%zmm25
  addmod512 %zmm24,%zmm22,%zmm25
  vprold $16,%zmm23,%zmm25
  addmod512 %zmm23,%zmm23,%zmm25
  vmovdqa64 %zmm24,%zmm26
  vmovdqu16 %zmm23,%zmm26{%k3}
  vmovdqu16 %zmm26,128(%rdi)
  vprold $16,%zmm11,%zmm21
  dmul512 %zmm22,%zmm3,%zmm11
  dmul512 %zmm23,%zmm3,%zmm21
  vmovdqa64 384(%r8),%zmm18
  vmovdqa64 448(%r8),%zmm19
  montmul512 %zmm24,%zmm22
  vprold $16,%zmm24,%zmm25
  addmod512 %zmm24,%zmm22,%zmm25
  vprold $16,%zmm23,%zmm25
  addmod512 %zmm23,%zmm23,%zmm25
  vmovdqa64 %zmm24,%zmm26
  vmovdqu16 %zmm23,%zmm26{%k3}
  vmovdqu16 %zmm26,192(%rdi)
  vprold $16,%zmm12,%zmm21
  dmul512 %zmm22,%zmm4,%zmm12
  dmul512 %zmm23,%zmm4,%zmm21
  vmovdqa64 512(%r8),%zmm18
  vmovdqa64 576(%r8),%zmm19
  montmul512 %zmm24,%zmm22
  vprold $16,%zmm24,%zmm25
  addmod512 %zmm24,%zmm22,%zmm25
  vprold $16,%zmm23,%zmm25
  addmod512 %zmm23,%zmm23,%zmm25
  vmovdqa64 %zmm24,%zmm26
  vmovdqu16 %zmm23,%zmm26{%k3}
  vmovdqu16 %zmm26,256(%rdi)
  vprold $16,%zmm13,%zmm21
  dmul512 %zmm22,%zmm5,%zmm13
  dmul512 %zmm23,%zmm5,%zmm21
  vmovdqa64 640(%r8),%zmm18
  vmovdqa64 704(%r8),%zmm19
  montmul512 %zmm24,%zmm22
  vprold $16,%zmm24,%zmm25
  addmod512 %zmm24,%zmm22,%zmm25
  vprold $16,%zmm23,%zmm25
  addmod512 %zmm23,%zmm23,%zmm25
  vmovdqa64 %zmm24,%zmm26
  vmovdqu16 %zmm23,%zmm26{%k3}
  vmovdqu16 %zmm26,320(%rdi)
  vprold $16,%zmm14,%zmm21
  dmul512 %zmm22,%zmm6,%zmm14
  dmul512 %zmm23,%zmm6,%zmm21
  vmovdqa64 768(%r8),%zmm18
  vmovdqa64 832(%r8),%zmm19
  montmul512 %zmm24,%zmm22
  vprold $16,%zmm24,%zmm25
  addmod512 %zmm24,%zmm22,%zmm25
  vprold $16,%zmm23,%zmm25
  addmod512 %zmm23,%zmm23,%zmm25
  vmovdqa64 %zmm24,%zmm26
  vmovdqu16 %zmm23,%zmm26{%k3}
  vmovdqu16 %zmm26,384(%rdi)
  vprold $16,%zmm15,%zmm21
  dmul512 %zmm22,%zmm7,%zmm15
  dmul512 %zmm23,%zmm7,%zmm21
  vmovdqa64 896(%r8),%zmm18
  vmovdqa64 960(%r8),%zmm19
  montmul512 %zmm24,%zmm22
  vprold $16,%zmm24,%zmm25
  addmod512 %zmm24,%zmm22,%zmm25
  vprold $16,%zmm23,%zmm25
  addmod512 %zmm23,%zmm23,%zmm25
  vmovdqa64 %zmm24,%zmm26
  vmovdqu16 %zmm23,%zmm26{%k3}
  vmovdqu16 %zmm26,448(%rdi)
  vzeroupper
  ret


/* nttunpack: standard-order int32 samples -> avx512 NTT layout.
 * Kyber-style forward shuffle ladder, no butterflies; narrows int32 [0,q) -> int16.
 *   rdi = dst (int16[n]), rsi = src (int32[n]). */
.global LITHIUM_NAMESPACE(nttunpack_avx512)
LITHIUM_NAMESPACE(nttunpack_avx512):
  vmovdqu32 0(%rsi),%zmm20
  vpmovdw %zmm20,%ymm0
  vmovdqu32 64(%rsi),%zmm20
  vpmovdw %zmm20,%ymm21
  vinserti64x4 $1,%ymm21,%zmm0,%zmm0
  vmovdqu32 128(%rsi),%zmm20
  vpmovdw %zmm20,%ymm1
  vmovdqu32 192(%rsi),%zmm20
  vpmovdw %zmm20,%ymm21
  vinserti64x4 $1,%ymm21,%zmm1,%zmm1
  vmovdqu32 256(%rsi),%zmm20
  vpmovdw %zmm20,%ymm2
  vmovdqu32 320(%rsi),%zmm20
  vpmovdw %zmm20,%ymm21
  vinserti64x4 $1,%ymm21,%zmm2,%zmm2
  vmovdqu32 384(%rsi),%zmm20
  vpmovdw %zmm20,%ymm3
  vmovdqu32 448(%rsi),%zmm20
  vpmovdw %zmm20,%ymm21
  vinserti64x4 $1,%ymm21,%zmm3,%zmm3
  vmovdqu32 512(%rsi),%zmm20
  vpmovdw %zmm20,%ymm4
  vmovdqu32 576(%rsi),%zmm20
  vpmovdw %zmm20,%ymm21
  vinserti64x4 $1,%ymm21,%zmm4,%zmm4
  vmovdqu32 640(%rsi),%zmm20
  vpmovdw %zmm20,%ymm5
  vmovdqu32 704(%rsi),%zmm20
  vpmovdw %zmm20,%ymm21
  vinserti64x4 $1,%ymm21,%zmm5,%zmm5
  vmovdqu32 768(%rsi),%zmm20
  vpmovdw %zmm20,%ymm6
  vmovdqu32 832(%rsi),%zmm20
  vpmovdw %zmm20,%ymm21
  vinserti64x4 $1,%ymm21,%zmm6,%zmm6
  vmovdqu32 896(%rsi),%zmm20
  vpmovdw %zmm20,%ymm7
  vmovdqu32 960(%rsi),%zmm20
  vpmovdw %zmm20,%ymm21
  vinserti64x4 $1,%ymm21,%zmm7,%zmm7
  shuf16_512 %zmm0,%zmm1
  shuf16_512 %zmm2,%zmm3
  shuf16_512 %zmm4,%zmm5
  shuf16_512 %zmm6,%zmm7
  shuf8_512 %zmm0,%zmm1
  shuf8_512 %zmm2,%zmm3
  shuf8_512 %zmm4,%zmm5
  shuf8_512 %zmm6,%zmm7
  shuf4_512 %zmm0,%zmm1
  shuf4_512 %zmm2,%zmm3
  shuf4_512 %zmm4,%zmm5
  shuf4_512 %zmm6,%zmm7
  shuf2_512 %zmm0,%zmm1
  shuf2_512 %zmm2,%zmm3
  shuf2_512 %zmm4,%zmm5
  shuf2_512 %zmm6,%zmm7
  vmovdqu16 %zmm0,0(%rdi)
  vmovdqu16 %zmm1,64(%rdi)
  vmovdqu16 %zmm2,128(%rdi)
  vmovdqu16 %zmm3,192(%rdi)
  vmovdqu16 %zmm4,256(%rdi)
  vmovdqu16 %zmm5,320(%rdi)
  vmovdqu16 %zmm6,384(%rdi)
  vmovdqu16 %zmm7,448(%rdi)
  vzeroupper
  ret


/* nttunpack16: standard-order uint16 samples -> avx512 NTT layout.
 * Same shuffle ladder as nttunpack_avx512 but src is ALREADY int16
 * (in [0,q)): load directly, no int32 widen / vpmovdw.  In-place safe
 * (whole poly loaded before any store).
 *   rdi = dst (int16[n]), rsi = src (int16[n]). */
.global LITHIUM_NAMESPACE(nttunpack16_avx512)
LITHIUM_NAMESPACE(nttunpack16_avx512):
  vmovdqu16 0(%rsi),%zmm0
  vmovdqu16 64(%rsi),%zmm1
  vmovdqu16 128(%rsi),%zmm2
  vmovdqu16 192(%rsi),%zmm3
  vmovdqu16 256(%rsi),%zmm4
  vmovdqu16 320(%rsi),%zmm5
  vmovdqu16 384(%rsi),%zmm6
  vmovdqu16 448(%rsi),%zmm7
  shuf16_512 %zmm0,%zmm1
  shuf16_512 %zmm2,%zmm3
  shuf16_512 %zmm4,%zmm5
  shuf16_512 %zmm6,%zmm7
  shuf8_512 %zmm0,%zmm1
  shuf8_512 %zmm2,%zmm3
  shuf8_512 %zmm4,%zmm5
  shuf8_512 %zmm6,%zmm7
  shuf4_512 %zmm0,%zmm1
  shuf4_512 %zmm2,%zmm3
  shuf4_512 %zmm4,%zmm5
  shuf4_512 %zmm6,%zmm7
  shuf2_512 %zmm0,%zmm1
  shuf2_512 %zmm2,%zmm3
  shuf2_512 %zmm4,%zmm5
  shuf2_512 %zmm6,%zmm7
  vmovdqu16 %zmm0,0(%rdi)
  vmovdqu16 %zmm1,64(%rdi)
  vmovdqu16 %zmm2,128(%rdi)
  vmovdqu16 %zmm3,192(%rdi)
  vmovdqu16 %zmm4,256(%rdi)
  vmovdqu16 %zmm5,320(%rdi)
  vmovdqu16 %zmm6,384(%rdi)
  vmovdqu16 %zmm7,448(%rdi)
  vzeroupper
  ret

.text

/* ============ NTT input conversion helpers (standard order) ============
 * C ABI:
 *   poly_freeze_to_u16_avx512(dst:int16_t*, src:int32_t*)
 *   poly_freeze_small_to_u16_avx512(dst:int16_t*, src:int32_t*)
 *   poly_cast_modq_avx512(dst:int16_t*, src:int32_t*)
 * These helpers feed ntt_avx512's standard-order input buffer; unlike
 * nttunpack_avx512, they do not apply the AVX-512 NTT slot shuffle ladder. */
.p2align 5
.global LITHIUM_NAMESPACE(poly_freeze_to_u16_avx512)
LITHIUM_NAMESPACE(poly_freeze_to_u16_avx512):
  vpbroadcastd .Lfreeze_q_avx512(%rip),%zmm7
  vpbroadcastd .Lfreeze_rec_lo_avx512(%rip),%zmm5
  vpbroadcastd .Lfreeze_rec_hi_avx512(%rip),%zmm8
  mov $16,%ecx
.Lfreeze_full_loop_avx512:
  vmovdqu32 (%rsi),%zmm0
  vpsrad    $31,%zmm0,%zmm1
  vpxord    %zmm1,%zmm0,%zmm2
  vpsubd    %zmm1,%zmm2,%zmm2
  vpmuludq  %zmm5,%zmm2,%zmm3
  vpsrlq    $47,%zmm3,%zmm3
  vpsrlq    $32,%zmm2,%zmm4
  vpmuludq  %zmm5,%zmm4,%zmm4
  vpsrlq    $47,%zmm4,%zmm4
  vpsllq    $32,%zmm4,%zmm4
  vpord     %zmm4,%zmm3,%zmm3
  vpsrld    $15,%zmm2,%zmm4
  vpmulld   %zmm8,%zmm4,%zmm4
  vpaddd    %zmm4,%zmm3,%zmm3
  vpmulld   %zmm7,%zmm3,%zmm3
  vpsubd    %zmm3,%zmm2,%zmm2
  vpsubd    %zmm7,%zmm2,%zmm2
  vpsrad    $31,%zmm2,%zmm4
  vpandd    %zmm7,%zmm4,%zmm4
  vpaddd    %zmm4,%zmm2,%zmm2
  vpxord    %zmm1,%zmm2,%zmm2
  vpsubd    %zmm1,%zmm2,%zmm2
  vpsrad    $31,%zmm2,%zmm4
  vpandd    %zmm7,%zmm4,%zmm4
  vpaddd    %zmm4,%zmm2,%zmm2
  vpmovdw   %zmm2,%ymm0
  vmovdqu   %ymm0,(%rdi)
  add       $64,%rsi
  add       $32,%rdi
  dec       %ecx
  jnz       .Lfreeze_full_loop_avx512
  ret

.p2align 5
.global LITHIUM_NAMESPACE(poly_freeze_small_to_u16_avx512)
LITHIUM_NAMESPACE(poly_freeze_small_to_u16_avx512):
  vpbroadcastd .Lfreeze_q_avx512(%rip),%zmm7
  mov $16,%ecx
.Lfreeze_small_loop_avx512:
  vmovdqu32 (%rsi),%zmm0
  vpsrad    $31,%zmm0,%zmm1
  vpandd    %zmm7,%zmm1,%zmm1
  vpaddd    %zmm1,%zmm0,%zmm0
  vpmovdw   %zmm0,%ymm0
  vmovdqu   %ymm0,(%rdi)
  add       $64,%rsi
  add       $32,%rdi
  dec       %ecx
  jnz       .Lfreeze_small_loop_avx512
  ret

.p2align 5
.global LITHIUM_NAMESPACE(poly_cast_modq_avx512)
LITHIUM_NAMESPACE(poly_cast_modq_avx512):
  mov $16,%ecx
.Lfreeze_cast_loop_avx512:
  vmovdqu32 (%rsi),%zmm0
  vpmovdw   %zmm0,%ymm0
  vmovdqu   %ymm0,(%rdi)
  add       $64,%rsi
  add       $32,%rdi
  dec       %ecx
  jnz       .Lfreeze_cast_loop_avx512
  ret

.section .rodata
.p2align 4
.Lfreeze_q_avx512:       .long 57601
.Lfreeze_rec_lo_avx512:  .long 2443316754
.Lfreeze_rec_hi_avx512:  .long 0
.text
.section .note.GNU-stack,"",@progbits

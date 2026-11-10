/* Public labels retain the frozen bit packing. Unlab acts on W blocks. */
int copa_unpack_pk(uint8_t rho[32],int32_t *b,const uint8_t *pk,size_t len,unsigned mode){
 if(len!=copa_pk_bytes(mode))return 0;
 uint8_t padded[1480+8];memcpy(padded,pk,len);memset(padded+len,0,8);pk=padded;memcpy(rho,pk,32);
 unsigned n=copa_label_width(mode),shift[16],odd[16],s=sigma(mode),q=copa_modulus(mode),low_total=0;copa_pk_wide power=1;size_t pos=256;
 for(unsigned i=0;i<n;++i){unsigned p=copa_label_radix(mode,i);shift[i]=lowbits(p);odd[i]=p>>shift[i];power*=odd[i];low_total+=shift[i];}
 unsigned mixed_bits=bitlength(power-1);
 for(unsigned block=0;block<64;block+=COPAV_W){uint32_t labels[16][COPAV_W];
  for(unsigned lane=0;lane<COPAV_W;++lane){uint64_t low=(uint64_t)get(pk,&pos,low_total);for(unsigned i=0;i<n;++i){labels[i][lane]=(uint32_t)low&((1u<<shift[i])-1);low>>=shift[i];}copa_pk_wide mixed=get(pk,&pos,mixed_bits);if(mixed>=power)return 0;
   uint32_t row[16];for(unsigned i=0;i<n;++i)row[i]=labels[i][lane];unpack_digits(row,shift,n,odd[0],mixed);for(unsigned i=0;i<n;++i)labels[i][lane]=row[i];}

  CP_BEGIN(punlab);copav_vec rm[16];
  for(unsigned i=0;i<n;++i)rm[i]=cv_mul32(cv_load(labels[i]),cv_set(s*diagonal(mode,i)));
  /* basis[i,j]=[j subset i]*diag[j]: exact subset-zeta transform. */
  for(unsigned bit=1;bit<n;bit<<=1)for(unsigned i=0;i<n;++i)if(i&bit)rm[i]=cv_add(rm[i],rm[i^bit]);
  for(unsigned i=0;i<n;++i)rm[i]=cv_sub(rm[i],cv_mul32(qv_div_positive(rm[i],q),cv_set(q)));
  for(unsigned i=0;i<n;++i){cv_store(b+(i/4)*256+(i%4)*64+block,rm[pi(mode,i)]);}CP_END(CP_UNLAB,punlab);
 }
 if(mode==6)for(unsigned i=0;i<256;++i){unsigned v=(unsigned)get(pk,&pos,bitlength(q/3-1));if(v>=q/3)return 0;b[512+i]=3*(int32_t)v;}
 assert(pos==8*len);return 1;
}

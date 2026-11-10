/* Pack each butterfly's lower/upper inputs from adjacent vectors, so
 * small levels compute one product per butterfly instead of two. */
static inline copav_vec ntt_zip_low(copav_vec a,copav_vec b,unsigned len){
#if COPAV_W==8
 if(len==4)return _mm256_permute2x128_si256(a,b,0x20);
 if(len==2)return _mm256_unpacklo_epi64(a,b);
 return _mm256_castps_si256(_mm256_shuffle_ps(_mm256_castsi256_ps(a),_mm256_castsi256_ps(b),0x88));
#else
 if(len==8)return _mm512_shuffle_i32x4(a,b,0x44);
 if(len==4)return _mm512_shuffle_i32x4(a,b,0x88);
 if(len==2)return _mm512_unpacklo_epi64(a,b);
 return _mm512_castps_si512(_mm512_shuffle_ps(_mm512_castsi512_ps(a),_mm512_castsi512_ps(b),0x88));
#endif
}
static inline copav_vec ntt_zip_high(copav_vec a,copav_vec b,unsigned len){
#if COPAV_W==8
 if(len==4)return _mm256_permute2x128_si256(a,b,0x31);
 if(len==2)return _mm256_unpackhi_epi64(a,b);
 return _mm256_castps_si256(_mm256_shuffle_ps(_mm256_castsi256_ps(a),_mm256_castsi256_ps(b),0xdd));
#else
 if(len==8)return _mm512_shuffle_i32x4(a,b,0xee);
 if(len==4)return _mm512_shuffle_i32x4(a,b,0xdd);
 if(len==2)return _mm512_unpackhi_epi64(a,b);
 return _mm512_castps_si512(_mm512_shuffle_ps(_mm512_castsi512_ps(a),_mm512_castsi512_ps(b),0xdd));
#endif
}
static inline copav_vec ntt_unzip_low(copav_vec a,copav_vec b,unsigned len){
#if COPAV_W==8
 if(len==4)return _mm256_permute2x128_si256(a,b,0x20);
 if(len==2)return _mm256_unpacklo_epi64(a,b);
 return _mm256_unpacklo_epi32(a,b);
#else
 if(len==8)return _mm512_shuffle_i32x4(a,b,0x44);
 if(len==4)return _mm512_permutex2var_epi64(a,_mm512_setr_epi64(0,1,8,9,2,3,10,11),b);
 if(len==2)return _mm512_unpacklo_epi64(a,b);
 return _mm512_unpacklo_epi32(a,b);
#endif
}
static inline copav_vec ntt_unzip_high(copav_vec a,copav_vec b,unsigned len){
#if COPAV_W==8
 if(len==4)return _mm256_permute2x128_si256(a,b,0x31);
 if(len==2)return _mm256_unpackhi_epi64(a,b);
 return _mm256_unpackhi_epi32(a,b);
#else
 if(len==8)return _mm512_shuffle_i32x4(a,b,0xee);
 if(len==4)return _mm512_permutex2var_epi64(a,_mm512_setr_epi64(4,5,12,13,6,7,14,15),b);
 if(len==2)return _mm512_unpackhi_epi64(a,b);
 return _mm512_unpackhi_epi32(a,b);
#endif
}



/****************************************************************************************\
*                                     RGB <-> L*a*b*                                     *
\****************************************************************************************/

/** https://www.image-engineering.de/library/technotes/958-how-to-convert-between-srgb-and-ciexyz **/
#define  labXr_32f  0.4124564f /* = xyzXr_32f */
#define  labXg_32f  0.3575761f /* = xyzXg_32f */
#define  labXb_32f  0.1804375f /* = xyzXb_32f */

#define  labYr_32f  0.212671f /* = xyzYr_32f */
#define  labYg_32f  0.715160f /* = xyzYg_32f */ 
#define  labYb_32f  0.072169f /* = xyzYb_32f */ 

#define  labZr_32f  0.0193339f /* = xyzZr_32f */
#define  labZg_32f  0.1191920f /* = xyzZg_32f */
#define  labZb_32f  0.9503041f /* = xyzZb_32 */

#define  labT_32f   0.008856f

/* https://en.wikipedia.org/wiki/CIELAB_color_space#From_CIEXYZ_to_CIELAB */
#define labXScale_32f 0.01052090029f /* 1/95.0489 (D65 illuminant) */
#define labYScale_32f .01f
#define labZScale_32f 0.009184085816f /* 1/108.8840 (D65 illuminant) */

#define labSmallScale_32f  7.787f
#define labSmallShift_32f  0.13793103448275862f  /* 16/116 */
#define labLScale_32f      116.f
#define labLShift_32f      16.f
#define labLScale2_32f     903.3f

#define rgbT_32f 0.04045f
#define rgbScale_32f 0.0773993808f /* 1/12.92 */
#define rgbLShift_32f 0.055f
#define rgbLScale_32f 0.9478672986f /* 1/1.055 */
#define rgbLPow_32f 2.4f

/** @brief Converts RGB (0-255) values to CIELAB colorspace
 * @param[in] R red channel (0-255)
 * @param[in] G green channel (0-255)
 * @param[in] B blue channel (0-255)
 * @param[out] L lightness channel (0-100)
 * @param[out] a axis (-128 - 128)
 * @param[out] b axis (-128 - 128)
 **/
void rgb2lab(int R, int G, int B, float *L, float *a, float *b);
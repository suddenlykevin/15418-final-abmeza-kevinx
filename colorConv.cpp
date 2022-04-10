#include <math.h>
#include <stdio.h>
#include "colorConv.h"

/** https://www.image-engineering.de/library/technotes/958-how-to-convert-between-srgb-and-ciexyz **/
float linearize(int V) {
    float Vf = ((float) V) / 255.0f;
    if (Vf > rgbT_32f) {
        Vf = (Vf + rgbLShift_32f)*rgbLScale_32f;
        return pow(Vf, rgbLPow_32f) * 100;
    } else {
        return Vf * rgbScale_32f * 100;
    }
}

float RGBtoX(float R, float G, float B) {
    return labXr_32f * R + labXg_32f * G + labZb_32f * B;
}

float RGBtoY(float R, float G, float B) {
    return labYr_32f * R + labYg_32f * G + labYb_32f * B;
}

float RGBtoZ(float R, float G, float B) {
    return labZr_32f * R + labZg_32f * G + labZb_32f * B;
}

float XYZtoL(float X, float Y, float Z) {
    if (Y > labT_32f) {
        Y = cbrt(Y);
        return Y*labLScale_32f - labLShift_32f;
    } else {
        return Y*labSmallScale_32f + labSmallShift_32f;
    }
}

float XYZtoA(float X, float Y, float Z) {
    if (X > labT_32f) {
        X = cbrt(X);
    } else {
        X = X*labSmallScale_32f + labSmallShift_32f;
    }

    if (Y > labT_32f) {
        Y = cbrt(Y);
    } else {
        Y = Y*labSmallScale_32f + labSmallShift_32f;
    }

    return 500.f*(X - Y);
}

float XYZtoB(float X, float Y, float Z) {
    if (Y > labT_32f) {
        Y = cbrt(Y);
    } else {
        Y = Y*labSmallScale_32f + labSmallShift_32f;
    }

    if (Z > labT_32f) {
        Z = cbrt(Z);
    } else {
        Z = Z*labSmallScale_32f + labSmallShift_32f;
    }
    return 200.f*(Y - Z);
}

void rgb2lab(int R, int G, int B, float *L, float *a, float *b) {
    float Rf = linearize(R);
    float Gf = linearize(G);
    float Bf = linearize(B);

    float X = RGBtoX(Rf, Gf, Bf) * labXScale_32f;
    float Y = RGBtoY(Rf, Gf, Bf) * labYScale_32f;
    float Z = RGBtoZ(Rf, Gf, Bf) * labZScale_32f;

    *L = XYZtoL(X, Y, Z);
    *a = XYZtoA(X, Y, Z);
    *b = XYZtoB(X, Y, Z);
}
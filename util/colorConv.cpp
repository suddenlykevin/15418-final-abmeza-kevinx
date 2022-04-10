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

unsigned char delin(float V) {
    V = V/100;
    float Vf;
    if (V > rgbTinv_32f) {
        Vf = rgbInvScale_32f * pow(V, rgbLPowInv_32f) - rgbLShift_32f;
    } else {
        Vf = rgbSmallScale_32f * V; 
    }
    return (unsigned char) (Vf * 255.f);
}

float RGBtoX(float R, float G, float B) {
    return labXr_32f * R + labXg_32f * G + labXb_32f * B;
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

float XYZtoR(float X, float Y, float Z) {
    return labRx_32f * X + labRy_32f * Y + labRz_32f * Z;
}

float XYZtoG(float X, float Y, float Z) {
    return labGx_32f * X + labGy_32f * Y + labGz_32f * Z;

}

float XYZtorgB(float X, float Y, float Z) {
    return labBx_32f * X + labBy_32f * Y + labBz_32f * Z;
}

/* http://www.brucelindbloom.com/index.html?Eqn_Lab_to_XYZ.html */
float LABtoX(float L, float a, float b) {
    float fy = (L + labLShift_32f)/labLScale_32f;
    float fx = a*labAScale_32f + fy;
    if (pow(fx, labPow_32f) > labT_32f) {
        return pow(fx, labPow_32f);
    } else {
        return (labLScale_32f*fx - labLShift_32f)/labLScale2_32f;
    }
}

float LABtoY(float L, float a, float b) {
    float fy = (L + labLShift_32f)/labLScale_32f;
    if (L > labLScale2_32f*labT_32f) {
        return pow(fy, labPow_32f);
    } else {
        return (labLScale_32f*fy - labLShift_32f)/labLScale2_32f;
    }
}

float LABtoZ(float L, float a, float b) {
    float fy = (L + labLShift_32f)/labLScale_32f;
    float fz = fy - b*labBScale_32f;
    if (pow(fz, labPow_32f) > labT_32f) {
        return pow(fz, labPow_32f);
    } else {
        return (labLScale_32f*fz - labLShift_32f)/labLScale2_32f;
    }
}

void rgb2lab(int R, int G, int B, float *L, float *a, float *b) {
    float Rf = linearize(R);
    float Gf = linearize(G);
    float Bf = linearize(B);

    float X = RGBtoX(Rf, Gf, Bf);
    float Y = RGBtoY(Rf, Gf, Bf);
    float Z = RGBtoZ(Rf, Gf, Bf);

    X *= labXScale_32f;
    Y *= labYScale_32f;
    Z *= labZScale_32f;

    *L = XYZtoL(X, Y, Z);
    *a = XYZtoA(X, Y, Z);
    *b = XYZtoB(X, Y, Z);
}

void lab2rgb(float L, float a, float b, unsigned char *R, unsigned char *G, unsigned char *B) {
    float X = LABtoX(L, a, b);
    float Y = LABtoY(L, a, b);
    float Z = LABtoZ(L, a, b);
    
    X *= labXScaleInv_32f;
    Y *= labYScaleInv_32f;
    Z *= labZScaleInv_32f;

    float Rf = XYZtoR(X, Y, Z);
    float Gf = XYZtoG(X, Y, Z);
    float Bf = XYZtorgB(X, Y, Z);

    *R = delin(Rf);
    *G = delin(Gf);
    *B = delin(Bf);
}
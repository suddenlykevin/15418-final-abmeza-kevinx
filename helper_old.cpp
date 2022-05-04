/**
 * @file helper.cpp
 * @author Kevin Xie (kevinx)
 *         Anthony Meza (abmeza)
 * @brief 
 * @version 0.1
 * @date 2022-05-02
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>

#include "helper.h"

//********************************************************//
//*******************  HELPER FUNCTIONS ******************//
//********************************************************//

int maxEigen3(float *matrix, float *value, LabColor *vector) {
    
    // extract unique values from top triangle
    float a = matrix[0];
    float b = matrix[4];
    float c = matrix[8];
    float d = matrix[1];
    float e = matrix[5];
    float f = matrix[3];

    float x_1 = a*a + b*b + c*c - a*b - a*c - b*c + 3*(d*d + f*f + e*e);
    float x_2 = -(2*a - b - c) * (2*b - a - c) * (2*c - a - b) + 
                9*((2*c - a - b)*d*d + (2*b - a - c)*f*f + (2*a - b - c)*e*e) - 54*d*e*f;
    
    float phi;
    if (x_2 > 0) {
        phi = atan(sqrt(4*(x_1*x_1*x_1) - x_2*x_2)/x_2);
    } else if (x_2 == 0) {
        phi = M_PI/2;
    } else {
        phi = atan(sqrt(4*(x_1*x_1*x_1) - x_2*x_2)/x_2) + M_PI;
    }

    float lam_1 = (a + b + c - 2*sqrt(x_1)*cos(phi/3))/3;
    float lam_2 = (a + b + c + 2*sqrt(x_1)*cos((phi - M_PI)/3))/3;
    float lam_3 = (a + b + c + 2*sqrt(x_1)*cos((phi + M_PI)/3))/3;

    float lam = (lam_2 > lam_1) ? lam_2 : ((lam_3 > lam_1) ? lam_3 : lam_1);

    // unlikely special case causes divide by zero
    if (f == 0 || f*(b-lam)-d*e == 0) {
        printf("(f: %f, b: %f, lam: %f, d: %f, e: %f)\n", f, b, lam, d, e);
        return -1;
    }

    // store highest eigenvalue
    * value = lam;

    float m = (d*(c-lam)-e*f)/(f*(b-lam)-d*e);

    float v_0 = (lam - c - e*m)/f;

    // store associated eigenvector
    vector->L = v_0;
    vector->a = m;
    vector->b = 1.f;

    return 0;
}


float dist_k(int m, float S, float l_k, float a_k, float b_k, int x_k, int y_k,
           float l_i, float a_i, float b_i, int x_i, int y_i) {
    float d_lab = sqrt(pow(l_k - l_i, 2.f) + pow(a_k - a_i, 2.f) + pow(b_k - b_i, 2.f));
    float d_xy = sqrt(pow((float) x_k - x_i, 2.f) + pow((float) y_k - y_i, 2.f));
    float k = ((float) m) / S;
    return d_lab + k * d_xy;
}


float gaussian(float x, float sigma, float mean) {
    return exp((x-mean)*(x-mean)/(-2.0f*sigma*sigma))/sqrt(6.28319*sigma*sigma);
}

void* wrp_malloc(size_t size){ 
    void* ptr = malloc(size);

    // Check that no error occured
    if (ptr == NULL) {
        printf("Unable to allocate memory with malloc of size %zu\n", size);
        exit(1);
    }

    return ptr;
}

void* wrp_calloc(size_t nitems, size_t size){ 
    void* ptr = calloc(nitems, size);

    // Check that no error occured
    if (ptr == NULL) {
        printf("Unable to allocate memory with calloc of items %zu and size %zu\n", nitems, size);
        exit(1);
    }

    return ptr;
}



/** https://www.image-engineering.de/library/technotes/958-how-to-convert-between-srgb-and-ciexyz **/
float linearize(int V) {
    float Vf = ((float) V) / 255.0f;
    if (Vf > rgbT_32f) {
        Vf = (Vf + rgbLShift_32f)*rgbLScale_32f;
        return pow(Vf, rgbLPow_32f) * 100;
    } else {
        return (Vf * rgbScale_32f) * 100;
    }
}

unsigned char delin(float V) {
    V = V / 100.f;
    float Vf;
    if (V > rgbTinv_32f) {
        Vf = rgbInvScale_32f * pow(V, rgbLPowInv_32f) - rgbLShift_32f;
    } else {
        Vf = rgbSmallScale_32f * V; 
    }
    //Vf = (Vf < 0) ? 0 : Vf;
    return (unsigned char) roundf(Vf * 255.f);
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

/**
 * @brief Peicewise function to help calculate conversion XYZ to Lab
 *   
 * @param t 
 * @return float output of function
 */
float f_func(float t){
    if (t > lab_delta3_32f) {
        return cbrt(t);
    } else {
        return t*labSmallScale_32f + labSmallShift_32f;
    }
}

float XYZtoL(float X, float Y, float Z) {
    Y = f_func(Y);

    return Y*labLScale_32f - labLShift_32f;

}

float XYZtoA(float X, float Y, float Z) {
    X = f_func(X);
    Y = f_func(Y);
    
    return 500.f*(X - Y);
}

float XYZtoB(float X, float Y, float Z) {
    Y = f_func(Y);
    Z = f_func(Z);

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
/**
 * @brief Peicewise function to help calculate conversion Lab to XYZ
 *   
 * @param t 
 * @return float output of function
 */
float fInv_func(float t){
    if ( t > lab_delta_32f) {
        return pow(t, labPow_32f);
    } else {
        return labSmallScaleInv_32f*(t - labSmallShift_32f);
    }
}

float LABtoX(float L, float a, float b) {
    float fy = (L + labLShift_32f)/labLScale_32f;
    float fx = fy + a*labAScale_32f;
    
    return fInv_func(fx);
}

float LABtoY(float L, float a, float b) {
    float fy = (L + labLShift_32f)/labLScale_32f;

    return fInv_func(fy);
}

float LABtoZ(float L, float a, float b) {
    float fy = (L + labLShift_32f)/labLScale_32f;
    float fz = fy - b*labBScale_32f;

    return fInv_func(fz);
}

void rgb2lab(int R, int G, int B, float *L, float *a, float *b) {
    // printf("_RGB2LAB_\n");
    // printf("(R:%d,G:%d,B:%d)->",R,G,B);
    float Rf = linearize(R);
    float Gf = linearize(G);
    float Bf = linearize(B);
    
    // printf("(Rf:%f,Bf:%f,Gf:%f)->",Rf,Gf,Bf);

    float X = RGBtoX(Rf, Gf, Bf);
    float Y = RGBtoY(Rf, Gf, Bf);
    float Z = RGBtoZ(Rf, Gf, Bf);
    
    // printf("(X:%f,Y:%f,Z:%f)->",X,Y,Z);

    X *= labXScale_32f;
    Y *= labYScale_32f;
    Z *= labZScale_32f;
    
    // printf("(XS:%f,YS:%f,ZS:%f)->",X,Y,Z);

    *L = XYZtoL(X, Y, Z);
    *a = XYZtoA(X, Y, Z);
    *b = XYZtoB(X, Y, Z);
    
    // printf("(L:%f,a:%f,b:%f)\n",*L,*a,*b);
}

void lab2rgb(float L, float a, float b, unsigned char *R, unsigned char *G, unsigned char *B) {
    // printf("_LAB2RGB_\n");
    // printf("(L:%f,a:%f,b:%f)->",L,a,b);
    float X = LABtoX(L, a, b);
    float Y = LABtoY(L, a, b);
    float Z = LABtoZ(L, a, b);
    
    // printf("(XS:%f,YS:%f,ZS:%f)->",X,Y,Z);
    X *= labXScaleInv_32f;
    Y *= labYScaleInv_32f;
    Z *= labZScaleInv_32f;

    // printf("(X:%f,Y:%f,Z:%f)->",X,Y,Z);

    float Rf = XYZtoR(X, Y, Z);
    float Gf = XYZtoG(X, Y, Z);
    float Bf = XYZtorgB(X, Y, Z);

    // printf("(Rf:%f,Bf:%f,Gf:%f)->",Rf,Gf,Bf);

    *R = delin(Rf);
    *G = delin(Gf);
    *B = delin(Bf);
    // printf("(R:%d,G:%d,B:%d)",*R,*G,*B);
}
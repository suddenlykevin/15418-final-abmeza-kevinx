/**
 * @file pixImage.cpp
 * @author Kevin Xie () 
 *          Anthony Meza (abmeza)
 * @brief Implementation of PixImage
 * 
 * @version 0.1
 * @date 2022-04-22
 * 
 * @copyright Copyright (c) 2022
 * 
 */


// Import util libraries
#include "util/colorConv.h"
#include "util/superPixel.h"
#include "pixImage.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <algorithm> 
#include <stack>
#include <cmath>
using namespace std;

// computes eigendecomposition of real 3x3 hermitian matrix and 
// returns maximum eigenvalue/vector
//
// closed-form solution based on https://hal.archives-ouvertes.fr/hal-01501221/document
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


/**
 * @brief Performs the SILC algorithm to help figure out the superpixel contents
 * 
 * @param m 
 * @param S 
 * @param l_k 
 * @param a_k 
 * @param b_k 
 * @param x_k 
 * @param y_k 
 * @param l_i 
 * @param a_i 
 * @param b_i 
 * @param x_i 
 * @param y_i 
 * @return float 
 */
float dist_k(int m, float S, float l_k, float a_k, float b_k, int x_k, int y_k,
           float l_i, float a_i, float b_i, int x_i, int y_i) {
    float d_lab = sqrt(pow(l_k - l_i, 2.f) + pow(a_k - a_i, 2.f) + pow(b_k - b_i, 2.f));
    float d_xy = sqrt(pow((float) x_k - x_i, 2.f) + pow((float) y_k - y_i, 2.f));
    float k = ((float) m) / S;
    return d_lab + k * d_xy;
}

//implements a gaussian function with std "sigma" and mean "mean"
float gaussian(float x, float sigma, float mean) {
    return exp((x-mean)*(x-mean)/(-2.0f*sigma*sigma))/sqrt(6.28319*sigma*sigma);
}

/**
 * @brief wrapper for malloc to check for errors in allocation
 * 
 * @param size   size of item being allocated
 * @return void* pointer to space allocated
 */
void* wrp_malloc(size_t size){ 
    void* ptr = malloc(size);

    // Check that no error occured
    if (ptr == NULL) {
        printf("Unable to allocate memory with malloc of size %zu\n", size);
        exit(1);
    }

    return ptr;
}

/**
 * @brief wrapper for calloc to check for errors in allocation
 * 
 * @param nitems number of items to allocated
 * @param size   size of item being allocated
 * @return void* pointer to space allocated
 */
void* wrp_calloc(size_t nitems, size_t size){ 
    void* ptr = calloc(nitems, size);

    // Check that no error occured
    if (ptr == NULL) {
        printf("Unable to allocate memory with calloc of items %zu and size %zu\n", nitems, size);
        exit(1);
    }

    return ptr;
}




PixImage :: PixImage(unsigned char* input_image, int in_w, int in_h, int out_w, int out_h, int K){
    // Store our new variables
    input_img = input_image;
    in_width = in_w;
    in_height = in_h;
    out_width = out_w;
    out_height = out_h;
    K_colors = K;

    // Get value for number of pixels
    M_pix = in_width * in_height;
    N_pix = out_width * out_height;
}


/**
 * @brief Initializes the superPixel_pos array and the region_map array
 */
void PixImage :: initSuperPixels(){
    // Get change in length of values
    float dx = (float) in_width/(float) out_width;
    float dy = (float) in_height/(float) out_height;

    // initialize superpixel positions (centers)
    for (int j = 0; j < out_height; ++j) {
        for (int i = 0; i < out_width; ++i) {

            // Calculate midpoint value
            float x = ((float) i + 0.5f) * dx;
            float y = ((float) j + 0.5f) * dy;
            FloatVec pos =  (FloatVec) {x,y};
            // Set value
            superPixel_pos[out_width * j + i] = pos;
            
            // ********** DEBUG ********** //
            //printf("superpixel %d: (%f, %f)\n", out_width * j + i, x, y);
            // ********** DEBUG ********** //
        }
    }

    // Initial assignment of pixels to a specific superpxel  
    for (int j = 0; j < in_height; ++j) {
        for (int i = 0; i < in_width; ++i) {
            // Calculate which superpixel to set
            int x = (int) ((float) i / dx);
            int y = (int) ((float) j / dy);

            // Set Value
            region_map[in_width * j + i] = out_width * y + x;

            // ********** DEBUG ********** //
            //printf("Pixel %d: (%d)\n", in_width * j + i, out_width * y + x);
            // ********** DEBUG ********** //
        }
    }
}

void PixImage :: updateSuperPixelMeans(){
    FloatVec sp_sums[N_pix];
    memset(sp_sums, 0, N_pix*sizeof(FloatVec));
    LabColor color_sums[N_pix];
    memset(color_sums, 0, N_pix*sizeof(LabColor));
    int sp_count[N_pix];
    memset(sp_count, 0, N_pix*sizeof(int));
    // Find the mean colors (from input image) for each superpixel
    for (int j = 0; j < in_height; j++) {
        for (int i = 0; i < in_width; i++) {
            int idx = j*in_width + i;
            int spidx = region_map[idx];
            sp_count[spidx] ++;
            sp_sums[spidx].x += i;
            sp_sums[spidx].y += j;

            color_sums[spidx].L += input_img_lab[idx].L;
            color_sums[spidx].a += input_img_lab[idx].a;
            color_sums[spidx].b += input_img_lab[idx].b;
        }
    }
    
    // Repostion superpixels and update the output color pallete
    for (int j = 0; j < out_height; j++) {
        for (int i = 0; i < out_width; i++) {
            // Index of superpixel
            int spidx = j*out_width + i;

            // Calculate new position for super pixel
            float x = sp_sums[spidx].x / sp_count[spidx];
            float y = sp_sums[spidx].y / sp_count[spidx];
            FloatVec newpos = {x, y};
            superPixel_pos[spidx] = newpos;

            // Set output_img_lab to new mean value
            sp_mean_lab[spidx].L = color_sums[spidx].L/sp_count[spidx];
            sp_mean_lab[spidx].a = color_sums[spidx].a/sp_count[spidx];
            sp_mean_lab[spidx].b = color_sums[spidx].b/sp_count[spidx];
        }
    }
}

void PixImage :: pushPaletteColor(LabColor color, float prob) {
    palette_lab[palette_size] = color;
    prob_c[palette_size] = prob;
    for (int idx = 0; idx < N_pix; idx ++) {
        prob_c_if_sp[palette_size*N_pix + idx] = prob;
    }
    palette_size ++; 
}

void PixImage :: pushPalettePair(int a, int b) {
    PalettePair newPair = {a, b};
    int idx = (palette_size >> 1) - 1;
    if (idx < 0 || idx > K_colors + 1) {
        return;
    }
    palette_pairs[idx] = newPair;
}

void PixImage :: getAveragedPalette(LabColor *avg_palette) {
    for (int i = 0; i < palette_size >> 1; i++) {
        PalettePair pair = palette_pairs[i];
        float weight_a = prob_c[pair.a];
        float weight_b = prob_c[pair.b];
        float total_weight = weight_a + weight_b;
        weight_a /= total_weight;
        weight_b /= total_weight;

        LabColor ca = palette_lab[pair.a];
        LabColor cb = palette_lab[pair.b];

        LabColor avg = {ca.L*weight_a + cb.L*weight_b,
                        ca.a*weight_a + cb.a*weight_b,
                        ca.b*weight_a + cb.b*weight_b};
        
        avg_palette[pair.a] = avg;
        avg_palette[pair.b] = avg;
    }
}

void PixImage :: initialize(){
    ///*** Allocate Array Space ***///
    input_img_lab = (LabColor *) wrp_calloc(M_pix, sizeof(LabColor));

    superPixel_pos = (FloatVec *) wrp_calloc(N_pix, sizeof(FloatVec)); 
    region_map = (int *) wrp_calloc(M_pix, sizeof(int));

    palette_lab = (LabColor *) wrp_calloc((K_colors + 1) * 2 , sizeof(float));
    palette_size = 0;
    palette_pairs = (PalettePair *) wrp_calloc(K_colors + 1, sizeof(PalettePair));
    palette_assign = (int *) wrp_calloc(N_pix, sizeof(int));
    prob_c = (float *) wrp_calloc((K_colors + 1) * 2 , sizeof(float)); 
    prob_c_if_sp = (float *) wrp_calloc((K_colors + 1) * 2 * N_pix, sizeof(float));
    prob_sp = 1.0f/(out_width*out_height); 

    buf_lab = (LabColor *) wrp_calloc(N_pix, sizeof(LabColor));
    output_img = (unsigned char *) wrp_malloc(N_pix * 3); 
    sp_mean_lab = (LabColor *) wrp_calloc(N_pix, sizeof(LabColor)); 

    ///*** Create input_img_lab version ***///
    unsigned char *p; LabColor *pl;
    for(p = input_img, pl = input_img_lab; p != input_img + (M_pix*3); p += 3, pl ++) 
        rgb2lab(*p, *(p+1), *(p+2), &(pl->L), &(pl->a), &(pl->b));

    ///*** Initialize Superpixel Values ***///
    initSuperPixels();
    updateSuperPixelMeans();

    ///*** Initialize Palette Values ***///
    // Find mean of all M_pix color inputs
    LabColor color_sum;
    //add all colors
    for (int p = 0; p < N_pix; p++) {
        color_sum.L += sp_mean_lab[p].L;
        color_sum.a += sp_mean_lab[p].a;
        color_sum.b += sp_mean_lab[p].b;
    }
    // divide all by M_pix
    color_sum.L = color_sum.L * prob_sp;
    color_sum.a = color_sum.a * prob_sp;
    color_sum.b = color_sum.b * prob_sp;

    // Store color and update prob to any
    pushPaletteColor(color_sum, 0.5f);
    LabColor majorAxis;
    float variance;
    getMajorAxis(0, &variance, &majorAxis);
    
    color_sum.L += majorAxis.L * kSubclusterPertubation;
    color_sum.a += majorAxis.a * kSubclusterPertubation;
    color_sum.b += majorAxis.b * kSubclusterPertubation;

    pushPaletteColor(color_sum, 0.5f);

    pushPalettePair(0, 1);

    T = sqrt(2*variance) * kT0SafetyFactor;
}

void PixImage :: iterate(){
    
    // Size of super pixel on input image
    float S = sqrt(((float) (M_pix))/((float) (N_pix)));

    // update superpixel segments
    for (int iter = 0; iter < 35; ++iter) {

        //*** ************************ ***//
        //*** (4.2) REFINE SUPERPIXELS ***//
        //*** ************************ ***//

        //*** ************************ ***//
        //*** (4.2.1) ASSOCIATE SUPERPIXELS ***//
        //*** ************************ ***//

        ///*** Update boundaries of pixels assosiated with super pixels ***///
        LabColor average_palette[(K_colors + 1)*2];
        getAveragedPalette(average_palette);

        float distance[M_pix];
        for (int i = 0; i < M_pix; i++) distance[i] = -1.0f;

        for (int j = 0; j < out_height; ++j) {
            for (int i = 0; i < out_width; ++i) {
                
                // get local region
                int idx = out_width * j + i;
                FloatVec center = superPixel_pos[idx];
                int min_x = std::max(0.0f, center.x - S);
                int min_y = std::max(0.0f, center.y - S);
                int max_x = std::min((float) (in_width - 1), center.x + S);
                int max_y = std::min((float) (in_height - 1), center.y + S);            
                //printf("iter %d superpixel %d: (%d, %d) -> (%d, %d)\n", iter, out_width * j + i, min_x, min_y, max_x, max_y);
                int x = (int) round(center.x);
                int y = (int) round(center.y);

                LabColor sp_color = palette_lab[palette_assign[idx]];

                // within region
                for (int yy = min_y; yy <= max_y; ++yy) {
                    for (int xx = min_x; xx <= max_x; ++xx) {
                        int curr_idx = yy * in_width + xx;

                        // check new distance
                        float dist_new = dist_k(m_gerstner, S, sp_color.L, sp_color.a, sp_color.b, 
                                                x, y, input_img_lab[curr_idx].L, input_img_lab[curr_idx].a, 
                                                input_img_lab[curr_idx].b, xx, yy);

                        // Check if the distance is less in order to minimize
                        if (distance[curr_idx] < 0 || dist_new < distance[curr_idx]) {
                            distance[curr_idx] = dist_new;
                            region_map[curr_idx] = out_width*j + i;
                        }
                    }
                }
            }
        }

        updateSuperPixelMeans();

        // smooth positions
        for (int j = 0; j < out_height; j++) {
            for (int i = 0; i < out_width; i++) {                
                int spidx = j*out_width + i;
                FloatVec sum = {0, 0};
                float count = 0.0f;
                if(i > 0) {
                    sum.x += superPixel_pos[j*out_width + i-1].x;
                    sum.y += superPixel_pos[j*out_width + i-1].y;
                    count += 1.0f;
                }
                if(i < out_width -1) {
                    sum.x += superPixel_pos[j*out_width + i+1].x;
                    sum.y += superPixel_pos[j*out_width + i+1].y;
                    count += 1.0f;
                }
                if(j > 0) {
                    sum.x += superPixel_pos[(j-1)*out_width + i].x;
                    sum.y += superPixel_pos[(j-1)*out_width + i].y;
                    count += 1.0f;
                }
                if(j < out_height - 1) {
                    sum.x += superPixel_pos[(j+1)*out_width + i].x;
                    sum.y += superPixel_pos[(j+1)*out_width + i].y;
                    count += 1.0f;
                }
                sum.x /= count;
                sum.y /= count;
                FloatVec pos = superPixel_pos[spidx];
                FloatVec newPos = {0, 0};
                if(i == 0 || i == out_width -1) {
                    newPos.x = pos.x;
                } else {
                    newPos.x = (0.6f)*pos.x + 0.4f*sum.x;
                }
                if(j == 0 || j == out_height - 1) {
                    newPos.y = pos.y;
                } else {
                    newPos.y = 0.6f*pos.y + 0.4f*sum.y;
                }
                // printf("pos: (%f, %f) -> (%f, %f)\n", pos.x, pos.y, newPos.x, newPos.y);
                superPixel_pos[spidx] = newPos;
            }
        }
    
        // smooth colors
        for(int j = 0; j < out_height; ++j) {
            for(int i = 0; i < out_width; ++i) {

            //get bounds of 3x3 kernel (make sure we don't go off the image)
            int min_x = std::max(0,i-1);
            int max_x = std::min(out_width-1,i+1);
            int min_y = std::max(0,j-1);
            int max_y = std::min(out_height-1,j+1);

            //Initialize
            LabColor sum = {0.f, 0.f, 0.f};
            float weight = 0.f;

            //get current SP color and (grid) position
            LabColor superpixel_color = sp_mean_lab[j*out_width + i];
            FloatVec p = {(float) j, (float) i};

            //get bilaterally weighted average color of SP neighborhood
            for(int ii = min_x; ii<= max_x; ++ii) {
                for(int jj = min_y; jj<=max_y; ++jj) {
                
                LabColor c_n = sp_mean_lab[jj*out_width + ii];
                float d_color = (float) (pow(superpixel_color.L - c_n.L, 2.f) +
                                             pow(superpixel_color.a - c_n.a, 2.f) +
                                             pow(superpixel_color.b - c_n.b, 2.f));
                float w_color = gaussian(d_color, 2.0f ,0.0f);
                float d_pos = (float) (pow(i-ii, 2.f) + pow(j-jj, 2.f));
                float w_pos = gaussian(d_pos, 0.97f, 0.0f);
                float w_total = w_color*w_pos;

                weight += w_total;
                sum.L += c_n.L*w_total;
                sum.a += c_n.a*w_total;
                sum.b += c_n.b*w_total;
                }
            }
            sum.L *= 1.0f/weight;
            sum.a *= 1.0f/weight;
            sum.b *= 1.0f/weight;
            buf_lab[j*out_width + i] = sum;
            }
        }

        
        printf("smoothed colors\n");
        //update the SP mean colors with the smoothed values
        memcpy(sp_mean_lab, buf_lab, N_pix * sizeof(LabColor));
        
        printf("memcpy\n");
        
        //*** ************************************** ***//
        //*** (4.3) ASSOSIATE SUPERPIXELS TO PALETTE ***//
        //*** ************************************** ***//
        // List of P(c_k|p_s) values for all superpixels
        // prob_c_if_sp = (float *) wrp_calloc(N_pix*K_colors, sizeof(float));   
        // // Update superpixel colors from color palette based on P(c_k|p_s) calculation
        // for(int p = 0; p < N_pix; p++) {
        //     // Get the best color value to update the superpixel color
        //     int best_c = -1;
        //     float best_norm_val = 0.0f;
        //     for (int c = 0; c < palette_size; c++){
        //         // m_s' - c_k TODO: MIGHT NOT WORK?
        //         LabColor pixDiff;
        //         pixDiff.L = output_img_lab[p].L - palette_lab[c].L;
        //         pixDiff.a = output_img_lab[p].a - palette_lab[c].a;
        //         pixDiff.b = output_img_lab[p].b - palette_lab[c].b;

        //         // || m_s' - c_k ||
        //         float norm_val = sqrt(pow(pixDiff.L, 2.f) + pow(pixDiff.a, 2.f) + pow(pixDiff.b, 2.f));

        //         //  - (|| m_s' - c_k ||/T)
        //         float pow_val = -1.0f*(norm_val/T);
        //         prob_c_if_sp[c + (p*K_colors)] = prob_c[c] * exp(pow_val);
        //         //Update if better value
        //         if (best_c == -1 || norm_val < best_norm_val){
        //             best_c = c;
        //             best_norm_val = norm_val;
        //         }
        //     } 
        //     // Store values to update output_img_lab to best pixel value
        //     buf_lab[p] = palette_lab[best_c];
        // }

        // // Update to probability color is assosiated to ANY superpixel (P(c_k))
        // // TODO: OPTIMIZATION EXISTS WHERE PREVIOUS LOOPS IS MIXED WITH THIS LOOP
        // for (int c = 0; c < palette_size; c++){
        //     prob_c[c] = 0.0f;
        //     for (int p = 0; p < N_pix; p++){
        //         prob_c[c] = prob_c[c] + prob_c_if_sp[c + (p*K_colors)] * prob_sp[p];
        //     }
        // }


        // //update the SP colors wiht updated values
        // memcpy(output_img_lab, buf_lab, N_pix * sizeof(LabColor));
        
        // //*** ******************** ***//
        // //*** (4.3) REFINE PALETTE ***//
        // //*** ******************** ***//

        // //TODO: DIFF FROM THERE IMPLEMENTATION? CHECK?
        // for (int c = 0; c< palette_size; c++){

        //     LabColor c_sum = {0.0f,0.0f,0.0f};
        //     // Observe all superpixels to get sum of equation
        //     for (int p = 0; p < N_pix; p++){
        //         c_sum.L = c_sum.L + output_img_lab[p].L * prob_c_if_sp[c + (p*K_colors)] * prob_sp[p];
        //         c_sum.a = c_sum.a + output_img_lab[p].a * prob_c_if_sp[c + (p*K_colors)] * prob_sp[p];
        //         c_sum.b = c_sum.b + output_img_lab[p].b * prob_c_if_sp[c + (p*K_colors)] * prob_sp[p];
        //     }

        //     //Update palette color
        //     palette_lab[c].L = c_sum.L/prob_c[c];
        //     palette_lab[c].a = c_sum.a/prob_c[c];
        //     palette_lab[c].b = c_sum.b/prob_c[c];

        // }

        // free(prob_c_if_sp);
        //*** ******************** ***//
        //*** (4.3) EXPAND PALETTE ***//
        //*** ******************** ***//
    
}
    

    //*** ******************** ***//
    //*** PROCESS OUTPUT IMAGE ***//
    //*** ******************** ***//

    // Create output image in rgb color values
    for (int j = 0; j < out_height; j++) {
        for (int i = 0; i < out_width; i++) {
            int idx = j*out_width + i;
            lab2rgb(sp_mean_lab[idx].L, sp_mean_lab[idx].a, sp_mean_lab[idx].b, 
                    &(output_img[3*idx]), &(output_img[3*idx + 1]), &(output_img[3*idx + 2]));
        }
    }

}

void PixImage :: getMajorAxis(int palette_index, float *value, LabColor *vector) {
    float covariance[9];
    memset(covariance, 0, 9*sizeof(float));
    float sum = 0;

    // compute covariance matrix
    for (int j = 0; j < out_height; j++) {
        for (int i = 0; i < out_height; i++) {
            int idx = j*out_width + i;

            // probability of superpixel given palette color
            float prob_oc = prob_c_if_sp[palette_index * N_pix + idx] 
                            * prob_sp / prob_c[palette_index];
            sum += prob_oc;
            
            // find color error with current superpixel
            LabColor pl_color = palette_lab[palette_index];    
            LabColor sp_color = sp_mean_lab[idx];
            float L_error = fabs(pl_color.L - sp_color.L);
            float a_error = fabs(pl_color.a - sp_color.a);
            float b_error = fabs(pl_color.b - sp_color.b);

            // update covariance
            covariance[0] += prob_oc*L_error*L_error;
            covariance[1] += prob_oc*a_error*L_error;
            covariance[2] += prob_oc*b_error*L_error;
            covariance[3] += prob_oc*L_error*a_error;
            covariance[4] += prob_oc*a_error*a_error;
            covariance[5] += prob_oc*b_error*a_error;
            covariance[6] += prob_oc*L_error*b_error;
            covariance[7] += prob_oc*a_error*b_error;
            covariance[8] += prob_oc*b_error*b_error;
        }
    }

    LabColor eVec;
    float eVal;

    int error = maxEigen3(covariance, &eVal, &eVec);
    if (error < 0) {
        printf("maxEigen3 special case\n");
    }

    float len = sqrt(eVec.L*eVec.L + eVec.a*eVec.a + eVec.b*eVec.b);
    if (len > 0) {
        eVec.L *= (1.0f/len);
        eVec.a *= (1.0f/len);
        eVec.b *= (1.0f/len);
    }

    *value = eVal;
    *vector = eVec;
}

void PixImage :: freeAll(){

    free(input_img_lab);

    free(output_img); 
    free(sp_mean_lab); 

    free(superPixel_pos);
    free(region_map);

    free(palette_lab); 
    free(prob_c); 

    free(buf_lab);
}

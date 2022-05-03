/**
 * @file pixImage.cu
 * @author Kevin Xie (kevinx) 
 *         Anthony Meza (abmeza)
 * @brief Implementation of PixImage. Also a replica off pixImage.cpp, but 
 *        created in order to implement cuda. Some of the structure was 
 *        inspired by assignment 2 in 15418 which also used cuda to render
 *        circles on a grid.
 * 
 * @version 0.1
 * @date 2022-04-22
 * 
 * @copyright Copyright (c) 2022
 * 
 */

// Constants to regulate what we print
//#define DEBUG     // misc. debug statements
#define RUN_DEBUG // debug statements that check running progress
#define TIMING // Calculate and print timing information


// Import util libraries
#include "CycleTimer.h"
#include "pixImage.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <algorithm> 
#include <stack>
#include <cmath>
using namespace std;


//********************************************************//
//*******************  GLOBAL VARIABLES ******************//
//********************************************************//

// storage for global constants that usble by GPU
struct GlobalConstants {
    // Input Image Content 
    int in_width, in_height;  //<- width and height of input_img, gives pixel dimensions
    int M_pix;                //<- # of pixels from the input image (M from paper)
    unsigned char *input_img; //<- input image loaded, uses rgb values for pixels (0-255)
    LabColor *input_img_lab;  //<- input image, using cielab values 

    // Output Image 
    int out_width, out_height; //<- output version of width, height
    int N_pix;                 //<- # of pixels in the output image (N from paper)
    unsigned char *output_img; //<- output version of input_img
    unsigned char *spoutput_img; //<- debug output for superpixels
    LabColor *buf_lab;         //<- buffer for smoothing and palette refinement

    // Superpixel calculation 
    FloatVec *superPixel_pos; //<- Super pixel coordinate positions "on input image"
    LabColor *sp_mean_lab;  //<- superpixel mean color value
    int *region_map;      //<- array with values for pixels associated with a specific superpixel 

    // Palette 
    int K_colors;          //<- number of colors we aim to use in the pallette
    int * palette_size;      //<- POINTER SO WE CAN MODIFY Current # of colors stored in palette_lab
    PalettePair *palette_pairs;
    int *palette_assign; //<- palette assignment for each superpixel
    LabColor *palette_lab; //<- palette array with color values in palette
    LabColor *average_palette; //<- average palette array with average color values
   
    bool *palette_complete; //<- POINTER SO WE CAN MODIFY

    float *prob_c;         //<- array of probabiities that a color in the palette is set to ANY super pixel
    float prob_sp;         //<- array of probabiities of each super pixel TODO:EDIT
    float *prob_c_if_sp;   //<- List of P(c_k|p_s) values for all superpixels
    
    // Temperature
    float *T;   //<- POINTER SO WE CAN MODIFY Current temperature

    //Bool thing
    bool *converged
};

// Constant for GPU
__constant__ GlobalConstants cuGlobalConsts;

//********************************************************//
//*******************  KERNAL FUNCTIONS ******************//
//********************************************************//

/**
 * @brief creates lab version of input_img
 */
__global__ void kernelCreateInputLAB() {

    int index = blockIdx.x * blockDim.x + threadIdx.x;

    //TODO: RUN ON ONE KERNAL FOR NOW
    if (index == 0){
    
    // *** TODO TRANSFER OVER CONSTANTS ***//
    int M_pix = cuGlobalConsts.M_pix;
    unsigned char *input_img = cuGlobalConsts.input_img;
    LabColor *input_img_lab = cuGlobalConsts.input_img_lab;


    unsigned char *p; 
    LabColor *pl;
    for(p = input_img, pl = input_img_lab; p != input_img + (M_pix*3); p += 3, pl ++) 
        rgb2lab(*p, *(p+1), *(p+2), &(pl->L), &(pl->a), &(pl->b));

    
    }
}

/**
 * @brief kernal that runs initSuperPixels on Device
 */
__global__ void kernelInitSuperPixels() {

    int index = blockIdx.x * blockDim.x + threadIdx.x;

    //TODO: RUN ON ONE KERNAL FOR NOW
    if (index == 0){
    
    // *** TODO TRANSFER OVER CONSTANTS ***//
    int in_width = cuGlobalConsts.in_width;
    int in_height = cuGlobalConsts.in_height;
    int out_width = cuGlobalConsts.out_width;
    int out_height = cuGlobalConsts.out_height;
    FloatVec *superPixel_pos = cuGlobalConsts.superPixel_pos;
    FloatVec *region_map = cuGlobalConsts.region_map;


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


        }
    }

    }
}
/**
 * @brief kernal that runs updateSuperPixelMeans on Device
 */
__global__ void kernelUpdateSuperPixelMeans() {
    
    int index = blockIdx.x * blockDim.x + threadIdx.x;

    //TODO: RUN ON ONE KERNAL FOR NOW
    if (index == 0){
    
    // *** TODO TRANSFER OVER CONSTANTS ***//
    int N_pix = cuGlobalConsts.N_pix;
    int in_width = cuGlobalConsts.in_width;
    int in_height = cuGlobalConsts.in_height;
    int out_width = cuGlobalConsts.out_width;
    int out_height = cuGlobalConsts.out_height;

    int *region_map = cuGlobalConsts.region_map;
    LabColor *input_img_lab = cuGlobalConsts.input_img_lab;
    FloatVec *superPixel_pos = cuGlobalConsts.superPixel_pos;


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
}

/**
 * @brief kernal that runs initializes palette values on Device
 */
__global__ void kernelInitPaletteValues() {
    
    int index = blockIdx.x * blockDim.x + threadIdx.x;

    //TODO: RUN ON ONE KERNAL FOR NOW
    if (index == 0){
    

    // *** TODO TRANSFER OVER CONSTANTS ***//
    int N_pix = cuGlobalConsts.N_pix;

    float *T = cuGlobalConsts.T;
    LabColor *sp_mean_lab = cuGlobalConsts.sp_mean_lab;
    float *prob_c_if_sp = cuGlobalConsts.prob_c_if_sp;


    // Find mean of all M_pix color inputs
    LabColor color_sum = {0.f, 0.f, 0.f};

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

    #ifdef DEBUG
    printf("color_init: (%f, %f, %f)\n", color_sum.L, color_sum.a, color_sum.b);
    #endif

    // Store color and update prob to any
    pushPaletteColor(color_sum, 0.5f);
    for (int idx = 0; idx < N_pix; idx ++) {
        prob_c_if_sp[idx] = 0.5f;
    }
    LabColor majorAxis;
    float variance;
    getMajorAxis(0, &variance, &majorAxis);
    
    color_sum.L += majorAxis.L * kSubclusterPertubation;
    color_sum.a += majorAxis.a * kSubclusterPertubation;
    color_sum.b += majorAxis.b * kSubclusterPertubation;

    pushPaletteColor(color_sum, 0.5f);
    for (int idx = 0; idx < N_pix; idx ++) {
        prob_c_if_sp[N_pix + idx] = 0.5f;
    }

    pushPalettePair(0, 1);

    *T = sqrt(2*variance) * kT0SafetyFactor;

    }
}

/**
 * @brief Update the average palette
 */
__global__ void kernelGetAveragedPalette() {
    
    int index = blockIdx.x * blockDim.x + threadIdx.x;

    //TODO: RUN ON ONE KERNAL FOR NOW
    if (index == 0){
    

    // *** TODO TRANSFER OVER CONSTANTS ***//

    LabColor *average_palette = cuGlobalConsts.average_palette;
    LabColor *palette_lab = cuGlobalConsts.palette_lab;
    int *palette_size = cuGlobalConsts.palette_size;

    PalettePair *palette_pairs = cuGlobalConsts.palette_pairs;
    float *prob_c = cuGlobalConsts.prob_c;

    bool *palette_complete = cuGlobalConsts.palette_complete;


    if (*palette_complete) {
        memcpy(average_palette, palette_lab, (*palette_size) *sizeof(LabColor));
        return;
    }
    for (int i = 0; i < (*palette_size) >> 1; i++) {
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
        
        average_palette[pair.a] = avg;
        average_palette[pair.b] = avg;
    }
    
    }
}

/**
 * @brief Assosiate Pixels to specific super pixels
 */
__global__ void kernelAssosiatetoSuperPixels() {
    
    int index = blockIdx.x * blockDim.x + threadIdx.x;

    //TODO: RUN ON ONE KERNAL FOR NOW
    if (index == 0){
    

    // *** TODO TRANSFER OVER CONSTANTS ***//
    int M_pix = cuGlobalConsts.M_pix;
    int out_height = cuGlobalConsts.out_height;
    int out_width = cuGlobalConsts.out_width;

    int *region_map = cuGlobalConsts.region_map;
    FloatVec *superPixel_pos = cuGlobalConsts.superPixel_pos;
    LabColor *palette_lab = cuGlobalConsts.palette_lab;
    int *palette_assign = cuGlobalConsts.palette_assign;
    LabColor *input_img_lab = cuGlobalConsts.input_img_lab;

    //Global bois
    float *distance = (float *)wrp_calloc(M_pix, sizeof(float));


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

    free(distance);
    }
}

/**
 * @brief Assosiate Pixels to specific super pixels
 */
__global__ void kernelSmoothPositions() {
    
    int index = blockIdx.x * blockDim.x + threadIdx.x;

    //TODO: RUN ON ONE KERNAL FOR NOW
    if (index == 0){
    

    // *** TODO TRANSFER OVER CONSTANTS ***//
    int M_pix = cuGlobalConsts.M_pix;
    int out_height = cuGlobalConsts.out_height;
    int out_width = cuGlobalConsts.out_width;

    FloatVec *superPixel_pos = cuGlobalConsts.superPixel_pos;
    LabColor *palette_lab = cuGlobalConsts.palette_lab;
    int *palette_assign = cuGlobalConsts.palette_assign;
    LabColor *input_img_lab = cuGlobalConsts.input_img_lab;


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
                newPos.x = (0.55f)*pos.x + 0.45f*sum.x;
            }
            if(j == 0 || j == out_height - 1) {
                newPos.y = pos.y;
            } else {
                newPos.y = 0.55f*pos.y + 0.45f*sum.y;
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
            float d_color = (float) sqrt(pow(superpixel_color.L - c_n.L, 2.f) +
                                            pow(superpixel_color.a - c_n.a, 2.f) +
                                            pow(superpixel_color.b - c_n.b, 2.f));
            float w_color = gaussian(d_color, 2.0f ,0.0f);
            float d_pos = (float) sqrt(pow(i-ii, 2.f) + pow(j-jj, 2.f));
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
    
    //update the SP mean colors with the smoothed values
    memcpy(sp_mean_lab, buf_lab, N_pix * sizeof(LabColor));

    }
}

/**
 * @brief Assosiate Pixels to specific super pixels
 */
__global__ void kernelAssosiateToPalette() {
    
    int index = blockIdx.x * blockDim.x + threadIdx.x;

    //TODO: RUN ON ONE KERNAL FOR NOW
    if (index == 0){
    

    // *** TODO TRANSFER OVER CONSTANTS ***//
    int N_pix = cuGlobalConsts.N_pix;
    int K_colors = cuGlobalConsts.K_colors;
    int prob_sp = cuGlobalConsts.prob_sp;
    

    int *palette_size = cuGlobalConsts.palette_size;
    float *prob_c_if_sp = cuGlobalConsts.prob_c_if_sp;
    LabColor *sp_mean_lab = cuGlobalConsts.sp_mean_lab;
    LabColor *palette_lab = cuGlobalConsts.palette_lab;
    int *palette_assign = cuGlobalConsts.palette_assign;
    float *prob_c = cuGlobalConsts.prob_c;


    float new_prob_c[(*palette_size)];
    memset(new_prob_c, 0, (*palette_size)*sizeof(float));
    memset(prob_c_if_sp, 0, K_colors * 2 * N_pix *sizeof(float));
    
    // Update superpixel colors from color palette based on P(c_k|p_s) calculation
    for(int p = 0; p < N_pix; p++) {
        // Get the best color value to update the superpixel color
        int best_c = -1;
        float best_norm_val = 0.0f;
        double probs[(*palette_size)];
        double sum_prob = 0;

        for (int c = 0; c < (*palette_size); c++){

            // m_s' - c_k TODO: MIGHT NOT WORK?
            LabColor pixDiff;
            pixDiff.L = sp_mean_lab[p].L - palette_lab[c].L;
            pixDiff.a = sp_mean_lab[p].a - palette_lab[c].a;
            pixDiff.b = sp_mean_lab[p].b - palette_lab[c].b;

            // || m_s' - c_k ||
            float norm_val = sqrt(pow(pixDiff.L, 2.f) + pow(pixDiff.a, 2.f) + pow(pixDiff.b, 2.f));

            //  - (|| m_s' - c_k ||/T)
            float pow_val = -1.0f*(norm_val/T);
            float prob = prob_c[c] * exp(pow_val);
            
            probs[c] = prob;
            sum_prob += prob;

            //Update if better value
            if (best_c < 0 || norm_val < best_norm_val){
                best_c = c;
                best_norm_val = norm_val;
            }
        } 

        // update palette assignment
        palette_assign[p] = best_c;

        for(int c = 0; c < (*palette_size); c++) {
            double p_norm = probs[c]/sum_prob;
            prob_c_if_sp[c*(N_pix) + p] = p_norm;
            new_prob_c[c] += prob_sp*p_norm;
        }
    }
    
    
    // update color probabilities
    memcpy(prob_c, new_prob_c, (*palette_size)*sizeof(float));
        
    }


}

/**
 * @brief Assosiate Pixels to specific super pixels
 */
__global__ void kernelAssosiateToPalette() {
    
    int index = blockIdx.x * blockDim.x + threadIdx.x;

    //TODO: RUN ON ONE KERNAL FOR NOW
    if (index == 0){
    

    // *** TODO TRANSFER OVER CONSTANTS ***//
    int N_pix = cuGlobalConsts.N_pix;
    int K_colors = cuGlobalConsts.K_colors;
    int prob_sp = cuGlobalConsts.prob_sp;
    

    int *palette_size = cuGlobalConsts.palette_size;
    float *prob_c_if_sp = cuGlobalConsts.prob_c_if_sp;
    LabColor *sp_mean_lab = cuGlobalConsts.sp_mean_lab;
    LabColor *palette_lab = cuGlobalConsts.palette_lab;
    int *palette_assign = cuGlobalConsts.palette_assign;
    float *prob_c = cuGlobalConsts.prob_c;


    float new_prob_c[(*palette_size)];
    memset(new_prob_c, 0, (*palette_size)*sizeof(float));
    memset(prob_c_if_sp, 0, K_colors * 2 * N_pix *sizeof(float));
    
    // Update superpixel colors from color palette based on P(c_k|p_s) calculation
    for(int p = 0; p < N_pix; p++) {
        // Get the best color value to update the superpixel color
        int best_c = -1;
        float best_norm_val = 0.0f;
        double probs[(*palette_size)];
        double sum_prob = 0;

        for (int c = 0; c < (*palette_size); c++){

            // m_s' - c_k TODO: MIGHT NOT WORK?
            LabColor pixDiff;
            pixDiff.L = sp_mean_lab[p].L - palette_lab[c].L;
            pixDiff.a = sp_mean_lab[p].a - palette_lab[c].a;
            pixDiff.b = sp_mean_lab[p].b - palette_lab[c].b;

            // || m_s' - c_k ||
            float norm_val = sqrt(pow(pixDiff.L, 2.f) + pow(pixDiff.a, 2.f) + pow(pixDiff.b, 2.f));

            //  - (|| m_s' - c_k ||/T)
            float pow_val = -1.0f*(norm_val/T);
            float prob = prob_c[c] * exp(pow_val);
            
            probs[c] = prob;
            sum_prob += prob;

            //Update if better value
            if (best_c < 0 || norm_val < best_norm_val){
                best_c = c;
                best_norm_val = norm_val;
            }
        } 

        // update palette assignment
        palette_assign[p] = best_c;

        for(int c = 0; c < (*palette_size); c++) {
            double p_norm = probs[c]/sum_prob;
            prob_c_if_sp[c*(N_pix) + p] = p_norm;
            new_prob_c[c] += prob_sp*p_norm;
        }
    }
    
    
    // update color probabilities
    memcpy(prob_c, new_prob_c, (*palette_size)*sizeof(float));
        
    }


}


/**
 * @brief Refine the palette
 */
__global__ void kernelRefinePalette() {
    
    int index = blockIdx.x * blockDim.x + threadIdx.x;

    //TODO: RUN ON ONE KERNAL FOR NOW
    if (index == 0){
    

    // *** TODO TRANSFER OVER CONSTANTS ***//
    int N_pix = cuGlobalConsts.N_pix;
    int K_colors = cuGlobalConsts.K_colors;
    int prob_sp = cuGlobalConsts.prob_sp;
    
    int *palette_size = cuGlobalConsts.palette_size;
    float *prob_c_if_sp = cuGlobalConsts.prob_c_if_sp;
    LabColor *sp_mean_lab = cuGlobalConsts.sp_mean_lab;
    float *prob_c = cuGlobalConsts.prob_c;
    LabColor *palette_lab = cuGlobalConsts.palette_lab;
    PalettePair *palette_pairs = cuGlobalConsts.palette_pairs;
    bool *palette_complete = cuGlobalConsts.palette_complete;

    bool *converged = cuGlobalConsts.converged;


    float palette_error = 0.f;

    //TODO: DIFF FROM THERE IMPLEMENTATION? CHECK?
    for (int c = 0; c < (*palette_size); c++){

        LabColor c_sum = {0.0f,0.0f,0.0f};
        // Observe all superpixels to get sum of equation
        for (int p = 0; p < N_pix; p++){
            c_sum.L += sp_mean_lab[p].L * prob_c_if_sp[c*N_pix + p] * prob_sp;
            c_sum.a += sp_mean_lab[p].a * prob_c_if_sp[c*N_pix + p] * prob_sp;
            c_sum.b += sp_mean_lab[p].b * prob_c_if_sp[c*N_pix + p] * prob_sp;
        }
        
        if (prob_c[c] > 0) {
            LabColor last = palette_lab[c];
            //Update palette color
            palette_lab[c].L = c_sum.L/prob_c[c];
            palette_lab[c].a = c_sum.a/prob_c[c];
            palette_lab[c].b = c_sum.b/prob_c[c];
            LabColor curr = palette_lab[c];

            palette_error += sqrt(pow(last.L-curr.L, 2.0f) + pow(last.a-curr.a, 2.0f) + pow(last.b-curr.b, 2.0f));
        }
    }
   
    #ifdef RUN_DEBUG
    printf("expand... ");
    #endif

    if (palette_error < kPaletteErrorTolerance) {
        // check for convergence, lower temperature
        if (T <= kTF) {
            (*converged) = true;
        } else {
            T = std::max(T*kDT, kTF);
        }
        
        // if palette is incomplete
        if (!(*palette_complete)) {
            int splits[K_colors];
            int curr = 0;
            for (int i = 0; i < (*palette_size) >> 1; i++) {
                #ifdef DEBUG
                // printf("(%d, %d)\n", palette_pairs[i].a, palette_pairs[i].b);
                #endif

                LabColor color_a = palette_lab[palette_pairs[i].a];
                LabColor color_b = palette_lab[palette_pairs[i].b];

                float error = sqrt(pow(color_a.L-color_b.L, 2.0f) + 
                                pow(color_a.a-color_b.a, 2.0f) + 
                                pow(color_a.b-color_b.b, 2.0f));
                // printf("%f, (%f, %f, %f), (%f, %f, %f)\n", error, color_a.L,color_a.a,color_a.b, color_b.L, color_b.a, color_b.b);
                // determine if split or simply perturb 
                if (error > kSubclusterTolerance) {
                    splits[curr] = i;
                    curr ++;
                } else {
                    float value;
                    LabColor majorAxis;
                    getMajorAxis(palette_pairs[i].a, &value, &majorAxis);
                    color_b.L += majorAxis.L * kSubclusterPertubation;
                    color_b.a += majorAxis.a * kSubclusterPertubation;
                    color_b.b += majorAxis.b * kSubclusterPertubation;
                    // printf("perturbed by: (%f, %f, %f)\n", majorAxis.L, majorAxis.a, majorAxis.b);

                    palette_lab[palette_pairs[i].b] = color_b;
                }
            }

            // should sort splits by distance here.
            if (curr > 0) {
                #ifdef RUN_DEBUG
                printf("expanding... %d, %d", (*palette_size), curr);
                #endif
            }

            for (int i = 0; i < curr; i++) {
                splitColor(splits[i]);

                // if full, seal palette
                if ((*palette_size) >= 2 * K_colors) {            
                    #ifdef RUNF_DEBUG
                    printf("COMPLETE\n");                        
                    #endif

                    condensePalette();
                    break;
                }
            }
        }
    }



    }


}


/**
 * @brief Process the output image
 */
__global__ void kernelProcessOutputImage() {
    
    int index = blockIdx.x * blockDim.x + threadIdx.x;

    //TODO: RUN ON ONE KERNAL FOR NOW
    if (index == 0){
    

    // *** TODO TRANSFER OVER CONSTANTS ***//
    int out_width = cuGlobalConsts.out_width;
    int out_height = cuGlobalConsts.out_height;
    int in_width = cuGlobalConsts.in_width;
    int in_height = cuGlobalConsts.in_height;

    
    LabColor *average_palette = cuGlobalConsts.average_palette;
    int *palette_assign = cuGlobalConsts.palette_assign;
    unsigned char  *output_img = cuGlobalConsts.output_img;
    unsigned char  *spoutput_img = cuGlobalConsts.spoutput_img;
    int *region_map = cuGlobalConsts.region_map;

    // Create the output image
    for (int j = 0; j < out_height; j++) {
        for (int i = 0; i < out_width; i++) {
            int idx = j*out_width + i;
            LabColor color = averaged_palette[palette_assign[idx]];

            lab2rgb(color.L, color.a, color.b, 
                    &(output_img[3*idx]), &(output_img[3*idx + 1]), &(output_img[3*idx + 2]));

        }
    }

    // Create superpixel output image
    for (int j = 0; j < in_height; j++) {
        for (int i = 0; i < in_width; i++) {
            int idx = j*in_width + i;
            if ((region_map[idx]/out_width % 2 == 0 && region_map[idx] % 2 == 0) ||
                (region_map[idx]/out_width % 2 == 1 && region_map[idx] % 2 == 1)) {
                spoutput_img[idx*3] = 0;
                spoutput_img[idx*3 + 1] = 0;
                spoutput_img[idx*3 + 2] = 0;
            } else {
                spoutput_img[idx*3] = 255;
                spoutput_img[idx*3 + 1] = 255;
                spoutput_img[idx*3 + 2] = 255;
            }
        }
    }



    }


}


//********************************************************//
//*******************  PIXEL FUNCTIONS  ******************//
//********************************************************//


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

    // Initialze basic values
    palette_complete = false;
    palette_size = 0; 
    
    #ifdef TIMING
        //Timing variables
        startAllTime=0.f;
        endAllTime=0.f; 
    #endif

    // Init Them Arrays
    input_img = NULL; 
    input_img_lab = NULL;

    output_img = NULL; 
    spoutput_img = NULL; 
    buf_lab = NULL;  

    superPixel_pos = NULL; 
    sp_mean_lab = NULL; 
    region_map = NULL;  

    palette_pairs = NULL;
    palette_assign = NULL; 
    palette_lab = NULL;  
    average_palette = NULL; 

    prob_c = NULL;      
    prob_sp = 0.0f;      
    prob_c_if_sp = NULL;   

    T = 0.0f;  


    // Cuda Device versions of values               
    cuDev_input_img = NULL; 
    cuDev_input_img_lab = NULL;  

    cuDev_output_img = NULL; 
    cuDev_spoutput_img = NULL; 
    cuDev_buf_lab = NULL;        

    cuDev_superPixel_pos = NULL; 
    cuDev_sp_mean_lab = NULL;  
    cuDev_region_map = NULL;     
    
    cuDev_palette_size = NULL;  //single pointer value
    cuDev_palette_pairs = NULL;
    cuDev_palette_assign = NULL;
    cuDev_palette_lab = NULL;
    cuDev_average_palette = NULL; 
    cuDev_palette_complete = NULL;  //single pointer value

    cuDev_prob_c = NULL;         
    cuDev_prob_c_if_sp = NULL; 

    cuDev_T = NULL;  
}


/**
 * @brief Initializes the superPixel_pos array and the region_map array
 */
void PixImage :: initSuperPixels(){
    
    // Intialize size of kernal
    dim3 blockDim(in_width, in_height, 1);
    dim3 gridDim(1,1);


    kernelInitSuperPixels<<<gridDim, blockDim>>>();
    cudaDeviceSynchronize();

}

void PixImage :: updateSuperPixelMeans(){
    
    // Intialize size of kernal
    dim3 blockDim(in_height,in_width, 1);
    dim3 gridDim(1,1);


    kernelUpdateSuperPixelMeans<<<gridDim, blockDim>>>();
    cudaDeviceSynchronize();
}

void PixImage :: pushPaletteColor(LabColor color, float prob) {
    palette_lab[palette_size] = color;
    prob_c[palette_size] = prob;
    palette_size ++; 
}

void PixImage :: pushPalettePair(int a, int b) {
    PalettePair newPair = {a, b};
    int idx = (palette_size >> 1) - 1;
    if (idx < 0 || idx > K_colors) {
        return;
    }
    palette_pairs[idx] = newPair;
}

void PixImage :: getAveragedPalette() {

    // Intialize size of kernal
    dim3 blockDim(1, 1, 1);
    dim3 gridDim(1,1);


    kernelGetAveragedPalette<<<gridDim, blockDim>>>();
    cudaDeviceSynchronize();

}

void PixImage :: splitColor(int pair_index) {
    int index_a = palette_pairs[pair_index].a;
    int index_b = palette_pairs[pair_index].b;

    int next_a = palette_size;
    int next_b = next_a + 1;

    // perturb a
    LabColor color_a = palette_lab[index_a];
    LabColor color_a_b = color_a;
    LabColor majorAxis;
    float variance;
    getMajorAxis(index_a, &variance, &majorAxis);
    color_a_b.L += majorAxis.L * kSubclusterPertubation;
    color_a_b.a += majorAxis.a * kSubclusterPertubation;
    color_a_b.b += majorAxis.b * kSubclusterPertubation;
    
    // reconstruct pair a and copy probabilities
    prob_c[index_a] *= 0.5f;
    pushPaletteColor(color_a_b, prob_c[index_a]);
    int index_a_b = palette_size - 1;
    memcpy(&prob_c_if_sp[index_a_b*N_pix], &prob_c_if_sp[index_a*N_pix], N_pix*sizeof(float));
    palette_pairs[pair_index].b = index_a_b;

    // perturb b
    LabColor color_b = palette_lab[index_b];
    LabColor color_b_b = color_b;
    getMajorAxis(index_b, &variance, &majorAxis);
    color_b_b.L += majorAxis.L * kSubclusterPertubation;
    color_b_b.a += majorAxis.a * kSubclusterPertubation;
    color_b_b.b += majorAxis.b * kSubclusterPertubation;

    // reconstruct pair b and copy probabilities
    prob_c[index_b] *= 0.5f;
    pushPaletteColor(color_b_b, prob_c[index_b]);
    int index_b_b = palette_size - 1;
    memcpy(&prob_c_if_sp[index_b_b*N_pix], &prob_c_if_sp[index_b*N_pix], N_pix*sizeof(float));
    pushPalettePair(index_b, index_b_b);
}

void PixImage :: condensePalette() {
    LabColor new_palette[K_colors * 2];
    float new_prob_c[K_colors * 2];
    float new_prob_c_if_sp[K_colors * 2 * N_pix];
    int new_palette_assign[N_pix];
    
    #ifdef RUN_DEBUG
    printf("averaging... %d ", palette_size);
    #endif

    getAveragedPalette();

    #ifdef RUN_DEBUG
    printf("averaged... ");
    #endif

    // for each pair, condense to average
    for(int j = 0; j < palette_size >> 1; ++j) {
        int index_a = palette_pairs[j].a;
        int index_b = palette_pairs[j].b;
        new_palette[j] = average_palette[index_a]; //TODO: SHOULD BE IN KERNALS
        //update the probability of the single superpixel
        new_prob_c[j] = prob_c[index_a] + 
                        prob_c[index_b];

        // reassign superpixels
        for(int i = 0; i < N_pix; i++) {
            new_prob_c_if_sp[j*N_pix + i] = prob_c_if_sp[index_a*N_pix + i];
            if (palette_assign[i] == index_a || palette_assign[i] == index_b) {
                new_palette_assign[i] = j;
            }
        }
    }

    // copy new values
    memcpy(palette_lab, new_palette, K_colors * 2 * sizeof(LabColor));
    memcpy(palette_assign, new_palette_assign, N_pix * sizeof(int));
    
    // TODO: could be wrong?? wtf is prob_oc_
    memcpy(prob_c_if_sp, new_prob_c_if_sp, K_colors * 2 * N_pix * sizeof(float));

    palette_size = K_colors;
    palette_complete = true;

    
    #ifdef RUN_DEBUG
    printf("DONE \n");
    #endif
}

void PixImage :: initVariables(){
    ///*** Allocate Array Space ***///
    input_img_lab = (LabColor *) wrp_calloc(M_pix, sizeof(LabColor));

    superPixel_pos = (FloatVec *) wrp_calloc(N_pix, sizeof(FloatVec)); 
    region_map = (int *) wrp_calloc(M_pix, sizeof(int));

    palette_lab = (LabColor *) wrp_calloc(K_colors * 2 , sizeof(LabColor));
    average_palette = (LabColor *) wrp_calloc(K_colors * 2 , sizeof(LabColor)); 
    palette_size = 0;
    palette_pairs = (PalettePair *) wrp_calloc(K_colors, sizeof(PalettePair));
    palette_assign = (int *) wrp_calloc(N_pix, sizeof(int));
    prob_c = (float *) wrp_calloc(K_colors * 2 , sizeof(float)); 
    prob_c_if_sp = (float *) wrp_calloc(K_colors * 2 * N_pix, sizeof(float));
    prob_sp = 1.0f/(out_width*out_height); 

    buf_lab = (LabColor *) wrp_calloc(N_pix, sizeof(LabColor));
    output_img = (unsigned char *) wrp_malloc(N_pix * 3); 
    spoutput_img = (unsigned char *) wrp_calloc(M_pix*3, sizeof(unsigned char));
    sp_mean_lab = (LabColor *) wrp_calloc(N_pix, sizeof(LabColor)); 

    // Allocate space for Device, place in global for easy access

    cudaMalloc(&cuDev_input_img, 3*M_pix*sizeof(unsigned char));
    cudaMalloc(&cuDev_input_img_lab, M_pix* sizeof(LabColor));
    cudaMalloc(&cuDev_output_img, N_pix * 3);
    cudaMalloc(&cuDev_spoutput_img, M_pix*3*sizeof(unsigned char));
    cudaMalloc(&cuDev_buf_lab, N_pix*sizeof(LabColor));
    cudaMalloc(&cuDev_superPixel_pos, N_pix * sizeof(FloatVec));
    cudaMalloc(&cuDev_sp_mean_lab, N_pix * sizeof(LabColor));
    cudaMalloc(&cuDev_region_map, M_pix * sizeof(int));
    cudaMalloc(&cuDev_palette_pairs, K_colors*sizeof(PalettePair));
    cudaMalloc(&cuDev_palette_assign, N_pix *sizeof(int));
    cudaMalloc(&cuDev_palette_lab, K_colors * 2 *sizeof(LabColor));
    cudaMalloc(&cuDev_average_palette, K_colors * 2 *sizeof(LabColor));
    cudaMalloc(&cuDev_prob_c_if_sp, K_colors * 2 * N_pix * sizeof(float));
        
    cudaMalloc(&cuDev_palette_size , sizeof(int));
    cudaMalloc(&cuDev_palette_complete, sizeof(bool));
    cudaMalloc(&cuDev_T , sizeof(float));
    cudaMalloc(&cuDev_converged , sizeof(bool));
    
    // Set values to cuda variables (if calloc simply memset)

    cudaMemcpy(cuDev_input_img, input_img, 3*M_pix*sizeof(unsigned char), cudaMemcpyHostToDevice);
    cudaMemset(cuDev_input_img_lab, 0, M_pix* sizeof(LabColor));
    cudaMemset(cuDev_output_img, 0, N_pix * 3);
    cudaMemset(cuDev_spoutput_img, 0, M_pix*3*sizeof(unsigned char));
    cudaMemset(cuDev_buf_lab, 0, N_pix*sizeof(LabColor));
    cudaMemset(cuDev_superPixel_pos, 0, N_pix * sizeof(FloatVec));
    cudaMemset(cuDev_sp_mean_lab, 0, N_pix * sizeof(LabColor));
    cudaMemset(cuDev_region_map, 0, M_pix * sizeof(int));
    cudaMemset(cuDev_palette_pairs, 0, K_colors*sizeof(PalettePair));
    cudaMemset(cuDev_palette_assign,0, N_pix *sizeof(int));
    cudaMemset(cuDev_palette_lab, 0, K_colors * 2 *sizeof(LabColor));
    cudaMemset(cuDev_average_palette, 0, K_colors * 2 *sizeof(LabColor));
    cudaMemset(cuDev_prob_c_if_sp,0, K_colors * 2 * N_pix * sizeof(float));
        
    cudaMemset(cuDev_palette_size , 0, sizeof(int));
    cudaMemset(cuDev_palette_complete,false, sizeof(bool));
    cudaMemset(cuDev_T , 0.0f, sizeof(float));
    cudaMemset(cuDev_converged, false, sizeof(bool));
    

    // Initialize parameters in constant memory in order to take 
    // advantge of cuda optimizations. 

    GlobalConstants params;

    // NON pointer values
    params.in_width = in_width;
    params.in_height = in_height;
    params.M_pix = M_pix; 

    params.out_width = out_width;
    params.out_height = out_height;
    params.N_pix = N_pix;

    params.K_colors = K_colors;
    params.palette_size = palette_size;

    params.prob_sp = prob_sp;

    // Pointer Values
    params.input_img = cuDev_input_img;
    params.input_img_lab = cuDev_input_img_lab;

    params.output_img = cuDev_output_img;
    params.spoutput_img = cuDev_spoutput_img;
    params.buf_lab = cuDev_buf_lab;

    params.superPixel_pos = cuDev_superPixel_pos;
    params.sp_mean_lab = cuDev_sp_mean_lab;
    params.region_map = cuDev_region_map;

    params.palette_pairs = cuDev_palette_pairs;
    params.palette_assign = cuDev_palette_assign;
    params.palette_lab = cuDev_palette_lab;
    params.average_palette = cuDev_average_palette;
    params.palette_complete = cuDev_palette_complete;
    
    params.prob_c = cuDev_prob_c;
    params.prob_c_if_sp = cuDev_prob_c_if_sp;

    params.T = cuDev_T; 
    params.converged = cuDev_converged;   
    
    cudaMemcpyToSymbol(cuGlobalConsts, &params, sizeof(GlobalConstants)); 

}
void PixImage :: initialize(){
    initVariables();

    ///*** Create input_img_lab version ***///
    dim3 blockDim(1, 1, 1);
    dim3 gridDim(1,1);

    kernelCreateInputLAB<<<gridDim, blockDim>>>();
    cudaDeviceSynchronize();

    ///*** Initialize Superpixel Values ***///
    initSuperPixels();
    updateSuperPixelMeans();


    ///*** Initialize Palette Values ***///
    dim3 blockDim(in_width, in_height, 1);
    dim3 gridDim(1,1);

    kernelInitPaletteValues<<<gridDim, blockDim>>>();
    cudaDeviceSynchronize();

}

void PixImage :: runPixelate(){
    
    #ifdef TIMING
    startAllTime = CycleTimer::currentSeconds();
    #endif

    //*** ******************** ***//
    //*** (4.1) INITIALIZATION ***//
    //*** ******************** ***//
    #ifdef TIMING
    startInitializeTime = CycleTimer::currentSeconds();
    #endif
    
    initialize();
    
    #ifdef TIMING
    endInitializeTime = CycleTimer::currentSeconds();
    #endif
    //*** ******************* ***//
    //*** CORE ALGORITHM LOOP ***//
    //*** ******************* ***//

    // Size of super pixel on input image
    float S = sqrt(((float) (M_pix))/((float) (N_pix)));


    
    bool *converged = wrp_malloc(sizeof(bool));
    *converged = false;
    
    int iter = 0;

    // update superpixel segments
    while (!(*converged) && iter < maxIter) {
        
        #ifdef RUN_DEBUG
        printf("iter %d, %f\n", iter, T);
        #endif

        //*** ************************ ***//
        //*** (4.2) REFINE SUPERPIXELS ***//
        //*** ************************ ***//

        ///*** Update boundaries of pixels assosiated with super pixels ***///
        
        #ifdef RUN_DEBUG
        printf("average...");
        #endif

        ///*** Get average colors for palette ***///
        getAveragedPalette();
        
        #ifdef RUN_DEBUG
        printf("DONE\n");
        #endif'

        #ifdef RUN_DEBUG
        printf("associate...\n");
        #endif

        ///*** Assosiate to superpixels ***///

        dim3 blockDim(1, 1, 1);
        dim3 gridDim(1,1);

        kernelAssosiatetoSuperPixels<<<gridDim, blockDim>>>();
        cudaDeviceSynchronize();

        
        #ifdef RUN_DEBUG 
        printf("udpate means...");
        #endif
        
        updateSuperPixelMeans();

        #ifdef RUN_DEBUG
        printf("DONE\n");
        #endif

        #ifdef RUN_DEBUG
        printf("smooth...");
        #endif
        
        ///*** Smooth positions of Superpixel and pixel ***///

        dim3 blockDim(1, 1, 1);
        dim3 gridDim(1,1);

        kernelSmoothPositions<<<gridDim, blockDim>>>();
        cudaDeviceSynchronize();

        #ifdef RUN_DEBUG
        printf("DONE\n");
        #endif


        //*** ************************************** ***//
        //*** (4.3) ASSOSIATE SUPERPIXELS TO PALETTE ***//
        //*** ************************************** ***//
        #ifdef RUN_DEBUG
        printf("associate...");
        #endif

        dim3 blockDim(1, 1, 1);
        dim3 gridDim(1,1);

        kernelAssosiateToPalette<<<gridDim, blockDim>>>();
        cudaDeviceSynchronize();


        #ifdef RUN_DEBUG
        printf("DONE\n");
        #endif

        // //*** ***************************** ***//
        // //*** (4.3) REFINE + EXPAND PALETTE ***//
        // //*** ***************************** ***//

        
        #ifdef RUN_DEBUG
        printf("refine...");
        #endif
        
        dim3 blockDim(1, 1, 1);
        dim3 gridDim(1,1);

        kernelRefinePalette<<<gridDim, blockDim>>>();
        cudaDeviceSynchronize();

        #ifdef RUN_DEBUG
        printf("DONE\n");
        #endif

        //Transfer converged from device (TODO: maybe bug)
        cudaMemcpy(converged, cuDev_converged, sizeof(bool), cudaMemcpyDeviceToHost);
    
  
        iter ++;
    }

    free(converged);

    //*** ******************** ***//
    //*** PROCESS OUTPUT IMAGE ***//
    //*** ******************** ***//
    // palette_complete = false;
    // Create output image in rgb color values
    getAveragedPalette();


    dim3 blockDim(1, 1, 1);
    dim3 gridDim(1,1);

    kernelProcessOutputImage<<<gridDim, blockDim>>>();
    cudaDeviceSynchronize();

    //Transfer new image stuff from device (TODO: maybe bug)
    cudaMemcpy(output_img, cuDev_output_img, N_pix * 3, cudaMemcpyDeviceToHost);
    cudaMemcpy(spoutput_img, cuDev_spoutput_img, M_pix*3, cudaMemcpyDeviceToHost);
    

    #ifdef TIMING
    endAllTime = CycleTimer::currentSeconds();
    #endif


    #ifdef TIMING
    //*** ******************** ***//
    //*** PRINT TIMING RESULTS ***//
    //*** ******************** ***//
    // Print Total Time to run algorithm (not include get image file and create image file)
    printf("Overall: %.3f s\n", 1000.f * (endAllTime - startAllTime));
    
    printf("\t- Initialize: %.3f s\n", 1000.f * (endInitializeTime - startInitializeTime));
    #endif

}

void PixImage :: getMajorAxis(int palette_index, float *value, LabColor *vector) {
    float covariance[9];
    memset(covariance, 0, 9*sizeof(float));
    float sum = 0;
    #ifdef DEBUG
    //printf("fsds %d, %f, %f\n", palette_index, prob_sp, prob_c[palette_index]);
    #endif

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

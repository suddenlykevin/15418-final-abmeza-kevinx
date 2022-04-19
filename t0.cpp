/**
 * @file t0.cpp
 * @authors Kevin Xie () 
 *          Anthony Meza (abmeza)
 * 
 * @brief  This has been taken from an online tutorial to figure out how to use
 *         stb_image files for our project. Link found: 
 *         https://solarianprogrammer.com/2019/06/10/c-programming-reading-writing-images-stb_image-libraries/
 * @version 0.1
 * @date 2022-04-09
 * 
 * @copyright Copyright (c) 2022
 * 
 * NOTE: Running the files will require the following prompt:
 *       gcc -Wall -pedantic t0.cpp -o t0 -lm
 *       
 * 
 */


// Import util libraries
#include "util/colorConv.h"
#include "util/superPixel.h"

// Import stb_image libraries
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image/stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image/stb_image_write.h"


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <algorithm> 
#include <stack>


/**
 * @brief position vector data structure
 */
typedef struct {
    float x;   //< x coor
    float y;   //< y coor
} FloatVec;

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




int main(void) {
    
    //*** ********* ***//
    //*** VARIABLES ***//
    //*** ********* ***//
    
    // Values of input image
    int width, height, channels;
    // Values of output image
    int out_width = 16;
    int out_height = 16;
    // Some thingy TODO
    int m_gerstner = 45;
 
    //*** ********** ***//
    //*** LOAD IMAGE ***//
    //*** ********** ***//

    unsigned char *img = stbi_load("SonicFlower.jpeg", &width, &height, &channels, 0);
    if(img == NULL) {
        printf("Error in loading the image\n");
        exit(1);
    }

    //********** DEGUB **********//
    printf("Loaded image with a width of %dpx, a height of %dpx and %d channels\n", width, height, channels);
    //********** DEGUB **********//
    
    
    //*** ******************** ***//
    //*** (4.1) INITIALIZATION ***//
    //*** ******************** ***//

    // Establish superpixel data structures 
    size_t img_size = width * height * channels;
    size_t out_img_size = out_width * out_height * channels;
    unsigned char *out_img = (unsigned char *) malloc(out_img_size);

    // Generate array of super pixels position 
    //    Specifically creating an array with the center value on a superpixel on the original image
    FloatVec *superpixel_pos = (FloatVec *) calloc(out_width * out_height, sizeof(FloatVec)); 
    if (superpixel_pos == NULL) {
        printf("Unable to allocate memory for the superpixel_pos image.\n");
        exit(1);
    }
    
    // initialize superpixel positions (centers)
    for (int j = 0; j < out_height; ++j) {
        for (int i = 0; i < out_width; ++i) {
            float x = (i + 0.5f) * width / out_width;
            float y = (j + 0.5f) * height / out_height;
            FloatVec pos = {x, y};
            superpixel_pos[out_width * j + i] = pos;
            
            // ********** DEBUG ********** //
            //printf("superpixel %d: (%f, %f)\n", out_width * j + i, x, y);
            // ********** DEBUG ********** //
        }
    }

    // Initial assignment of pixels to a specific superpxel
    int *superpixel_img = (int *) calloc(width * height, sizeof(int));  
    for (int j = 0; j < height; ++j) {
        for (int i = 0; i < width; ++i) {
            float dx = (float) width/(float) out_width;
            float dy = (float) height/(float) out_height;
            int x = (int) ((float) i / dx);
            int y = (int) ((float) j / dy);
            superpixel_img[width * j + i] = out_width * y + x;
            
            // ********** DEBUG ********** //
            //printf("Pixel %d: (%d)\n", width * j + i, out_width * y + x);
            // ********** DEBUG ********** //
        }
    }

    // Create and convert images into lab
    // Create image pixels in the lab space for input and output
    float *lab_image = (float *) calloc(img_size, sizeof(float));   // input lab image
    float *lab_out = (float *) calloc(out_img_size, sizeof(float)); // output lab image 
    if (lab_image == NULL || lab_out == NULL) {
        printf("Unable to allocate memory for the cielab image.\n");
        exit(1);
    }
    // Loop through pixels to save them in lab space
    unsigned char *p;
    float *pl;
    for(p = img, pl = lab_image; p != img + img_size; p += channels, pl += channels) {
        rgb2lab(*p, *(p + 1), *(p + 2), pl, pl + 1, pl + 2);
    }

    
    //*** ******************* ***//
    //*** CORE ALGORITHM LOOP ***//
    //*** ******************* ***//

    // Size of super pixel on input image
    float S = sqrt(((float) (width * height))/((float) (out_width * out_height)));
    
    // update superpixel segments
    for (int iter = 0; iter < 35; ++iter) {

        // update boundaries
        for (int j = 0; j < out_height; ++j) {
            for (int i = 0; i < out_width; ++i) {
                
                // get local region
                FloatVec center = superpixel_pos[out_width * j + i];
                int min_x = std::max(0.0f, center.x - S);
                int min_y = std::max(0.0f, center.y - S);
                int max_x = std::min((float) (width - 1), center.x + S);
                int max_y = std::min((float) (height - 1), center.y + S);            
                //printf("iter %d superpixel %d: (%d, %d) -> (%d, %d)\n", iter, out_width * j + i, min_x, min_y, max_x, max_y);
                int x = (int) round(center.x);
                int y = (int) round(center.y);
                int idx = y*width + x;
                
                // within region
                for (int yy = min_y; yy <= max_y; ++yy) {
                    for (int xx = min_x; xx <= max_x; ++xx) {
                        int curr_idx = yy * width + xx;
                        int curr_spidx = superpixel_img[curr_idx];
                        FloatVec curr_spixel = superpixel_pos[curr_spidx];
                        int curr_spx = (int) round(curr_spixel.x);
                        int curr_spy = (int) round(curr_spixel.y);
                        curr_spidx = curr_spy * width + curr_spx;
                        float dist_curr = dist_k(m_gerstner, S, lab_image[curr_spidx*channels], 
                                                lab_image[curr_spidx*channels + 1], lab_image[curr_spidx*channels + 2], 
                                                curr_spx, curr_spy, lab_image[curr_idx*channels], lab_image[curr_idx*channels + 1], 
                                                lab_image[curr_idx*channels + 2], xx, yy);
                        float dist_new = dist_k(m_gerstner, S, lab_image[idx*channels], 
                                                lab_image[idx*channels + 1], lab_image[idx*channels + 2], 
                                                x, y, lab_image[curr_idx*channels], lab_image[curr_idx*channels + 1], 
                                                lab_image[curr_idx*channels + 2], xx, yy);
                        if (dist_new < dist_curr) {
                            superpixel_img[curr_idx] = out_width*j + i;
                        }
                    }
                }
            }
        }


        // Variables to find mean color values
        FloatVec *sp_sums = (FloatVec *) calloc(out_width * out_height, sizeof(FloatVec));
        float *color_sums = (float *) calloc(out_img_size, sizeof(float));
        int *sp_count = (int *) calloc(out_width * out_height, sizeof(int)); 

        // Find the mean colors (from input image) for each superpixel
        for (int j = 0; j < height; j++) {
            for (int i = 0; i < width; i++) {
                int idx = j*width + i;
                int spidx = superpixel_img[idx];
                sp_count[spidx] ++;
                sp_sums[spidx].x += i;
                sp_sums[spidx].y += j;

                color_sums[3*spidx] += lab_image[3*idx];
                color_sums[3*spidx + 1] += lab_image[3*idx + 1];
                color_sums[3*spidx + 2] += lab_image[3*idx + 2];
            }
        }
        
        // Repostion superpixels and update the output color pallete
        for (int j = 0; j < out_height; j++) {
            for (int i = 0; i < out_width; i++) {
                // Index of superpixel
                int spidx = j*out_width + i;

                // Calcualte new position for super pixel
                float x = sp_sums[spidx].x / sp_count[spidx];
                float y = sp_sums[spidx].y / sp_count[spidx];
                FloatVec newpos = {x, y};
                superpixel_pos[spidx] = newpos;

                // Set lab_out to new mean value
                lab_out[3*spidx] = color_sums[3*spidx]/sp_count[spidx];
                lab_out[3*spidx + 1] = color_sums[3*spidx + 1]/sp_count[spidx];
                lab_out[3*spidx + 2] = color_sums[3*spidx + 2]/sp_count[spidx];
                
            }
        }

        // smooth
        for (int j = 0; j < out_height; j++) {
            for (int i = 0; i < out_width; i++) {                
                int spidx = j*out_width + i;
                FloatVec sum = {0, 0};
                float count = 0.0f;
                if(i > 0) {
                    sum.x += superpixel_pos[j*out_width + i-1].x;
                    sum.y += superpixel_pos[j*out_width + i-1].y;
                    count += 1.0f;
                }
                if(i < out_width -1) {
                    sum.x += superpixel_pos[j*out_width + i+1].x;
                    sum.y += superpixel_pos[j*out_width + i+1].y;
                    count += 1.0f;
                }
                if(j > 0) {
                    sum.x += superpixel_pos[(j-1)*out_width + i].x;
                    sum.y += superpixel_pos[(j-1)*out_width + i].y;
                    count += 1.0f;
                }
                if(j < out_height - 1) {
                    sum.x += superpixel_pos[(j+1)*out_width + i].x;
                    sum.y += superpixel_pos[(j+1)*out_width + i].y;
                    count += 1.0f;
                }
                sum.x /= count;
                sum.y /= count;
                FloatVec pos = superpixel_pos[spidx];
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
                superpixel_pos[spidx] = newPos;
            }
        }
    }

    // for (int j = 0; j < height; j++) {
    //     for (int i = 0; i < width; i++) {
    //         int idx = j*width + i;
    //         int spidx = superpixel_img[idx];
    //         lab2rgb(lab_out[3*spidx], lab_out[3*spidx + 1], lab_out[3*spidx + 2], &out_img[3*idx], &out_img[3*idx + 1], &out_img[3*idx + 2]);
    //         // if ((superpixel_img[idx]/out_width % 2 == 0 && superpixel_img[idx] % 2 == 0) ||
    //         //     (superpixel_img[idx]/out_width % 2 == 1 && superpixel_img[idx] % 2 == 1)) {
    //         //     out_img[idx*3] = 0;
    //         //     out_img[idx*3 + 1] = 0;
    //         //     out_img[idx*3 + 2] = 0;
    //         // } else {
    //         //     out_img[idx*3] = 255;
    //         //     out_img[idx*3 + 1] = 255;
    //         //     out_img[idx*3 + 2] = 255;
    //         // }
    //     }
    // }

    // Create output image in rgb color values
    for (int j = 0; j < out_height; j++) {
        for (int i = 0; i < out_width; i++) {
            int idx = j*out_width + i;
            lab2rgb(lab_out[3*idx], lab_out[3*idx + 1], lab_out[3*idx + 2], &out_img[3*idx], &out_img[3*idx + 1], &out_img[3*idx + 2]);
        }
    }

    // Passed in the correct number of channels, in this case the number desired
    stbi_write_jpg("SonicFlower_output.jpeg", out_width, out_height, channels, out_img, 100);
    

    /// TESTING TO SEE WHAT IS DIFFERENT WITH THE TWO OUTPUT IMAGES ///
    
    // image_diff("SonicFlower.jpeg" , "SonicFlower_gray.jpeg", 1);
    
    /// TESTING TO SEE WHAT IS DIFFERENT WITH THE TWO OUTPUT IMAGES ///
 
    stbi_image_free(img);
}




// DEBUG HELPER FUNCTION STUFF


/**
 * @brief Check how in range a float is to consider it "equal"
 * 
 * @param p1 first float value
 * @param p2 second pixel value
 * @return int 
 */
int in_range(unsigned char* p1, unsigned char* p2){
    float value = ( ((float)*p1) - ((float) *p2));
    
    // The range we will consider it to be "equal"
    if (abs(value) > 2){
        return 1;
    }
    return 0;
}


/**
 * @brief find what pixels are different between two images
 * 
 * @param image_1_name file name of first image
 * @param image_2_name file name of second image
 * @param print_content state whether or not we want to print differnt pixels
 */
int image_diff(char* image_1_name, char* image_2_name, int print_content){
    printf("DEBUG: Running image_diff\n");
    // Initilize variables
    int widthg, heightg, channelsg;
    int widthb, heightb, channelsb;

    // Get content of first image
    unsigned char *img_g = stbi_load(image_1_name, &widthg, &heightg, &channelsg, 0);
    if(img_g == NULL) {
        printf("Error in loading the image\n");
        exit(1);
    }
    printf("Loaded image G with a width of %dpx, a height of %dpx and %d channels\n", widthg, heightg, channelsg);
    
    //Get content of second image
    unsigned char *img_b = stbi_load(image_2_name, &widthb, &heightb, &channelsb, 0);
    if(img_b == NULL) {
        printf("Error in loading the image\n");
        exit(1);
    }
    printf("Loaded image B with a width of %dpx, a height of %dpx and %d channels\n", widthb, heightb, channelsb);
    
    size_t img_sizeg = widthg * heightg * channelsg;
    size_t img_sizeb = widthb * heightb * channelsb;

    // Loop through pixels to make them grey scale
    int pix = 0;
    int count =0;
    for(unsigned char *pg = img_g, *pb = img_b;  
        pg != img_g + img_sizeg; 
        pg += channelsg, pb += channelsb){

        if ( in_range(pg, pb) || in_range(pg+1, pb+1) || in_range(pg+2, pb+2) ) {
            if (print_content){
                printf("At %d: pg = (%d,%d,%d) pb = (%d,%d,%d) dif: (%d,%d,%d) \n", 
                        pix, 
                        *pg, *(pg+1), *(pg+2), 
                        *pb, *(pb+1), *(pb+2),
                        *pg - *pb, *(pg+1) - *(pb+1), *(pg+2) - *(pb+2));
            }
            count++;
        }
        pix++;
    }
    printf("Count: %d\n",count);
    printf("DEBUG: DONE image_diff\n");
    return 0;
}

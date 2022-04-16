/**
 * @file t0.cpp
 * @author Anthony Meza (you@domain.com)
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

#include "util/colorConv.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image/stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image/stb_image_write.h"

typedef struct {
    float x;
    float y;
} FloatVec;

//DEBUG: FUNCTION USED TO CHECK DIFF BETWEEN PIXELS IN IMAGE
int in_range(unsigned char* p1, unsigned char* p2){
    float value = ( ((float)*p1) - ((float) *p2));
    
    if (abs(value) > 2){
        return 1;
    }
    return 0;
}

//DEBUG: FUNCTION USED TO FIND WHAT PIXELS ARE DIFFRENT BETWEEN TWO IMAGES
/**
 * @brief 
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
}


int main(void) {
    
    int width, height, channels;

    int out_width = 16;
    int out_height = 16;

    float L, a, b;
    unsigned char R, G, B;

    /// TESTING TO SEE WHAT IS DIFFERENT WITH THE TWO OUTPUT IMAGES ///
    
    image_diff("SonicFlower.jpeg" , "SonicFlower_gray.jpeg", 1);
    
    /// TESTING TO SEE WHAT IS DIFFERENT WITH THE TWO OUTPUT IMAGES ///
 
    /// Load the image into img
    unsigned char *img = stbi_load("SonicFlower.jpeg", &width, &height, &channels, 0);
    if(img == NULL) {
        printf("Error in loading the image\n");
        exit(1);
    }
    printf("Loaded image with a width of %dpx, a height of %dpx and %d channels\n", width, height, channels);
    
    ///*** Set up superpixel data structures ***///
    size_t img_size = width * height * channels;
    size_t out_img_size = out_width * out_height * channels;
    FloatVec *superpixel_pos = (FloatVec *) calloc(out_width * out_height, sizeof(FloatVec));
    int *superpixel_img = (int *) calloc(width * height, sizeof(int));
    if (superpixel_pos == NULL) {
        printf("Unable to allocate memory for the cielab image.\n");
        exit(1);
    }

    // initialize superpixel positions
    for (int j = 0; j < out_height; ++j) {
        for (int i = 0; i < out_width; ++i) {
            float x = (i + 0.5f) * width / out_width;
            float y = (j + 0.5f) * height / out_height;
            superpixel_pos[out_width * j + i] = {x, y};
        }
    }

    // initialize superpixel assignments
    for (int j = 0; j < height; ++j) {
        for (int i = 0; i < width; ++i) {
            float dx = (float) width/(float) out_width;
            float dy = (float) height/(float) out_height;
            int x = (int) ((float) i / dx);
            int y = (int) ((float) j / dy);
            superpixel_img[width * j + i] = out_width * y + x;
        }
    }

    float *lab_image = (float *) calloc(img_size, sizeof(float));
    if (lab_image == NULL) {
        printf("Unable to allocate memory for the cielab image.\n");
        exit(1);
    }
    
    unsigned char *p;
    float *pl;

    // Loop through pixels to save them in lab space
    for(p = img, pl = lab_image; p != img + img_size; p += channels, pl += channels) {
        rgb2lab(*p, *(p + 1), *(p + 2), pl, pl + 1, pl + 2);
    }



    // Passed in the correct number of channels, in this case the number desired
    stbi_write_jpg("SonicFlower_output.jpeg", width, height, channels, img, 100);
    

    /// TESTING TO SEE WHAT IS DIFFERENT WITH THE TWO OUTPUT IMAGES ///
    
    image_diff("SonicFlower.jpeg" , "SonicFlower_gray.jpeg", 1);
    
    /// TESTING TO SEE WHAT IS DIFFERENT WITH THE TWO OUTPUT IMAGES ///
 
    stbi_image_free(img);
}
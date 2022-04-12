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

    float L, a, b;
    unsigned char R, G, B;

    /// TESTING TO SEE WHAT IS DIFFERENT WITH THE TWO OUTPUT IMAGES ///
    
    image_diff("SonicFlower.jpeg" , "SonicFlower_gray.jpeg", 1);
    
    printf("\n\n _My test Cases_ \n\n");
    rgb2lab(228, 243, 255, &L, &a, &b);
    lab2rgb(L, a, b, &R, &G, &B);
    printf("(%d, %d, %d)\n\n", R, G, B);

    rgb2lab(47, 95, 159, &L, &a, &b);
    lab2rgb(L, a, b, &R, &G, &B);
    printf("(%d, %d, %d)\n\n", R, G, B);
  
    
    /// TESTING TO SEE WHAT IS DIFFERENT WITH THE TWO OUTPUT IMAGES ///
 
    /// Load the image into img
    unsigned char *img = stbi_load("SonicFlower.jpeg", &width, &height, &channels, 0);
    if(img == NULL) {
        printf("Error in loading the image\n");
        exit(1);
    }
    printf("Loaded image with a width of %dpx, a height of %dpx and %d channels\n", width, height, channels);
    
    ///*** Conversion to Grey image ***///
    size_t img_size = width * height * channels;
    int gray_channels = channels;
    size_t gray_img_size = width * height * gray_channels;

    // Allocate image for grey
    unsigned char *gray_img = (unsigned char *) malloc(gray_img_size);
    if(gray_img == NULL) {
        printf("Unable to allocate memory for the gray image.\n");
        exit(1);
    }
 
    // Loop through pixels to make them grey scale
    for(unsigned char *p = img, *pg = gray_img; p != img + img_size; p += channels, pg += gray_channels) {
        float L, a, b;
        rgb2lab(*p, *(p + 1), *(p + 2), &L, &a, &b);
        lab2rgb(L, a, b, pg, pg + 1, pg + 2);
        if(channels == 4) {
            *(pg + 1) = *(p + 3);
        }
    }

    //float L, a, b;
    //unsigned char R, G, B;
    printf("\n\n _My test Cases_ \n\n");
    rgb2lab(255, 255, 0, &L, &a, &b);
    lab2rgb(L, a, b, &R, &G, &B);
    printf("(%d, %d, %d)\n\n", R, G, B);
    
    rgb2lab(228, 243, 255, &L, &a, &b);
    lab2rgb(L, a, b, &R, &G, &B);
    printf("(%d, %d, %d)\n\n", R, G, B);

    rgb2lab(47, 95, 159, &L, &a, &b);
    lab2rgb(L, a, b, &R, &G, &B);
    printf("(%d, %d, %d)\n\n", R, G, B);

    // Passed in the correct number of channels, in this case the number desired
    stbi_write_jpg("SonicFlower_output.jpeg", width, height, channels, img, 100);
    stbi_write_jpg("SonicFlower_gray.jpeg", width, height, gray_channels, gray_img, 100);

    

    /// TESTING TO SEE WHAT IS DIFFERENT WITH THE TWO OUTPUT IMAGES ///
    
    image_diff("SonicFlower.jpeg" , "SonicFlower_gray.jpeg", 1);
    
    /// TESTING TO SEE WHAT IS DIFFERENT WITH THE TWO OUTPUT IMAGES ///
 
    stbi_image_free(img);
    free(gray_img);
}
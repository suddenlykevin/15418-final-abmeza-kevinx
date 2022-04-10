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

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image/stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image/stb_image_write.h"

int main(void) {
    
    /// Variables for image processing
    int width, height, channels;
    //int desired_no_channels = 4;
    
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

    float L, a, b;
    unsigned char R, G, B;

    rgb2lab(255, 255, 0, &L, &a, &b);
    lab2rgb(L, a, b, &R, &G, &B);
    printf("(%d, %d, %d)\n", R, G, B);
    
    rgb2lab(127, 2, 212, &L, &a, &b);
    lab2rgb(L, a, b, &R, &G, &B);
    printf("(%d, %d, %d)\n", R, G, B);

    rgb2lab(8, 25, 60, &L, &a, &b);
    lab2rgb(L, a, b, &R, &G, &B);
    printf("(%d, %d, %d)\n", R, G, B);

    rgb2lab(5, 100, 215, &L, &a, &b);
    lab2rgb(L, a, b, &R, &G, &B);
    printf("(%d, %d, %d)\n", R, G, B);

    // Passed in the correct number of channels, in this case the number desired
    stbi_write_jpg("SonicFlower_output.jpeg", width, height, channels, img, 100);
    stbi_write_jpg("SonicFlower_gray.jpeg", width, height, gray_channels, gray_img, 100);
 
    stbi_image_free(img);
    free(gray_img);
}
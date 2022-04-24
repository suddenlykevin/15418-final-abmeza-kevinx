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


#include "pixImage.h"

// Import stb_image libraries
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image/stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image/stb_image_write.h"


#include <stdio.h>
#include <stdlib.h>
#include <algorithm> 
#include <stack>



//*************************************************************//
//**************** WRAPPER HELPER FUNCTIONS *******************//
//*************************************************************//

/**
 * @brief Checks to make sure no error occured in stbi_load
 * 
 * @param filename             name of file being used
 * @param x_ptr                pointer of x parameter
 * @param y_ptr                pointer of y parameter
 * @param channels_in_file_ptr pointer to channels in file
 * @param desired_channels     int to desired channels
 * @return stbi_uc* the array of rgb values for image
 */
unsigned char* wrp_stbi_load(char const *filename, int *x_ptr, int *y_ptr, int *channels_in_file_ptr, int desired_channels){
    unsigned char *img = stbi_load(filename, x_ptr, y_ptr, channels_in_file_ptr, desired_channels);
    if(img == NULL) {
        printf("Error in loading the image with name %s\n",filename);
        exit(1);
    }
    return img;
}


//*************************************************************//
//************************ MAIN CODE **************************//
//*************************************************************//

/**
 * @brief 
 * 
 * @return int 
 */
int main(void) {   

    // INITAILZE VARIABLES
    // Input Image 

    int width, height;        //<- width and height of input_img, gives pixel dimensions
    int channels;             //<- number of channels for input_img (3:rgb or 4:rgba)
    unsigned char *input_img;      //<- input image loaded, uses rgb values for pixels (0-255)

    // Output Image 
    int out_width, out_height; //<- output version of width, height   
    //Palette
    int K_colors;          //<- number of colors we aim to use in the pallette
 


    // SET SOME VARIABLES 
    
    //(TODO: Dynamic through cmdline)
    out_width = 32; out_height = 32;
    K_colors = 16;

    //*** ******************* ***//
    //*** PROCESS INPUT IMAGE ***//
    //*** ******************* ***//

    // Load input image (always 3 channels)
    input_img = wrp_stbi_load("SonicFlower.jpeg", &width, &height, &channels, 3);

  

    // Create pixImage object!
    PixImage pixImage(input_img, width, height, out_width, out_height, K_colors);


    //*** ******************** ***//
    //*** (4.1) INITIALIZATION ***//
    //*** ******************** ***//
    pixImage.initialize();

    //*** ******************* ***//
    //*** CORE ALGORITHM LOOP ***//
    //*** ******************* ***//

    pixImage.iterate();



    //*** ******************** ***//
    //*** PROCESS OUTPUT IMAGE ***//
    //*** ******************** ***//
    // Passed in the correct number of channels, in this case the number desired
    stbi_write_jpg("SonicFlower_output.jpeg", pixImage.out_width, pixImage.out_height, channels, pixImage.output_img, 100);
    
    
    //*** ********** ***//
    //*** FREE STUFF ***//
    //*** ********** ***//

    stbi_image_free(input_img);

    pixImage.freeAll();
}




//****************************************************************//
//**************** DEBUG HELPER FUNCTION STUFF *******************//
//****************************************************************//

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
int image_diff(char const *image_1_name, char const *image_2_name, int print_content){
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

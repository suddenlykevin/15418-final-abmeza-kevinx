/**
 * @file t0.cpp
 * @authors Kevin Xie (kevinx) 
 *          Anthony Meza (abmeza)
 * 
 * @brief  main file for out project
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

#include "util/CycleTimer.h"

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <algorithm> 
#include <stack>



//*************************************************************//
//**************** WRAPPER/HELPER FUNCTIONS *******************//
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


/**
 * @brief Prints out how to input values into the command line
 *        Most of code borrowed from 15418 assignments
 */
void usage(const char *progname){
    printf("Usage: %s [options] scenename\n", progname);
    printf("Program Inputs:\n");
    printf("  -f  --filename   <FILENAME>  name of image file being inputted to be pixelated\n");
    printf("  -x  --out_width  <INT>       pixel width of output image\n");
    printf("  -y  --out_height <INT>       pixel height of output image\n");
    printf("  -k  --K_COLORS   <INT>       colors present in the output image\n");
    printf("Program Options:\n");
    printf("  -?  --help                   open this usage help page\n");
}


//*************************************************************//
//************************ MAIN CODE **************************//
//*************************************************************//

int main(int argc, char** argv) {   

    // INITAILZE VARIABLES  
    // Input Image
    int width, height;         //<- width and height of input_img, gives pixel dimensions
    int channels;              //<- number of channels for input_img (3:rgb or 4:rgba)
    unsigned char *input_img;  //<- input image loaded, uses rgb values for pixels (0-255)
    char *filename = NULL;     //<- Holder of input file name (default: NULL)

    // Output Image 
    int out_width = 30;        //<- user provided output width for image (default: 15)
    int out_height = 30;       //<- user provided output height for image (default: 15)  
    char *out_name;            //<- output image file name 
    char *spout_name;          //<- superpixel output image file namee
    char *input_basename;       //<- base of inputimage name used to create output image names

    //Palette
    int K_colors = 8;          //<- number of colors we aim to use in the pallette (default: 8)
     


    // PARSE + CHECK USER INPUT
    int opt; 
    while ((opt = getopt(argc, argv, "f:x:y:k:")) != -1){
        switch (opt) {
        case 'f':
            filename = optarg;
            break;
        case 'x':
            out_width = atoi(optarg);
            break;
        case 'y':
            out_height = atoi(optarg);
            break;
        case 'k':
            K_colors = atoi(optarg);
            break;
    
        case '?':
        default:
            usage(argv[0]);
            return 1;
            break;
        }
    }

    if (filename == NULL) {
        printf("ERROR: No input image was given\n");
        usage(argv[0]);
        return -1;
    }
    
    //*** ******************* ***//
    //*** PROCESS INPUT IMAGE ***//
    //*** ******************* ***//


    // Load input image (always 3 channels)
    input_img = wrp_stbi_load(filename, &width, &height, &channels, 3);
  

    // Create pixImage object!
    PixImage pixImage(input_img, width, height, out_width, out_height, K_colors);


    //*** ********************** ***//
    //*** CREATE PIXELATED IMAGE ***//
    //*** ********************** ***//
    
    pixImage.runPixelate();

    //*** ******************** ***//
    //*** PROCESS OUTPUT IMAGE ***//
    //*** ******************** ***//

    // Create output image names
    input_basename = strtok(filename, ".");
    asprintf(&out_name, "%s_%dx%d_%d.png", input_basename, out_width, out_height, K_colors);
    asprintf(&spout_name, "%s_%dx%d_%d_sp.png", input_basename, out_width, out_height, K_colors);

    // Create the output image files in current directory
    stbi_write_png(out_name, pixImage.out_width, pixImage.out_height, channels, pixImage.output_img, pixImage.out_width * channels);
    stbi_write_png(spout_name, pixImage.in_width, pixImage.in_height, channels, pixImage.spoutput_img, pixImage.in_width * channels);
    

    //*** ********** ***//
    //*** FREE STUFF ***//
    //*** ********** ***//

    stbi_image_free(input_img);
    pixImage.freeAll();
    
}
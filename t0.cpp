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
 * @brief data structure for Lab color scheme
 */
typedef struct {
    float L;   //< L value, represents light
    float a;   //< a value, represents TODO
    float b;   //< b value, repereents TODO
} LabColor;


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


//*********************************************//
//**************** colorConv? *******************//
//**********************************************//


//*********************************************//
//**************** superPixe? *******************//
//**********************************************//


/**
 * @brief Initializes the superPixel_pos array and the superPixel_img array
 * 
 * @param[out] superPixel_pos Array with superpixel center positions on the input image
 * @param[out] superPixel_img holds info which super pixel all input image pixels have
 * @param[in]  in_width       pixel width of input image
 * @param[in]  in_height      pixel height of input image
 * @param[in]  out_width      pixel width of output image
 * @param[in]  out_height     pixel height of output image
 */
void init_superPixels(FloatVec* superPixel_pos, int* superPixel_img, int in_width, int in_height, int out_width, int out_height){
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
            superPixel_img[in_width * j + i] = out_width * y + x;

            // ********** DEBUG ********** //
            //printf("Pixel %d: (%d)\n", in_width * j + i, out_width * y + x);
            // ********** DEBUG ********** //
        }
    }
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
    LabColor *input_img_lab;     //<- input image, using cielab values 
    int M_pix;                //<- # of pixels from the input image (M from paper)

    // Output Image 
    int out_width, out_height; //<- output version of width, height
    unsigned char *output_img;      //<- output version of input_img
    LabColor *output_img_lab;     //<- output version of input_img_lab
    int N_pix;                 //<- # of pixels in the output image (N from paper)

    // Superpixel calculation 
    FloatVec *superPixel_pos; //<- Super pixel coordinate positions "on input image"
    int *superPixel_img;      //<- array with values for pixels assosiated with a specific superpixel
    int m_gerstner = 45;       

    // Palette 
    LabColor *palette_lab;   //<- palette array with 
    int k_count;          //<- Current # of colors stored in palette_lab
    int K_colors;         //<- number of colors we aim to use in the pallette
 
    // Temperature
    int T;   //<- Current temperature
    int T_c; //<- T critical, what will determine convergence and increase in palette (TODO: Edit)
    int T_f; //<- T final, dictates when we finish our core algorithm (TODO: Edit)


    // SET SOME VARIABLES 
    
    //(TODO: Embedd some of these with the MAKE during file generation)
    out_width = 16; out_height = 16;
    
    K_colors = 8;
    k_count = 0;
     
    N_pix = out_width * out_height;
    T_f = 1;    

    //*** ******************* ***//
    //*** PROCESS INPUT IMAGE ***//
    //*** ******************* ***//

    // Load input image
    
    input_img= wrp_stbi_load("SonicFlower.jpeg", &width, &height, &channels, 0);

    // Get M pixels count for input image
    M_pix = width * height;

    // Create input_img_lab version
    input_img_lab = (LabColor *) wrp_calloc(M_pix, sizeof(LabColor));

    unsigned char *p; LabColor *pl;
    for(p = input_img, pl = input_img_lab; p != input_img + (M_pix*channels); p += channels, pl ++) {
        rgb2lab(*p, *(p+1), *(p+2), &(pl->L), &(pl->a), &(pl->b));
    }

    //*** ******************** ***//
    //*** (4.1) INITIALIZATION ***//
    //*** ******************** ***//


    ///*** Allocate Array Space ***///
    superPixel_pos = (FloatVec *) wrp_calloc(N_pix, sizeof(FloatVec)); 
    superPixel_img = (int *) wrp_calloc(M_pix, sizeof(int));
    palette_lab = (LabColor *) wrp_calloc(M_pix *channels, sizeof(float)); \
    output_img = (unsigned char *) wrp_malloc(N_pix * channels); 
    output_img_lab = (LabColor *) wrp_calloc(N_pix, sizeof(LabColor)); 

    ///*** Initialize Superpixel Contents ***///
    init_superPixels(superPixel_pos, superPixel_img, width, height, out_width, out_height);


    ///*** Initialize Palette Contents ***///
    // Find mean of all M_pix color inputs
    LabColor color_sum;
    //add all colors
    for (int p = 0; p < M_pix; p++) {
        color_sum.L += input_img_lab[p].L;
        color_sum.a += input_img_lab[p].a;
        color_sum.b += input_img_lab[p].b;
    }
    // divide all by M_pix
    color_sum.L = color_sum.L / (float) M_pix;
    color_sum.a = color_sum.a / (float) M_pix;
    color_sum.b = color_sum.b / (float) M_pix;
    // update all superpixel color values
    for (int p = 0; p < N_pix; p++) {
        output_img_lab[p].L = color_sum.L;
        output_img_lab[p].a = color_sum.a;
        output_img_lab[p].b = color_sum.b;
    }
    // Store this color
    palette_lab[0] = color_sum;
    k_count ++;



    //*** ******************* ***//
    //*** CORE ALGORITHM LOOP ***//
    //*** ******************* ***//

    // Size of super pixel on input image
    float S = sqrt(((float) (M_pix))/((float) (N_pix)));
    
    // update superpixel segments
    for (int iter = 0; iter < 35; ++iter) {

        // update boundaries
        for (int j = 0; j < out_height; ++j) {
            for (int i = 0; i < out_width; ++i) {
                
                // get local region
                FloatVec center = superPixel_pos[out_width * j + i];
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
                        int curr_spidx = superPixel_img[curr_idx];
                        FloatVec curr_spixel = superPixel_pos[curr_spidx];
                        int curr_spx = (int) round(curr_spixel.x);
                        int curr_spy = (int) round(curr_spixel.y);
                        curr_spidx = curr_spy * width + curr_spx;
                        float dist_curr = dist_k(m_gerstner, S, input_img_lab[curr_spidx].L, 
                                                input_img_lab[curr_spidx].a, input_img_lab[curr_spidx].b, 
                                                curr_spx, curr_spy, input_img_lab[curr_idx].L, input_img_lab[curr_idx].a, 
                                                input_img_lab[curr_idx].b, xx, yy);
                        float dist_new = dist_k(m_gerstner, S, input_img_lab[idx].L, 
                                                input_img_lab[idx].a, input_img_lab[idx].b, 
                                                x, y, input_img_lab[curr_idx].L, input_img_lab[curr_idx].a, 
                                                input_img_lab[curr_idx].b, xx, yy);
                        if (dist_new < dist_curr) {
                            superPixel_img[curr_idx] = out_width*j + i;
                        }
                    }
                }
            }
        }


        // Variables to find mean color values
        FloatVec *sp_sums = (FloatVec *) calloc(N_pix, sizeof(FloatVec));
        float *color_sums = (float *) calloc(N_pix * channels, sizeof(float));
        int *sp_count = (int *) calloc(N_pix, sizeof(int)); 

        // Find the mean colors (from input image) for each superpixel
        for (int j = 0; j < height; j++) {
            for (int i = 0; i < width; i++) {
                int idx = j*width + i;
                int spidx = superPixel_img[idx];
                sp_count[spidx] ++;
                sp_sums[spidx].x += i;
                sp_sums[spidx].y += j;

                color_sums[3*spidx] += input_img_lab[idx].L;
                color_sums[3*spidx + 1] += input_img_lab[idx].a;
                color_sums[3*spidx + 2] += input_img_lab[idx].b;
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
                superPixel_pos[spidx] = newpos;

                // Set output_img_lab to new mean value
                output_img_lab[spidx].L = color_sums[3*spidx]/sp_count[spidx];
                output_img_lab[spidx].a = color_sums[3*spidx + 1]/sp_count[spidx];
                output_img_lab[spidx].b = color_sums[3*spidx + 2]/sp_count[spidx];
                
            }
        }

        // smooth
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
    }

    // for (int j = 0; j < height; j++) {
    //     for (int i = 0; i < width; i++) {
    //         int idx = j*width + i;
    //         int spidx = superPixel_img[idx];
    //         lab2rgb(output_img_lab[3*spidx], output_img_lab[3*spidx + 1], output_img_lab[3*spidx + 2], &output_img[3*idx], &output_img[3*idx + 1], &output_img[3*idx + 2]);
    //         // if ((superPixel_img[idx]/out_width % 2 == 0 && superPixel_img[idx] % 2 == 0) ||
    //         //     (superPixel_img[idx]/out_width % 2 == 1 && superPixel_img[idx] % 2 == 1)) {
    //         //     output_img[idx*3] = 0;
    //         //     output_img[idx*3 + 1] = 0;
    //         //     output_img[idx*3 + 2] = 0;
    //         // } else {
    //         //     output_img[idx*3] = 255;
    //         //     output_img[idx*3 + 1] = 255;
    //         //     output_img[idx*3 + 2] = 255;
    //         // }
    //     }
    // }

    // Create output image in rgb color values
    for (int j = 0; j < out_height; j++) {
        for (int i = 0; i < out_width; i++) {
            int idx = j*out_width + i;
            lab2rgb(output_img_lab[idx].L, output_img_lab[idx].a, output_img_lab[idx].b, 
                    &(output_img[3*idx]), &(output_img[3*idx + 1]), &(output_img[3*idx + 2]));
        }
    }

    // Passed in the correct number of channels, in this case the number desired
    stbi_write_jpg("SonicFlower_output.jpeg", out_width, out_height, channels, output_img, 100);
    

    /// TESTING TO SEE WHAT IS DIFFERENT WITH THE TWO OUTPUT IMAGES ///
    
    // image_diff("SonicFlower.jpeg" , "SonicFlower_gray.jpeg", 1);
    
    /// TESTING TO SEE WHAT IS DIFFERENT WITH THE TWO OUTPUT IMAGES ///
 
    // FREE STUFF
    stbi_image_free(input_img);
    free(input_img_lab);

    free(output_img);
    free(output_img_lab);
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

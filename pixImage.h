/**
 * @file pixImage.h
 * @author Kevin Xie () 
 *          Anthony Meza (abmeza)
 * @brief Header file with the contents of the pixImage object
 * 
 * @version 0.1
 * @date 2022-04-22
 * 
 * @copyright Copyright (c) 2022
 * 
 */



// GLOBAL CONSTANTS
const float T_f = 1.0f; //<- T final, dictates when we finish our core algorithm (TODO: Edit)
const int m_gerstner = 45; //<- 
const float kSubclusterPertubation = .8f; //<- amount to perturb subcluster of each palette color
const float kT0SafetyFactor = 1.1f; //<- initial temperature is 1.1*T_c
const float kPaletteErrorTolerance = 1.0f;
const float kTF = 1.0f;
const float kDT = .7f;
const float kSubclusterTolerance = 1.6f;
const int maxIter = 1000;

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

typedef struct {
    int a;
    int b;
} PalettePair;

class PixImage{
  private:
    // Timing variables for our reference
    double startAllTime, endAllTime;
    double startInitializeTime, endInitializeTime;

  public:
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
    int *superPixel; //<- Map of superpixels
    int *region_map;      //<- array with values for pixels associated with a specific superpixel 

    // Palette 
    int K_colors;          //<- number of colors we aim to use in the pallette
    int palette_size;      //<- Current # of colors stored in palette_lab
    PalettePair *palette_pairs;
    int *palette_assign; //<- palette assignment for each superpixel
    LabColor *palette_lab; //<- palette array with color values in palette
    bool palette_complete;

    float *prob_c;         //<- array of probabiities that a color in the palette is set to ANY super pixel
    float prob_sp;        //<- array of probabiities of each super pixel TODO:EDIT
    float *prob_c_if_sp;   //<- List of P(c_k|p_s) values for all superpixels
    
    // Temperature
    float T;   //<- Current temperature
    
    /**
     * @brief Construct a new Pix Image object
     * 
     * @param input_image name of input image being used
     * @param in_w        width of input image
     * @param in_h        height of input image
     * @param out_w       desiered width of output image
     * @param out_h       desired height of output image
     * @param K           number of colors desired for output image
     */
    PixImage(unsigned char* input_image, int in_w, int in_h, int out_w, int out_h, int K);

    /**
     * @brief Initializes all values, run in initialize step
     */
    void initialize();


    /**
     * @brief 
     * 
     */
    void initSuperPixels();

    void updateSuperPixelMeans();

    void pushPaletteColor(LabColor color, float prob);

    void pushPalettePair(int a, int b);

    void getAveragedPalette(LabColor *avg_palette);

    /**
     * @brief Runs the entire algorithm that initializes variables and
     *        pixelates a given input image
     */
    void runPixelate();

    void getMajorAxis(int palette_index, float *value, LabColor *vector);

    void splitColor(int pair_index);

    void condensePalette();

    /**
     * @brief Frees everything that has been alocated after running runPixelate 
     */
    void freeAll();

};
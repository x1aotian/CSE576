#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "image.h"
#include "matrix.h"
#include <time.h>

// ADD helper: print image for debug
void print_image(image im){
    for (int x=0; x<im.w; x++){
        for (int y=0; y<im.h; y++){
            printf("%f ", get_pixel(im, x, y, 0));
        }
        printf("\n");
    }
}

// Frees an array of descriptors.
// descriptor *d: the array.
// int n: number of elements in array.
void free_descriptors(descriptor *d, int n)
{
    int i;
    for(i = 0; i < n; ++i){
        free(d[i].data);
    }
    free(d);
}

// Create a feature descriptor for an index in an image.
// image im: source image.
// int i: index in image for the pixel we want to describe.
// returns: descriptor for that index.
descriptor describe_index(image im, int i)
{
    int w = 5;
    descriptor d;
    d.p.x = i%im.w;
    d.p.y = i/im.w;
    d.data = calloc(w*w*im.c, sizeof(float));
    d.n = w*w*im.c;
    int c, dx, dy;
    int count = 0;
    // If you want you can experiment with other descriptors
    // This subtracts the central value from neighbors
    // to compensate some for exposure/lighting changes.
    for(c = 0; c < im.c; ++c){
        float cval = im.data[c*im.w*im.h + i];
        for(dx = -w/2; dx < (w+1)/2; ++dx){
            for(dy = -w/2; dy < (w+1)/2; ++dy){
                float val = get_pixel(im, i%im.w+dx, i/im.w+dy, c);
                d.data[count++] = cval - val;
            }
        }
    }
    return d;
}

// Marks the spot of a point in an image.
// image im: image to mark.
// ponit p: spot to mark in the image.
void mark_spot(image im, point p)
{
    int x = p.x;
    int y = p.y;
    int i;
    for(i = -9; i < 10; ++i){
        set_pixel(im, x+i, y, 0, 1);
        set_pixel(im, x, y+i, 0, 1);
        set_pixel(im, x+i, y, 1, 0);
        set_pixel(im, x, y+i, 1, 0);
        set_pixel(im, x+i, y, 2, 1);
        set_pixel(im, x, y+i, 2, 1);
    }
}

// Marks corners denoted by an array of descriptors.
// image im: image to mark.
// descriptor *d: corners in the image.
// int n: number of descriptors to mark.
void mark_corners(image im, descriptor *d, int n)
{
    int i;
    for(i = 0; i < n; ++i){
        mark_spot(im, d[i].p);
    }
}

// Creates a 1d Gaussian filter.
// float sigma: standard deviation of Gaussian.
// returns: single row image of the filter.
image make_1d_gaussian(float sigma)
{
    // TODO: make separable 1d Gaussian.
    int s = ceilf(6 * sigma);
    if (s%2 == 0) s++;
    float val = 0;
    image filter = make_image(s, 1, 1);
    for (int x=0; x<filter.w; x++){
        int x_ = x - s/2;
        val = 1. /sqrt(TWOPI) /sigma * exp( -pow(x_, 2) /2. /pow(sigma, 2) );
        set_pixel(filter, x, 0, 0, val);
    }
    l1_normalize(filter);
    return filter;
}

// ADD helper: Transpose image
image transpose_image(image im){
    image im_t = make_image(im.h, im.w, im.c);
    for (int x=0; x<im.w; x++){
        for (int y=0; y<im.h; y++){
            for (int c=0; c<im.c; c++){
                set_pixel(im_t, y, x, c, get_pixel(im, x, y, c));
            }
        }
    }
    return im_t;
}

// Smooths an image using separable Gaussian filter.
// image im: image to smooth.
// float sigma: std dev. for Gaussian.
// returns: smoothed image.
image smooth_image(image im, float sigma)
{
    // TODO: use two convolutions with 1d gaussian filter.
    // image gaussian_2d = make_gaussian_filter(1);
    // image img_smooth = convolve_image(im, gaussian_2d, 1);
    image f_g1d = make_1d_gaussian(sigma);
    image f_g1d_vertical = transpose_image(f_g1d);
    image img_smooth = convolve_image(im, f_g1d, 1);
    img_smooth = convolve_image(img_smooth, f_g1d_vertical, 1);
    free_image(f_g1d);
    free_image(f_g1d_vertical);
    return img_smooth;
}

// Calculate the structure matrix of an image.
// image im: the input image.
// float sigma: std dev. to use for weighted sum.
// returns: structure matrix. 1st channel is Ix^2, 2nd channel is Iy^2,
//          third channel is IxIy.
image structure_matrix(image im, float sigma)
{
    image S = make_image(im.w, im.h, 3);
    // TODO: calculate structure matrix for im.
    image f_gx = make_gx_filter();
    image f_gy = make_gy_filter();
    image img_gx = convolve_image(im, f_gx, 0);
    image img_gy = convolve_image(im, f_gy, 0);
    for (int x=0; x<im.w; x++){
        for (int y=0; y<im.h; y++){
            set_pixel(S, x, y, 0, pow(get_pixel(img_gx, x, y, 0), 2));
            set_pixel(S, x, y, 1, pow(get_pixel(img_gy, x, y, 0), 2));
            set_pixel(S, x, y, 2, get_pixel(img_gx, x, y, 0) * get_pixel(img_gy, x, y, 0));
        }
    }
    S = smooth_image(S, sigma);
    free_image(f_gx);
    free_image(f_gy);
    free_image(img_gx);
    free_image(img_gy);
    return S;
}

// Estimate the cornerness of each pixel given a structure matrix S.
// image S: structure matrix for an image.
// returns: a response map of cornerness calculations.
image cornerness_response(image S)
{
    image R = make_image(S.w, S.h, 1);
    // TODO: fill in R, "cornerness" for each pixel using the structure matrix.
    // We'll use formulation det(S) - alpha * trace(S)^2, alpha = .06.
    float det, trace;
    float alpha = .06;
    for (int x=0; x<S.w; x++){
        for (int y=0; y<S.h; y++){
            det = get_pixel(S, x, y, 0) * get_pixel(S, x, y, 1) - pow(get_pixel(S, x, y, 2), 2);
            trace = get_pixel(S, x, y, 0) + get_pixel(S, x, y, 1);
            set_pixel(R, x, y, 0, det - alpha * trace * trace);
        }
    }
    return R;
}

// Perform non-max supression on an image of feature responses.
// image im: 1-channel image of feature responses.
// int w: distance to look for larger responses.
// returns: image with only local-maxima responses within w pixels.
image nms_image(image im, int w)
{
    image r = copy_image(im);
    // TODO: perform NMS on the response map.
    // for every pixel in the image:
    //     for neighbors within w:
    //         if neighbor response greater than pixel response:
    //             set response to be very low (I use -999999 [why not 0??])
    float low = -999999;
    float p_val;
    for (int x=0; x<r.w; x++){
        for (int y=0; y<r.h; y++){
            p_val = get_pixel(im, x, y, 0);
            for (int i=x-w; i<=x+w; i++){
                for (int j=y-w; j<=y+2; j++){
                    if (get_pixel(r, i, j, 0) > p_val) set_pixel(r, x, y, 0, low);
                }
            }
        }
    }
    return r;
}

// Perform harris corner detection and extract features from the corners.
// image im: input image.
// float sigma: std. dev for harris.
// float thresh: threshold for cornerness.
// int nms: distance to look for local-maxes in response map.
// int *n: pointer to number of corners detected, should fill in.
// returns: array of descriptors of the corners in the image.
descriptor *harris_corner_detector(image im, float sigma, float thresh, int nms, int *n)
{
    // Calculate structure matrix
    image S = structure_matrix(im, sigma);

    // Estimate cornerness
    image R = cornerness_response(S);

    // Run NMS on the responses
    image Rnms = nms_image(R, nms);


    //TODO: count number of responses over threshold
    int count = 1; // change this
    for (int x=0; x<Rnms.w; x++) {
        for (int y=0; y<Rnms.h; y++) {
            if (get_pixel(Rnms, x, y, 0) > thresh) count++;
        }
    }
    
    *n = count; // <- set *n equal to number of corners in image.
    descriptor *d = calloc(count, sizeof(descriptor));
    //TODO: fill in array *d with descriptors of corners, use describe_index.
    int ptr = 0;
    for (int x=0; x<Rnms.w; x++) {
        for (int y=0; y<Rnms.h; y++) {
            if (get_pixel(Rnms, x, y, 0) > thresh) d[ptr++] = describe_index(im, x + y*im.w);
        }
    }

    free_image(S);
    free_image(R);
    free_image(Rnms);
    return d;
}

// Find and draw corners on an image.
// image im: input image.
// float sigma: std. dev for harris.
// float thresh: threshold for cornerness.
// int nms: distance to look for local-maxes in response map.
void detect_and_draw_corners(image im, float sigma, float thresh, int nms)
{
    int n = 0;
    descriptor *d = harris_corner_detector(im, sigma, thresh, nms, &n);
    mark_corners(im, d, n);
}

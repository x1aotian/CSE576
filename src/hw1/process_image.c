// Xiaotian Fang, Apr 13

#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "image.h"

float get_pixel(image im, int x, int y, int c)
{
    // TODO Fill this in
    x = fmax(x, 0);
    x = fmin(x, im.w-1);
    y = fmax(y, 0);
    y = fmin(y, im.h-1);
    c = fmax(c, 0);
    c = fmin(c, im.c-1);
    return im.data[c * im.w * im.h + y * im.w + x];
}

void set_pixel(image im, int x, int y, int c, float v)
{
    // TODO Fill this in
    if (x>=0 && y>=0 && c>=0 && x<im.w && y<im.h && c<im.c)
        im.data[c * im.w * im.h + y * im.w + x] = v;
}

image copy_image(image im)
{
    image copy = make_image(im.w, im.h, im.c);
    // TODO Fill this in
    memcpy(copy.data, im.data, im.w*im.h*im.c * sizeof(float));
    return copy;
}

image rgb_to_grayscale(image im)
{
    assert(im.c == 3);
    image gray = make_image(im.w, im.h, 1);
    // TODO Fill this in
    for (int x=0; x<im.w; x++){
        for (int y=0; y<im.h; y++){
            gray.data[y*im.w+x] = 0.299*im.data[y*im.w+x] + 0.587*im.data[im.w*im.h+y*im.w+x] + .114*im.data[2*im.w*im.h+y*im.w+x];
        }
    }
    return gray;
}

void shift_image(image im, int c, float v)
{
    // TODO Fill this in
    for (int x=0; x<im.w; x++){
        for (int y=0; y<im.h; y++){
            im.data[c*im.w*im.h + y*im.w + x] += v;
        }
    }
}

void clamp_image(image im)
{
    // TODO Fill this in
    for (int x=0; x<im.w; x++){
        for (int y=0; y<im.h; y++){
            for (int c=0; c<im.c; c++){
                if (im.data[c*im.w*im.h + y*im.w + x] > 1)
                    im.data[c*im.w*im.h + y*im.w + x] = 1;
                if (im.data[c*im.w*im.h + y*im.w + x] < 0)
                    im.data[c*im.w*im.h + y*im.w + x] = 0;
            }
        }
    }
}


// These might be handy
float three_way_max(float a, float b, float c)
{
    return (a > b) ? ( (a > c) ? a : c) : ( (b > c) ? b : c) ;
}

float three_way_min(float a, float b, float c)
{
    return (a < b) ? ( (a < c) ? a : c) : ( (b < c) ? b : c) ;
}

void rgb_to_hsv(image im)
{
    // TODO Fill this in
    for (int x=0; x<im.w; x++){
        for (int y=0; y<im.h; y++){
            float R = get_pixel(im, x, y, 0);
            float G = get_pixel(im, x, y, 1);
            float B = get_pixel(im, x, y, 2);
            float V = three_way_max(R, G, B);
            float m = three_way_min(R, G, B);
            float C = V - m;
            float S = 0;
            if (V > 0) S = C / V;
            float H_ = 0, H = 0;
            if (C != 0) {
                if (V == R) {
                    H_ = (G-B)/C;
                } else if (V == G) {
                    H_ = (B-R)/C + 2;
                } else {
                    H_ = (R-G)/C + 4;
                }
                if (H_ < 0) {
                    H = H_/6 + 1;
                } else {
                    H = H_/6;
                }
            }
            set_pixel(im, x, y, 0, H);
            set_pixel(im, x, y, 1, S);
            set_pixel(im, x, y, 2, V);
        }
    }
}

void hsv_to_rgb(image im)
{
    // TODO Fill this in
    for (int x=0; x<im.w; x++){
        for (int y=0; y<im.h; y++){
            float H = get_pixel(im, x, y, 0) * 360.;
            float S = get_pixel(im, x, y, 1);
            float V = get_pixel(im, x, y, 2);
            float C = V * S;
            float X = C * (1 - fabs(fmod(H/60., 2)-1));
            float m = V - C;
            float R = 0., G = 0., B = 0.;
            H /= 60.;
            if (H >= 0 && H < 1){
                R = C; G = X;
            } else if (H >= 1 && H < 2){
                R = X; G = C;
            } else if (H >= 2 && H < 3){
                G = C; B = X;
            } else if (H >= 3 && H < 4){
                G = X; B = C;
            } else if (H >= 4 && H < 5){
                R = X; B = C;
            } else {
                R = C; B = X;
            }
            R += m; G += m; B += m;
            set_pixel(im, x, y, 0, R);
            set_pixel(im, x, y, 1, G);
            set_pixel(im, x, y, 2, B);
        }
    }
}

// Extra credit
void scale_image(image im, int c, float v) {
    for (int x=0; x<im.w; x++){
        for (int y=0; y<im.h; y++){
            float a = get_pixel(im, x, y, c) * v;
            set_pixel(im, x, y, c, a);
        }
    }
}
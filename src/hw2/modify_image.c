#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "image.h"
#define TWOPI 6.2831853

/******************************** Resizing *****************************
  To resize we'll need some interpolation methods and a function to create
  a new image and fill it in with our interpolation methods.
************************************************************************/

float nn_interpolate(image im, float x, float y, int c)
{
    // TODO
    /***********************************************************************
      This function performs nearest-neighbor interpolation on image "im"
      given a floating column value "x", row value "y" and integer channel "c",
      and returns the interpolated value.
    ************************************************************************/
    float x_nn = round(x);
    float y_nn = round(y);

    return get_pixel(im,x_nn,y_nn,c);
}

image nn_resize(image im, int w, int h)
{
    // TODO Fill in (also fix the return line)
    /***********************************************************************
      This function uses nearest-neighbor interpolation on image "im" to a new
      image of size "w x h"
    ************************************************************************/
    image im_new = make_image(w, h, im.c);
    float w_scale = (float) im.w / w;
    float h_scale = (float) im.h / h;

    for (int x = 0; x < w; x++){
      for (int y = 0; y < h; y++){
        for (int ci = 0; ci < im.c; ci++){
          float x_ori = x * w_scale + 0.5*(w_scale -1);
          float y_ori = y * h_scale + 0.5*(h_scale -1);
          float val_ori = nn_interpolate(im, x_ori, y_ori, ci);
          set_pixel(im_new, x, y, ci, val_ori);
        }
      }
    }
    return im_new;
}

float bilinear_interpolate(image im, float x, float y, int c)
{
    // TODO
    /***********************************************************************
      This function performs bilinear interpolation on image "im" given
      a floating column value "x", row value "y" and integer channel "c".
      It interpolates and returns the interpolated value.
    ************************************************************************/
    float x1 = floorf(x);
    float x2 = ceilf(x);
    float y1 = floorf(y);
    float y2 = ceilf(y);
    float fq11 = get_pixel(im, x1, y1, c);
    float fq12 = get_pixel(im, x1, y2, c);
    float fq21 = get_pixel(im, x2, y1, c);
    float fq22 = get_pixel(im, x2, y2, c);
    float fxy1 = (x2-x)/(x2-x1) * fq11 + (x-x1)/(x2-x1) * fq21;
    float fxy2 = (x2-x)/(x2-x1) * fq12 + (x-x1)/(x2-x1) * fq22;
    float fxy = (y2-y)/(y2-y1) * fxy1 + (y-y1)/(y2-y1) * fxy2;
    return fxy;
}

image bilinear_resize(image im, int w, int h)
{
    // TODO
    /***********************************************************************
      This function uses bilinear interpolation on image "im" to a new image
      of size "w x h". Algorithm is same as nearest-neighbor interpolation.
    ************************************************************************/
    image im_new = make_image(w, h, im.c);
    float w_scale = (float) im.w / w;
    float h_scale = (float) im.h / h;

    for (int x = 0; x < w; x++){
      for (int y = 0; y < h; y++){
        for (int ci = 0; ci < im.c; ci++){
          float x_ori = x * w_scale + 0.5*(w_scale -1);
          float y_ori = y * h_scale + 0.5*(h_scale -1);
          float val_ori = bilinear_interpolate(im, x_ori, y_ori, ci);
          set_pixel(im_new, x, y, ci, val_ori);
        }
      }
    }
    return im_new;
}


/********************** Filtering: Box filter ***************************
  We want to create a box filter. We will only use square box filters.
************************************************************************/

void l1_normalize(image im)
{
    // TODO
    /***********************************************************************
      This function divides each value in image "im" by the sum of all the
      values in the image and modifies the image in place.
    ************************************************************************/
    float sum[3] = {0, 0, 0};
    for (int x = 0; x < im.w; x++){
      for (int y = 0; y < im.h; y++){
        for (int c = 0; c < im.c; c++){
          sum[c] += get_pixel(im, x, y, c);
        }
      }
    }
    for (int x = 0; x < im.w; x++){
      for (int y = 0; y < im.h; y++){
        for (int c = 0; c < im.c; c++){
          set_pixel(im, x, y, c, get_pixel(im, x, y, c)/sum[c]);
        }
      }
    }
}

image make_box_filter(int w)
{
    // TODO
    /***********************************************************************
      This function makes a square filter of size "w x w". Make an image of
      width = height = w and number of channels = 1, with all entries equal
      to 1. Then use "l1_normalize" to normalize your filter.
    ************************************************************************/
    image im = make_image(w, w, 1);
    for (int x = 0; x < im.w; x++){
      for (int y = 0; y < im.h; y++){
        set_pixel(im, x, y, 0, 1);
      }
    }
    l1_normalize(im);
    return im;
}

image convolve_image(image im, image filter, int preserve)
{
    // TODO
    /***********************************************************************
      This function convolves the image "im" with the "filter". The value
      of preserve is 1 if the number of input image channels need to be 
      preserved. Check the detailed algorithm given in the README.  
    ************************************************************************/
    assert (filter.c==1 || filter.c==im.c);

    image im_new = make_image(im.w, im.h, im.c);
    for (int x = 0; x < im.w; x++){
      for (int y = 0; y < im.h; y++){
        for (int ci = 0; ci < im.c; ci++){
          int ci_f = filter.c==0 ? 0:ci;
          float new_val = 0.;
          for (int i = 0; i < filter.w; i++){
            for (int j = 0; j < filter.h; j++){
              new_val += get_pixel(im, x+i-filter.w/2, y+j-filter.h/2, ci) * get_pixel(filter, i, j, ci_f);
            }
          }
          set_pixel(im_new, x, y, ci, new_val);
        }
      }
    }

    if (preserve == 0){
      image im_new_p = make_image(im.w, im.h, 1);
      for (int x = 0; x < im_new.w; x++){
        for (int y = 0; y < im_new.h; y++){
          float new_val = 0.;
          for (int ci = 0; ci < im_new.c; ci++){
            new_val += get_pixel(im_new, x, y, ci);
          }
          set_pixel(im_new_p, x, y, 0, new_val);
        }
      }
      return im_new_p;
    }

    return im_new;
    
}

image make_highpass_filter()
{
    // TODO
    /***********************************************************************
      Create a 3x3 filter with highpass filter values using image.data[]
    ************************************************************************/
    image filter = make_image(3, 3, 1);
    filter.data[0] = 0;
    filter.data[1] = -1;
    filter.data[2] = 0;
    filter.data[3] = -1;
    filter.data[4] = 4;
    filter.data[5] = -1;
    filter.data[6] = 0;
    filter.data[7] = -1;
    filter.data[8] = 0;
    return filter;
}

image make_sharpen_filter()
{
    // TODO
    /***********************************************************************
      Create a 3x3 filter with sharpen filter values using image.data[]
    ************************************************************************/
    image filter = make_image(3, 3, 1);
    filter.data[0] = 0;
    filter.data[1] = -1;
    filter.data[2] = 0;
    filter.data[3] = -1;
    filter.data[4] = 5;
    filter.data[5] = -1;
    filter.data[6] = 0;
    filter.data[7] = -1;
    filter.data[8] = 0;
    return filter;
}

image make_emboss_filter()
{
    // TODO
    /***********************************************************************
      Create a 3x3 filter with emboss filter values using image.data[]
    ************************************************************************/
    image filter = make_image(3, 3, 1);
    filter.data[0] = -2;
    filter.data[1] = -1;
    filter.data[2] = 0;
    filter.data[3] = -1;
    filter.data[4] = 1;
    filter.data[5] = 1;
    filter.data[6] = 0;
    filter.data[7] = 1;
    filter.data[8] = 2;
    return filter;
}

// Question 2.3.1: Which of these filters should we use preserve when we run our convolution and which ones should we not? Why?
// Answer: TODO
// High-pass filter is the only filters that we don't need to use preserve.
// Because the high-pass filter is used to find edges. Color information is not required.
// Box filter, sharpen filter, and enboss filter need to be (or can be) preserved.
// Box filter blurs images. Sharpen filter makes imgaes sharper. Emboss filter makes images embossing style.
// All of them needs color information (all channels) preserved.

// Question 2.3.2: Do we have to do any post-processing for the above filters? Which ones and why?
// Answer: TODO
// high-pass filter needs post-processing most.
// Beacuse imgae after high-pass filter is in grayscale. We can do clamping to prevent overflow of pixel values.
// We can also do post-processing to other filters, but not highly needed.

image make_gaussian_filter(float sigma)
{
    // TODO
    /***********************************************************************
      sigma: a float number for the Gaussian.
      Create a Gaussian filter with the given sigma. Note that the kernel size 
      is the next highest odd integer from 6 x sigma. Return the Gaussian filter.
    ************************************************************************/
    int s = ceilf(6 * sigma);
    if (s%2 == 0) s++;
    float val = 0;
    image filter = make_image(s, s, 1);
    for (int x=0; x<filter.w; x++){
      for (int y=0; y<filter.h; y++){
        int x_ = x - s/2;
        int y_ = y - s/2;
        val = 1. /TWOPI /pow(sigma, 2) * exp( -(pow(x_, 2)+pow(y_, 2)) /2. /pow(sigma, 2) );
        set_pixel(filter, x, y, 0, val);
      }
    }
    l1_normalize(filter);
    return filter;
}

image add_image(image a, image b)
{
    // TODO
    /***********************************************************************
      The input images a and image b have the same height, width, and channels.
      Sum the given two images and return the result, which should also have
      the same height, width, and channels as the inputs. Do necessary checks.
    ************************************************************************/
    assert (a.w==b.w && a.h==b.h && a.c==b.c);
    image im_add = make_image(a.w, a.h, a.c);
    for (int x = 0; x < im_add.w; x++){
      for (int y = 0; y < im_add.h; y++){
        for (int c = 0; c < im_add.c; c++){
          set_pixel(im_add, x, y, c, get_pixel(a, x, y, c) + get_pixel(b, x, y, c));
        }
      }
    }
    return im_add;
}

image sub_image(image a, image b)
{
    // TODO
    /***********************************************************************
      The input image a and image b have the same height, width, and channels.
      Subtract the given two images and return the result, which should have
      the same height, width, and channels as the inputs. Do necessary checks.
    ************************************************************************/
    assert (a.w==b.w && a.h==b.h && a.c==b.c);
    image im_sub = make_image(a.w, a.h, a.c);
    for (int x = 0; x < im_sub.w; x++){
      for (int y = 0; y < im_sub.h; y++){
        for (int c = 0; c < im_sub.c; c++){
          set_pixel(im_sub, x, y, c, get_pixel(a, x, y, c) - get_pixel(b, x, y, c));
        }
      }
    }
    return im_sub;
}

image make_gx_filter()
{
    // TODO
    /***********************************************************************
      Create a 3x3 Sobel Gx filter and return it
    ************************************************************************/
    image filter = make_image(3, 3, 1);
    filter.data[0] = -1;
    filter.data[1] = 0;
    filter.data[2] = 1;
    filter.data[3] = -2;
    filter.data[4] = 0;
    filter.data[5] = 2;
    filter.data[6] = -1;
    filter.data[7] = 0;
    filter.data[8] = 1;
    return filter;
}

image make_gy_filter()
{
    // TODO
    /***********************************************************************
      Create a 3x3 Sobel Gy filter and return it
    ************************************************************************/
    image filter = make_image(3, 3, 1);
    filter.data[0] = -1;
    filter.data[1] = -2;
    filter.data[2] = -1;
    filter.data[3] = 0;
    filter.data[4] = 0;
    filter.data[5] = 0;
    filter.data[6] = 1;
    filter.data[7] = 2;
    filter.data[8] = 1;
    return filter;
}

void feature_normalize(image im)
{
    // TODO
    /***********************************************************************
      Calculate minimum and maximum pixel values. Normalize the image by
      subtracting the minimum and dividing by the max-min difference.
    ************************************************************************/
    float min_ = im.data[0];
    float max_ = im.data[0];
    for (int x = 0; x < im.w; x++){
      for (int y = 0; y < im.h; y++){
        for (int c = 0; c < im.c; c++){
          min_ = get_pixel(im, x, y, c) < min_ ? get_pixel(im, x, y, c) : min_;
          max_ = get_pixel(im, x, y, c) > max_ ? get_pixel(im, x, y, c) : max_;
        }
      }
    }
    for (int x = 0; x < im.w; x++){
      for (int y = 0; y < im.h; y++){
        for (int c = 0; c < im.c; c++){
          set_pixel(im, x, y, c, get_pixel(im, x, y, c)/(max_-min_));
        }
      }
    }
}

image *sobel_image(image im)
{
    // TODO
    /***********************************************************************
      Apply Sobel filter to the input image "im", get the magnitude as sobelimg[0]
      and gradient as sobelimg[1], and return the result.
    ************************************************************************/
    image *sobelimg = calloc(2, sizeof(image));
    image im_gx = convolve_image(im, make_gx_filter(), 0);
    image im_gy = convolve_image(im, make_gy_filter(), 0);
    image grad_mag = make_image(im.w, im.h, 1);
    image grad_dir = make_image(im.w, im.h, 1);
    for (int x = 0; x < im.w; x++){
      for (int y = 0; y < im.h; y++){
          set_pixel(grad_mag, x, y, 0, sqrtf(pow(get_pixel(im_gx, x, y, 0), 2) + pow(get_pixel(im_gy, x, y, 0), 2)));
          set_pixel(grad_dir, x, y, 0, atan2f(get_pixel(im_gy, x, y, 0), get_pixel(im_gx, x, y, 0)));
        }
      }
    *sobelimg = grad_mag;
    *(sobelimg+1) = grad_dir;
    return sobelimg;
}

image colorize_sobel(image im)
{
  // TODO
  /***********************************************************************
    Create a colorized version of the edges in image "im" using the 
    algorithm described in the README.
  ************************************************************************/
  image *sobel_im = sobel_image(im);
  image grad_mag = sobel_im[0];
  image grad_dir = sobel_im[1];
  feature_normalize(grad_mag);
  feature_normalize(grad_dir);
  image im_new = make_image(im.w, im.h, 3);
  for (int x = 0; x < im.w; x++){
      for (int y = 0; y < im.h; y++){
        set_pixel(im_new, x, y, 0, get_pixel(grad_dir, x, y, 0));
        set_pixel(im_new, x, y, 1, get_pixel(grad_mag, x, y, 0));
        set_pixel(im_new, x, y, 2, get_pixel(grad_mag, x, y, 0));
      }
  }
  hsv_to_rgb(im_new);
  im_new = convolve_image(im_new, make_gaussian_filter(2), 1);
  return im_new;
}

// EXTRA CREDIT: Median filter

/*
image apply_median_filter(image im, int kernel_size)
{
  return make_image(1,1,1);
}
*/

// SUPER EXTRA CREDIT: Bilateral filter

/*
image apply_bilateral_filter(image im, float sigma1, float sigma2)
{
  return make_image(1,1,1);
}
*/

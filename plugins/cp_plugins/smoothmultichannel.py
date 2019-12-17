# coding=utf-8

"""
SmoothMultichannel
======

**SmoothMultichannel** smooths (i.e., blurs) images.

This module allows you to smooth (blur) images, which can be helpful to
remove small artifacts. Note that smoothing can be a time-consuming process
and that all channels of the image are smoothed individually.

|

============ ============ ===============
Supports 2D? Supports 3D? Respects masks?
============ ============ ===============
YES          NO           YES
============ ============ ===============

See also
^^^^^^^^

See also **Smooth** and several related modules in the *Advanced* category
(e.g., **MedianFilter** and **GaussianFilter**).
"""

import numpy as np
import scipy.ndimage as scind
from centrosome.filter import median_filter, bilateral_filter, circular_average_filter
from centrosome.smooth import circular_gaussian_kernel
from centrosome.smooth import fit_polynomial
from centrosome.smooth import smooth_with_function_and_mask
import skimage.restoration
from scipy import ndimage

import cellprofiler.image as cpi
import cellprofiler.module as cpm
import cellprofiler.setting as cps
from cellprofiler.modules._help import HELP_ON_MEASURING_DISTANCES, HELP_ON_PIXEL_INTENSITIES
from cellprofiler.setting import YES, NO

FIT_POLYNOMIAL = 'Fit Polynomial'
MEDIAN_FILTER = 'Median Filter'
MEDIAN_FILTER_SCIPY = 'Median Filter Scipy'
GAUSSIAN_FILTER = 'Gaussian Filter'
SMOOTH_KEEPING_EDGES = 'Smooth Keeping Edges'
CIRCULAR_AVERAGE_FILTER = 'Circular Average Filter'
SM_TO_AVERAGE = "Smooth to Average"
CLIP_HOT_PIXELS = "Remove single hot pixels"


class SmoothMultichannel(cpm.Module):
    module_name = 'Smooth Multichannel'
    category = "Image Processing"
    variable_revision_number = 5

    def create_settings(self):
        self.image_name = cps.ImageNameSubscriber('Select the input image', cps.NONE, doc="""Select the image to be smoothed.""")

        self.filtered_image_name = cps.ImageNameProvider('Name the output image', 'FilteredImage', doc="""Enter a name for the resulting image.""")

        self.smoothing_method = cps.Choice(
                'Select smoothing method',
                [CLIP_HOT_PIXELS, FIT_POLYNOMIAL, GAUSSIAN_FILTER, MEDIAN_FILTER, MEDIAN_FILTER_SCIPY, SMOOTH_KEEPING_EDGES,
                 CIRCULAR_AVERAGE_FILTER, SM_TO_AVERAGE], doc="""\
This module smooths images using one of several filters. Fitting a
polynomial is fastest but does not allow a very tight fit compared to
the other methods:

-  *%(FIT_POLYNOMIAL)s:* This method is fastest but does not allow
   a very tight “fit” compared to the other methods. Thus, it will usually be less
   accurate. The method treats the intensity of the image
   pixels as a polynomial function of the x and y position of each
   pixel. It fits the intensity to the polynomial, *A x* :sup:`2` *+ B
   y* :sup:`2` *+ C xy + D x + E y + F*. This will produce a smoothed
   image with a single peak or trough of intensity that tapers off
   elsewhere in the image. For many microscopy images (where the
   illumination of the lamp is brightest in the center of field of
   view), this method will produce an image with a bright central region
   and dimmer edges. But, in some cases the peak/trough of the
   polynomial may actually occur outside of the image itself.
-  *%(GAUSSIAN_FILTER)s:* This method convolves the image with a
   Gaussian whose full width at half maximum is the artifact diameter
   entered. Its effect is to blur and obscure features smaller than the
   specified diameter and spread bright or dim features larger than the
   specified diameter.
-  *%(MEDIAN_FILTER)s:* This method finds the median pixel value within
   the diameter you specify. It removes bright or dim features
   that are significantly smaller than the specified diameter.
-  *%(MEDIAN_FILTER_SCIPY)s:* This method finds the median pixel value within
   the diameter you specify. The scipy inplementation has been taken, as the centrosome
   one has problems with small diameter artifacts (=single pixel artefacts)
   Does NOT work with masked images!
-  *%(SMOOTH_KEEPING_EDGES)s:* This method uses a bilateral filter
   which limits Gaussian smoothing across an edge while applying
   smoothing perpendicular to an edge. The effect is to respect edges in
   an image while smoothing other features. *%(SMOOTH_KEEPING_EDGES)s*
   will filter an image with reasonable speed for artifact diameters
   greater than 10 and for intensity differences greater than 0.1. The
   algorithm will consume more memory and operate more slowly as you
   lower these numbers.
-  *%(CIRCULAR_AVERAGE_FILTER)s:* This method convolves the image with
   a uniform circular averaging filter whose size is the artifact
   diameter entered. This filter is useful for re-creating an
   out-of-focus blur to an image.
-  *%(SM_TO_AVERAGE)s:* Creates a flat, smooth image where every pixel
   of the image equals the average value of the original image.
-  *%(CLIP_HOT_PIXELS)s:* Clips hot pixels to the maximum local neighbor
   intensity. Hot pixels are identified by thresholding on the difference
   between the pixel intensity and it's maximum local neighbor intensity.

*Note, when deciding between %(MEDIAN_FILTER)s and %(GAUSSIAN_FILTER)s
we typically recommend
%(MEDIAN_FILTER)s over %(GAUSSIAN_FILTER)s because the
median is less sensitive to outliers, although the results are also
slightly less smooth and the fact that images are in the range of 0
to 1 means that outliers typically will not dominate too strongly
anyway.*
""" % globals())

        self.wants_automatic_object_size = cps.Binary(
                'Calculate artifact diameter automatically?', True, doc="""\
*(Used only if “%(GAUSSIAN_FILTER)s”, “%(MEDIAN_FILTER)s”, “%(SMOOTH_KEEPING_EDGES)s” or “%(CIRCULAR_AVERAGE_FILTER)s” is selected)*

Select *%(YES)s* to choose an artifact diameter based on the size of
the image. The minimum size it will choose is 30 pixels, otherwise the
size is 1/40 of the size of the image.

Select *%(NO)s* to manually enter an artifact diameter.
""" % globals())

        self.object_size = cps.Float(
                'Typical artifact diameter', 16.0, doc="""\
*(Used only if choosing the artifact diameter automatically is set to
“%(NO)s”)*

Enter the approximate diameter (in pixels) of the features to be blurred
by the smoothing algorithm. This value is used to calculate the size of
the spatial filter. %(HELP_ON_MEASURING_DISTANCES)s For most
smoothing methods, selecting a diameter over ~50 will take substantial
amounts of time to process.
""" % globals())

        self.sigma_range = cps.Float(
                'Edge intensity difference', 0.1, doc="""\
*(Used only if “%(SMOOTH_KEEPING_EDGES)s” is selected)*

Enter the intensity step (which indicates an edge in an image) that you
want to preserve. Edges are locations where the intensity changes
precipitously, so this setting is used to adjust the rough magnitude of
these changes. A lower number will preserve weaker edges. A higher
number will preserve only stronger edges. Values should be between zero
and one. %(HELP_ON_PIXEL_INTENSITIES)s
""" % globals())

        self.clip = cps.Binary(
                'Clip intensities to 0 and 1?', True, doc="""\
*(Used only if "%(FIT_POLYNOMIAL)s" is selected)*

The *%(FIT_POLYNOMIAL)s* method is the only smoothing option that can
yield an output image whose values are outside of the values of the
input image. This setting controls whether to limit the image
intensity to the 0 - 1 range used by CellProfiler.

Select *%(YES)s* to set all output image pixels less than zero to zero
and all pixels greater than one to one.

Select *%(NO)s* to allow values less than zero and greater than one in
the output image.
""" % globals())

        self.hp_filter_size = cps.Integer(
                'Neighborhood filter size', 3, minval=3, doc="""\
*(Used only if "%(CLIP_HOT_PIXELS)s" is selected)*)

Enter the size of the local neighborhood filter (recommended value: 3). This
value will be used to define the maximum distance around an individual pixel
within which other pixels will be considered as neighboring pixels. Note that
this value can be interpreted as diameter and thus has to be odd.
""" % globals())

        self.hp_threshold = cps.Float(
                'Hot pixel threshold', 50, minval=0, doc="""\
*(Used only if "%(CLIP_HOT_PIXELS)s" is selected)*)

Enter the absolute threshold on the difference between the individual pixel
intensity and it's maximum local neighbor intensity. This value will be used
to identify "hot pixels" whose intensity values will be clipped to the maximum
local neighbor intensity.
""" % globals())

        self.scale_hp_threshold = cps.Binary(
                'Scale hot pixel threshold to image scale?', True, doc="""\
*(Used only if "%(CLIP_HOT_PIXELS)s" is selected)*)

Specify whether the hot pixel threshold should be scaled to the image scale or
be taken as-is. In CellProfiler, the image values are usually rescaled to 0-1.
Thus, if the threshold should be selected based on the raw data values, it
needs to be scaled.

Example: If the image data type was uint16, the image is rescaled automatically
upon import. Thus all absolute values now have to be divided by 2^16. For
example, if one wants to set a threshold of 100 counts, a value of either
100 / 2^16 = 0.0015 (no scaling) or 100 (scaling) needs to be specified.
""" % globals())

    def settings(self):
        return [self.image_name, self.filtered_image_name,
                self.smoothing_method, self.wants_automatic_object_size,
                self.object_size, self.sigma_range, self.clip,
                self.hp_filter_size, self.hp_threshold,
                self.scale_hp_threshold]

    def upgrade_settings(self, setting_values, variable_revision_number,
                         module_name, from_matlab):
        if variable_revision_number < 2:
            setting_values += [3 , 20]  # hp_filter_size, hp_threshold
        if variable_revision_number < 4:
            setting_values.append(cps.NO)  # scale_hp_threshold
        return setting_values, variable_revision_number, from_matlab

    def visible_settings(self):
        result = [self.image_name, self.filtered_image_name,
                  self.smoothing_method]
        if self.smoothing_method.value not in [FIT_POLYNOMIAL, SM_TO_AVERAGE, CLIP_HOT_PIXELS]:
            result.append(self.wants_automatic_object_size)
            if not self.wants_automatic_object_size.value:
                result.append(self.object_size)
            if self.smoothing_method.value == SMOOTH_KEEPING_EDGES:
                result.append(self.sigma_range)
        if self.smoothing_method.value == FIT_POLYNOMIAL:
            result.append(self.clip)
        if self.smoothing_method.value == CLIP_HOT_PIXELS:
            result.append(self.hp_filter_size)
            result.append(self.hp_threshold)
            result.append(self.scale_hp_threshold)
        return result

    def run(self, workspace):
        image = workspace.image_set.get_image(self.image_name.value,
            must_be_grayscale=False)
        hp_threshold = self.hp_threshold.value
        if self.scale_hp_threshold.value is True:
            hp_threshold /= image.scale
        if len(image.pixel_data.shape) == 3:
            if self.smoothing_method.value == CLIP_HOT_PIXELS:
                # TODO support masks
                hp_filter_shape = (self.hp_filter_size.value, self.hp_filter_size.value, 1)
                output_pixels = SmoothMultichannel.clip_hot_pixels(image.pixel_data, hp_filter_shape, hp_threshold)
            else:
                output_pixels = image.pixel_data.copy()
                for channel in range(image.pixel_data.shape[2]):
                    output_pixels[:, :, channel] = self.run_grayscale(image.pixel_data[:, :, channel], image)
        else:
            if self.smoothing_method.value == CLIP_HOT_PIXELS:
                # TODO support masks
                hp_filter_shape = (self.hp_filter_size.value, self.hp_filter_size.value)
                output_pixels = SmoothMultichannel.clip_hot_pixels(image.pixel_data, hp_filter_shape, hp_threshold)
            else:
                output_pixels = self.run_grayscale(image.pixel_data, image)
        output_image = cpi.Image(output_pixels, parent_image=image)
        workspace.image_set.add(self.filtered_image_name.value, output_image)
        workspace.display_data.pixel_data = image.pixel_data
        workspace.display_data.output_pixels = output_pixels

    def run_grayscale(self, pixel_data, image):
        if self.wants_automatic_object_size.value:
            object_size = min(30, max(1, np.mean(pixel_data.shape) / 40))
        else:
            object_size = float(self.object_size.value)
        sigma = object_size / 2.35
        if self.smoothing_method.value == GAUSSIAN_FILTER:
            def fn(image):
                return scind.gaussian_filter(image, sigma,
                                             mode='constant', cval=0)

            output_pixels = smooth_with_function_and_mask(pixel_data, fn,
                                                          image.mask)
        elif self.smoothing_method.value == MEDIAN_FILTER:
            output_pixels = median_filter(pixel_data, image.mask,
                                          object_size / 2 + 1)
        elif self.smoothing_method.value == MEDIAN_FILTER_SCIPY:
            output_pixels = ndimage.median_filter(pixel_data, int(np.ceil(object_size / 2 + 1)))
        elif self.smoothing_method.value == SMOOTH_KEEPING_EDGES:
            sigma_range = float(self.sigma_range.value)

            output_pixels = skimage.restoration.denoise_bilateral(
                image=pixel_data,
                multichannel=image.multichannel,
                sigma_color=sigma_range,
                sigma_spatial=sigma
            )
        elif self.smoothing_method.value == FIT_POLYNOMIAL:
            output_pixels = fit_polynomial(pixel_data, image.mask,
                                           self.clip.value)
        elif self.smoothing_method.value == CIRCULAR_AVERAGE_FILTER:
            output_pixels = circular_average_filter(pixel_data, object_size / 2 + 1, image.mask)
        elif self.smoothing_method.value == SM_TO_AVERAGE:
            if image.has_mask:
                mean = np.mean(pixel_data[image.mask])
            else:
                mean = np.mean(pixel_data)
            output_pixels = np.ones(pixel_data.shape, pixel_data.dtype) * mean
        else:
            raise ValueError("Unsupported smoothing method: %s" %
                             self.smoothing_method.value)
        return output_pixels

    def display(self, workspace, figure):
        figure.set_subplots((2, 2))
        original = workspace.display_data.pixel_data
        if len(original.shape) == 3:
            original = np.sum(original / original.shape[2], axis=2)
        figure.subplot_imshow_grayscale(0, 0,
                                        original,
                                        "Original: %s" % self.image_name.value)
        filtered = workspace.display_data.output_pixels
        if len(filtered.shape) == 3:
            filtered = np.sum(filtered / filtered.shape[2], axis=2)
        figure.subplot_imshow_grayscale(1, 0,
                                        filtered,
                                        "Filtered: %s" % self.filtered_image_name.value,
                                        sharexy=figure.subplot(0, 0))
        difference = workspace.display_data.pixel_data - workspace.display_data.output_pixels
        if len(difference.shape) == 3:
            difference = np.sum(difference / difference.shape[2], axis=2)
        figure.subplot_imshow_grayscale(0, 1,
                                        difference,
                                        "Difference: original - filtered",
                                        sharexy=figure.subplot(0, 0))

    @staticmethod
    def clip_hot_pixels(img, hp_filter_shape, hp_threshold):
        if hp_filter_shape[0] % 2 != 1 or hp_filter_shape[1] % 2 != 1:
            raise ValueError("Invalid hot pixel filter shape: %s" % str(hp_filter_shape))
        hp_filter_footprint = np.ones(hp_filter_shape)
        hp_filter_footprint[int(hp_filter_shape[0] / 2), int(hp_filter_shape[1] / 2)] = 0
        max_img = scind.maximum_filter(img, footprint=hp_filter_footprint, mode='reflect')
        hp_mask = img - max_img > hp_threshold
        img = img.copy()
        img[hp_mask] = max_img[hp_mask]
        return img


# The image splitter for the streak finder algorithm
#
# Splits generated images into subimages for faster calculations
#
# Usage: "$ python image_splitter.py dataset x_parts y_parts"
# Author: Anatolii Zenin

import astropy.io.fits as pyfits
import numpy as np
import math
import os
import scipy.stats as stats
import shutil
import argparse

class streakImage:
    def __init__(self, width = 0, height = 0):
        self.im_w = width
        self.im_h = height
        self.hdr = 0
        self.noise_mu = 0
        self.noise_sigma = 1
        self.data = (np.zeros((width, height)))

    def get_wh(self):
        return (self.im_w, self.im_h)

    def openFile(self, fname):
        f1 = pyfits.open(fname)
        self.hdr = f1[0].header
        self.data = f1[0].data
        self.im_w = self.data.shape[0]
        self.im_h = self.data.shape[1]
        self.noise_mu = np.median(self.data)
        self.noise_sigma = stats.median_abs_deviation(self.data.flatten())
        print(self.noise_mu)
        print(self.noise_sigma)

    def savePart(self, fname, dirname, parts_x, parts_y):
        stepx = int(math.floor(self.data.shape[0] / parts_x))
        stepy = int(math.floor(self.data.shape[1] / parts_y))
        for i in range(parts_x):
            for j in range(parts_y):
                data = self.data[i * stepx : (i+1) * stepx,
                                        j * stepy : (j+1) * stepy]
                hdu = pyfits.PrimaryHDU(data)
                if(self.hdr == 0):
                    print('self.hdr == %d' % self.hdr)
                    hdu.header .set('NAXIS', 2)
                    hdu.header .set('NAXIS1', self.im_w)
                    hdu.header .set('NAXIS2', self.im_h)
                    hdu.header .set('BINX', 1)
                    hdu.header .set('BINY', 1)
                else:
                    hdu.header = self.hdr
                fname = fname.split('.')[0]
                partname = dirname + '/' + fname + '_' + 'x%d' % i + 'y%d' % j + '.fit'
                hdu.writeto(partname)
                f1 = pyfits.open(partname)
                hdr = f1[0].header
                f1.close()


def image_splitter(set_name, xparts, yparts):
    np.random.seed(0)
    setname = set_name.split('.')[0]
    dirname = "%s_part" % setname
    if(os.path.exists(dirname)):
        shutil.rmtree(dirname)
    os.mkdir(dirname)
    st = streakImage()
    f_name = "%s_bkg" % setname + '.fit'
    st.openFile(f_name)
    st.savePart(f_name, dirname, xparts, yparts)
    f_name = "%s_new" % setname + '.fit'
    st.openFile(f_name)
    st.savePart(f_name, dirname, xparts, yparts)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("setname")
    parser.add_argument("xparts")
    parser.add_argument("yparts")
    args = parser.parse_args()
    image_splitter(args.setname, int(args.xparts), int(args.yparts))

# The streakImage class for the streak finder algorithm
#
# Generates two images - without a streak and with a streak
# stars are distributed in the same way, noise is re-generated
#
# Usage: "$ python make_streak.py star_number streak_number star_snr streak_snr seed"
# Author: Anatolii Zenin

import astropy.io.fits as pyfits
import numpy as np
import math
import os
import scipy.stats as stats
from tqdm import tqdm
from functools import lru_cache
from scipy.ndimage import gaussian_filter
from skimage.draw import line
from contextlib import redirect_stdout
import argparse

class streakImage:
    def __init__(self, width = 0, height = 0, seed = 15, psf = 1.5):
        self.seed = seed
        self.initial_seed = seed
        self.addedNoise = 0
        self.psf = psf
        self.im_w = width
        self.im_h = height
        self.hdr = 0
        self.noise_mu = 100
        self.noise_sigma = 10
        self.data = (np.zeros((width, height)))
        self.streak_mask = np.zeros((self.im_w, self.im_h), dtype = np.short)
        self.streak_num = 0
        self.star_num = 0
        self.star_snr = 0
        self.streak_snr = 0

    def get_wh(self):
        return (self.im_w, self.im_h)

    def openFile(self, fname):
        f1 = pyfits.open(fname)
        self.hdr = f1[0].header
        self.data = f1[0].data
        self.im_w = self.data.shape[0]
        self.im_h = self.data.shape[1]
        self.noise_mu = np.median(self.data)
        self.noise_sigma = stats.median_abs_deviation(self.data.flatten()) * 1.4826
        print(self.noise_mu)
        print(self.noise_sigma)
        # print(repr(self.hdr))

    def getData(self):
        return self.data

    def readData(self, input):
        self.data = input

    def addNoise(self):
        rng = np.random.default_rng(self.seed)
        mean = 0
        if self.addedNoise == 0:
            mean = self.noise_mu
        noise = rng.normal(mean, self.noise_sigma, ((self.im_w, self.im_h)))
        self.data += noise
        self.seed += 1
        self.noise_mu = 0
        self.addedNoise = 1

    def addRandStars(self, number, max_snr):
        rng = np.random.default_rng(self.seed)
        self.star_snr = max_snr
        self.star_num = number
        x = rng.integers(0, self.im_w, (number))
        y = rng.integers(0, self.im_h, (number))
        star_img = np.zeros((self.im_w, self.im_h))
        for i in range(number):
            snr = rng.uniform(0, max_snr)
            star_img[x[i]][y[i]] = snr * self.noise_sigma * 2 * np.pi * self.psf**2
        star_img = gaussian_filter(star_img, sigma = self.psf)
        self.data += star_img

    def addStreak(self, a = 0.5, b = 0, snr = 5, start = 0, end = -1):
        if end == -1:
            end = self.im_h;
        line_img = np.zeros((self.im_w, self.im_h))
        if np.absolute(a) >= 1:
            for i in tqdm(range(start, end)):
                for j in range(math.ceil(np.absolute(a))): # to avoid line breaks when a > 1
                    if(int(a * i + b) + j > 0 and int(a * i + b) + j < self.im_w):
                        line_img[int(a * i + b) + j][i] += snr * self.noise_sigma

        if np.absolute(a) < 1:
            for i in tqdm(range(start, end)):
                if(int(a * i + b) > 0 and int(a * i + b) < self.im_w):
                    line_img[int(a * i + b)][i] += snr * self.noise_sigma

        gauss_line_img = gaussian_filter(line_img, sigma = self.psf)
        ratio = 1. / gauss_line_img[line_img>0].mean() * line_img[line_img>0].mean()
        print('ratio: %f' % ratio)
        self.streak_mask[gauss_line_img > 0.0001] = 1
        self.data += gauss_line_img * ratio

    def addStreakPoints(self, x0, y0, x1, y1, snr = 5):
        xl, yl = line(x0, y0, x1, y1)
        line_img = np.zeros((self.im_w, self.im_h))
        line_img[xl, yl] = snr * self.noise_sigma
        gauss_line_img = gaussian_filter(line_img, sigma = self.psf)
        ratio = 1. / gauss_line_img[line_img>0].max() * line_img[line_img>0].max()
        print('ratio: %f' % ratio)
        self.streak_mask[gauss_line_img > 0.0001] = 1
        self.data += gauss_line_img * ratio

    def saveFile(self, fname):
        if(os.path.isfile(fname)):
            os.remove(fname)
        hdu = pyfits.PrimaryHDU(self.data)
        if(self.hdr == 0):
            print('self.hdr == %d' % self.hdr)
            hdu.header.set('NAXIS', 2)
            hdu.header.set('NAXIS1', self.im_w)
            hdu.header.set('NAXIS2', self.im_h)
            hdu.header.set('BINX', 1)
            hdu.header.set('BINY', 1)
        else:
            hdu.header = self.hdr
        hdu.writeto(fname)
        f1 = pyfits.open(fname)
        hdr = f1[0].header
        print(repr(hdr))
        f1.close()
        with open('input_params.txt', 'w') as f:
            with redirect_stdout(f):
                print("# star_number\tstreak_number\tstar_snr\tstreak_snr\tseed")
                print("%d\t%d\t%.1f\t%.1f\t%d" % (self.star_num, self.streak_num, self.star_snr, self.streak_snr, self.initial_seed))

    def saveFileMask(self, fname):
        if(os.path.isfile(fname)):
            os.remove(fname)
        hdu = pyfits.PrimaryHDU(self.streak_mask)
        if(self.hdr == 0):
            print('self.hdr == %d' % self.hdr)
            hdu.header.set('NAXIS', 2)
            hdu.header.set('NAXIS', 2)
            hdu.header.set('NAXIS1', self.im_w)
            hdu.header.set('NAXIS2', self.im_h)
            hdu.header.set('BINX', 1)
            hdu.header.set('BINY', 1)
        else:
            hdu.header = self.hdr
        hdu.writeto(fname)
        f1 = pyfits.open(fname)
        hdr = f1[0].header
        print(repr(hdr))
        f1.close()

    def addRandStreaks(self, str_num, snr_low, snr_high):
        if(self.streak_num > 0):
            print("WARNING: adding streaks repeatedly")
        self.streak_num = str_num
        self.streak_snr = snr_high
        rng = np.random.default_rng(self.seed)
        for i in range(str_num):
            x0 = rng.integers(low = 0, high = self.im_w)
            y0 = rng.integers(low = 0, high = self.im_h)
            x1 = rng.integers(low = 0, high = self.im_w)
            y1 = rng.integers(low = 0, high = self.im_h)
            snr = rng.uniform(snr_low, snr_high)
            self.addStreakPoints(x0, y0, x1, y1, snr = snr)


# create a set of files with streaks
def makeRandStreakFile(star_number, streak_number, star_snr, streak_snr, seed = 15):
    image_w = 9576
    image_h = 6354
    fname = 'out'
    st = streakImage(image_w, image_h, seed = seed)
    st.addNoise()
    st.addRandStars(star_number, star_snr)
    st.saveFile("%s_bkg.fit" % fname)
    st.addNoise()
    st.addRandStreaks(streak_number, streak_snr, streak_snr)
    st.saveFile("%s_new.fit" % fname)
    st.saveFileMask("%s_mask.fit" % fname)

# create a file with random stars and streaks
def makeRandStreakFilesLAST(fname, seed = 15):
    image_w = 9576
    image_h = 6354
    rng = np.random.default_rng(seed)
    streak_snr = rng.uniform(0., 100)
    star_snr = rng.uniform(0., 1000)
    streak_number = rng.integers(1, 20)
    star_number = rng.integers(1, 500000)
    st = streakImage(image_w, image_h, seed = seed)
    st.addNoise()
    st.addRandStars(star_number, star_snr)
    st.saveFile("%s_bkg.fit" % fname)
    st.addNoise()
    st.addRandStreaks(streak_number, streak_snr, streak_snr)
    st.saveFile("%s_new.fit" % fname)
    st.saveFileMask("%s_mask.fit" % fname)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    # parser.add_argument("filename")
    parser.add_argument("star_num")
    parser.add_argument("streak_num")
    parser.add_argument("star_snr")
    parser.add_argument("streak_snr")
    parser.add_argument("seed")
    args = parser.parse_args()
    if args.star_num == '':
        makeRandStreakFilesLAST(args.seed,
                                seed = int(args.seed))
    else:
        makeRandStreakFile(int(args.star_num),
                           int(args.streak_num),
                           float(args.star_snr),
                           float(args.streak_snr),
                           seed = int(args.seed))

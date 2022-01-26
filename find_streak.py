# The StreakFinder class for the streak finder algorithm
#
# Finds streaks in images using probabilistic Hough and combines them
# into a mask. Compares the result with a true mask.
# Uses score maps from the proper subtraction paper:
# https://iopscience.iop.org/article/10.3847/0004-637X/830/1/27/pdf
#
# Usage: "$ python find_streak.py dataset x_parts y_parts"
# Author: Anatolii Zenin

import numpy as np
from astropy.io import fits
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
from astropy.io import fits
from skimage.transform import hough_line, hough_line_peaks
from skimage.transform import probabilistic_hough_line
from skimage.draw import line
from skimage import data
import glob
import os
import pandas as pd
from scipy import ndimage
from skimage.morphology import square, dilation
import math
import timeit
from contextlib import redirect_stdout
import argparse


class StreakFinder:
    def __init__(self, dset = '', x_parts = 0, y_parts = 0):
        self.psf = 1.5
        self.x_parts = x_parts # number of x-axis parts of split image
        self.y_parts = y_parts # number of x-axis parts of split image
        self.xmax = 0
        self.ymax = 0
        self.smap_xlen = 0 # x size of a score map
        self.smap_ylen = 0 # x size of a score map
        self.dset = dset # dataset name
        self.full_image = np.zeros(())
        self.true_mask = np.zeros(()) # true streak mask
        self.base_mask = np.zeros(()) # calculated streak mask
        self.df = pd.DataFrame()
        self.streak_counter = 0
        self.tp = 0
        self.fp = 1

        if self.dset != '':
            open_dataset(self.dset, x_parts = self.x_parts, y_parts= self.y_parts)

    # calculate the distance between two points
    def getdist(self, line):
        xd = line[0][0] - line[1][0]
        yd = line[0][1] - line[1][1]
        return np.sqrt(np.power(xd, 2) + np.power(yd, 2))

    # convert line (list) to coordinate tuple
    def dfToXY(self, df_streak):
        x0 = df_streak["x0"]
        x1 = df_streak["x1"]
        y0 = df_streak["y0"]
        y1 = df_streak["y1"]
        return x0, y0, x1, y1

    #get slope and intersect for dataframe
    def getAB(self, streak, transpose = 0):
        x0 = streak[0][1]
        y0 = streak[0][0]
        x1 = streak[1][1]
        y1 = streak[1][0]
        if(transpose):
            x0 = streak[0][0]
            y0 = streak[0][1]
            x1 = streak[1][0]
            y1 = streak[1][1]
        if x0 == x1: # if the line is vertical
            a = np.nan
            b = np.nan
            return a, b, x0, transpose
        else:
            a = (y1 - y0)/(x1 - x0)
            b = y0 - x0 * a
        return a, b, np.nan, transpose

    #get slope and intersect from a transposed image
    def getAB_t(self, streak):
        x0 = streak[0][1]
        y0 = streak[0][0]
        x1 = streak[1][1]
        y1 = streak[1][0]
        transpose = 0
        a = np.nan
        b = np.nan
        if x0 == x1: # if the line is vertical
            transpose = 1
        elif(abs((y1 - y0)/(x1 - x0)) > 1):
            transpose = 1
        if(transpose):
            x0 = streak[0][0]
            y0 = streak[0][1]
            x1 = streak[1][0]
            y1 = streak[1][1]
        a = (y1 - y0)/(x1 - x0)
        b = y0 - x0 * a
        return a, b, np.nan, transpose


    def open_dataset(self, dset, x_parts, y_parts):
        self.dset = dset
        lines_sample = []
        global_lines = []
        image_folder_smaps = self.dset + '_part/'
        full_image_path = self.dset + '_new.fit'
        true_mask_path = self.dset + '_mask.fit'
        full_image_file = fits.open(full_image_path)
        true_mask_file = fits.open(true_mask_path)
        self.full_image = full_image_file[0].data
        self.true_mask = true_mask_file[0].data
        full_image_file.close()
        true_mask_file.close()

        # create empty mask for full image
        mask = np.zeros((self.full_image.shape[0],self.full_image.shape[1]))
        self.xmax, self.ymax = mask.shape
        # properties of subimages and other internal stuff
        smap_file = fits.open(self.dset + '_part/' + 'x0y0_smap.fit')
        smap = smap_file[0].data
        self.smap_xlen, self.smap_ylen = smap.shape
        cutoff_smaps = np.zeros((x_parts, y_parts, self.smap_xlen, self.smap_ylen))
        smap_file.close()
        # apply a histogram cutoff to score maps
        line_array = np.zeros((x_parts, y_parts, 2))
        for i in range(x_parts):
            for j in range(y_parts):
                filename = self.dset + '_part/' + 'x%d' % i + 'y%d' % j + '_smap.fit'
                with fits.open(filename) as hdul:
                    image=hdul[0].data
                    image[np.isnan(image)] = 0
                    histogram = image.flatten()
                    try:
                        histogram=histogram[histogram < 0]
                        quantile = np.quantile(histogram, 0.32)
                    except:
                        quantile = 0
                    cutoff = abs(quantile)
                    image[image < cutoff] = 0
                    cutoff_smaps[i][j] = image

        for i in range(cutoff_smaps.shape[0]):
            for j in range(cutoff_smaps.shape[1]):
                lines = list(probabilistic_hough_line(cutoff_smaps[i][j],
                             threshold=0, line_length=int(20*self.psf), line_gap=0))
                for k in range(len(lines)):
                    lines[k] = list(lines[k])
                    lines[k][0] = list(lines[k][0])
                    lines[k][1] = list(lines[k][1])
                    lines[k][0][0] += j*smap.shape[1]
                    lines[k][1][0] += j*smap.shape[1]
                    lines[k][0][1] += i*smap.shape[0]
                    lines[k][1][1] += i*smap.shape[0]
                    global_lines.append(lines[k])
        print(len(global_lines))

        streaks = np.zeros((len(global_lines), 5))
        for i in range(streaks.shape[0]):
            streaks[i][0] = global_lines[i][0][0]
            streaks[i][1] = global_lines[i][0][1]
            streaks[i][2] = global_lines[i][1][0]
            streaks[i][3] = global_lines[i][1][1]
            streaks[i][4] = self.getdist(global_lines[i])
        df = pd.DataFrame(streaks, columns = ["x0", "y0", "x1", "y1", "length"])
        df = df.assign(a = np.nan)
        df = df.assign(b = np.nan)
        df = df.assign(x_vert = np.nan)
        df = df.assign(tp = np.nan)
        df = df.assign(group = -1)
        df = df.assign(group_len = 0)
        df = df.assign(streak_len = 0)
        df["x0"] = df["x0"].astype(np.int32)
        df["y0"] = df["y0"].astype(np.int32)
        df["x1"] = df["x1"].astype(np.int32)
        df["y1"] = df["y1"].astype(np.int32)
        df = df.sort_values("length", ascending = [False]).reset_index(drop=True)
        for i in range(len(df)):
            x0, y0, x1, y1 = self.dfToXY(df.loc[i])
            a, b, xx, tp = self.getAB_t([[x0, y0], [x1, y1]])
            df.at[i, "a"] = a
            df.at[i, "b"] = b
            df.at[i, "x_vert"] = xx
            df.at[i, "tp"] = tp
        self.df = df

    def drawLine(self, image, x0, y0, x1, y1):
        xl, yl = line(x0, y0, x1, y1)
        image[xl, yl] = 1

    def combine_streak(self, str_num = 0):
        newmask = np.zeros((self.full_image.shape[0], self.full_image.shape[1]), dtype=np.bool_)
        if(len(self.df.index[self.df['group'] == str_num].tolist())):
            print("streak number %d already exists" % str_num)
            return newmask, 0
        indices = self.df.index[self.df['group'] == -1].tolist()
        if(len(indices) == 0):
            print("no more lines to check")
            return newmask, -1
        print('\nlongest line: i = %d, len = %f' % (indices[0], self.df.at[indices[0], 'length']))

        chng_cntr = 0
        group_len = 0
        # calculate alignment using slope and intercept
        slope_tol_deg = 2
        slope_tol = slope_tol_deg * np.pi / 180
        intersect_tol = 200
        print('%d streaks to check' % len(indices))
        tochange = self.df[self.df['group'] == -1].index[np.isclose(np.arctan(self.df["a"].loc[indices]),
                                       np.arctan(self.df["a"].loc[indices[0]]),
                                       atol = slope_tol)
                            & np.isclose(self.df["b"].loc[indices],
                                         self.df["b"].loc[indices[0]],
                                         atol = intersect_tol)
                            & np.asarray(self.df["tp"].loc[indices] == 0)]
        if len(tochange) == 0:
            tochange = self.df[self.df['group'] == -1].index[np.isclose(np.arctan(self.df["a"].loc[indices]),
                                           np.arctan(self.df["a"].loc[indices[0]]),
                                           atol = slope_tol)
                                & np.isclose(self.df["b"].loc[indices],
                                             self.df["b"].loc[indices[0]],
                                             atol = intersect_tol)
                                & np.asarray(self.df["tp"].loc[indices] == 1)]

        self.df.at[tochange, "group"] = str_num
        self.df.at[tochange, "group_len"] = np.sum(self.df.loc[tochange]["length"])
        chng_cntr = len(tochange)
        print('got %d aligned streaks' % chng_cntr)


        # if no matches, exit
        if(chng_cntr == 0):
            self.df.at[indices[0], "group"]= -10
            return newmask, chng_cntr
        if(chng_cntr <= 10):
            print('too few streaks')
            indices = self.df.index[self.df['group'] == str_num].tolist()
            for i in range(len(indices)):
                self.df.at[indices[i], "group"] = -11
                self.df.at[indices[i], "group"] = -11
            return newmask, 0

        # make a new line
        sdf = self.df[self.df['group'] == str_num]
        min_y = np.min([sdf['y0'].min(), sdf['y1'].min()])
        max_y = np.max([sdf['y0'].max(), sdf['y1'].max()])
        if sdf['y0'].min() < sdf['y1'].min():
            min_x = np.mean(sdf[sdf["y0"] == min_y]['x0'])
            max_x = np.mean(sdf[sdf["y1"] == max_y]['x1'])
        else:
            min_x = np.mean(sdf[sdf["y1"] == min_y]['x1'])
            max_x = np.mean(sdf[sdf["y0"] == max_y]['x0'])
        self.drawLine(newmask, int(min_y), int(min_x), int(max_y), int(max_x))
        xd = max_x - min_x
        yd = max_y - min_y
        streak_len = np.sqrt(np.power(xd, 2) + np.power(yd, 2))

        self.df.at[tochange, "streak_len"] = streak_len

        actual_len = sdf["length"].sum()
        print('length ratio: %f' % (actual_len / streak_len))
        if actual_len / np.count_nonzero(newmask) < 0.5:
            print('inconsistent lengths')
            indices = self.df.index[self.df['group'] == str_num].tolist()
            for i in range(len(indices)):
                self.df.at[indices[i], "group"] = -12
            newmask = np.zeros((self.full_image.shape[0], self.full_image.shape[1]), dtype=np.bool_)
            return newmask, 0

        return newmask, chng_cntr

    def make_mask(self):
        self.base_mask = np.zeros((self.full_image.shape[0], self.full_image.shape[1]), dtype = np.bool_)
        while(True):
            newmask, chng_counter = self.combine_streak(self.streak_counter)
            if chng_counter == -1:
                break # stop if something went wrong
            elif chng_counter == 0:
                continue
            self.streak_counter += 1
            self.base_mask += newmask
        print(self.streak_counter)
        # make a streak wider
        struct2 = ndimage.generate_binary_structure(2, 2)
        self.base_mask = ndimage.binary_dilation(self.base_mask, structure=struct2,
                                                 iterations = 30).astype(self.base_mask.dtype)

    def get_mask(self):
        return self.base_mask, self.streak_counter

    def make_metrics(self):
        inv_true_mask = np.ones((self.full_image.shape[0], self.full_image.shape[1]))
        inv_true_mask -= self.true_mask
        # compare with the result
        targetpos = np.count_nonzero(self.true_mask > 0.1)
        targetneg = np.count_nonzero(inv_true_mask > 0.1)
        truepos = np.sum(self.base_mask[self.true_mask > 0.1])
        falsepos = np.sum(self.base_mask[self.true_mask < 0.1])
        self.tp = truepos / targetpos
        self.fp = falsepos / targetneg

    def save_result(self):
        with open('out_eff.txt', 'w') as f:
            with redirect_stdout(f):
                print("# true_positives\tfalse_positives")
                print("%.4f\t%.4f" % (self.tp, self.fp))

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("setname")
    parser.add_argument("xparts")
    parser.add_argument("yparts")
    args = parser.parse_args()

    sf = StreakFinder()
    sf.open_dataset(args.setname, int(args.xparts), int(args.yparts))
    sf.make_mask()
    sf.make_metrics()
    sf.save_result()
    # mask, cntr = sf.get_mask()
    # plt.imshow(mask)
    # plt.show()

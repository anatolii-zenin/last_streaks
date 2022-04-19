# LAST streak detector
A set of tools to make ROC curves of the streak detection algorithm

## Functionality
1. Generates images from scratch, randomly places stars and streaks. Both stars and streaks can have random SNR (signal-to-noise ratio). A true mask covering all the streaks in the generated image can also be produced.
2. Adds streaks to existing images and re-randomises background.
3. Searches for streaks in images and creates masks. Two input images are required - a backgound image without a streak and a new image with or without a streak. The naming convention is *dataset*_bkg.fit and *dataset*_new.fit.
4. Estimates the percentage of streak pixels covered by the calculated mask by comparing it with the true mask.
5. Generates sets of images with various star numbers, streak numbers and SNRs. Calculates the efficiency of the algorithm and plots a ROC curve.

## The algorithm
1. Take a pair of images - background image (w/o streaks) and new image (w/ streaks). Split them both into several (~20) parts.
2. Perform proper subtraction by E. Ofek (https://iopscience.iop.org/article/10.3847/0004-637X/830/1/27/pdf) on corresponding pairs of images. Obtain a set of score maps.
3. Clean up the score maps by removing 32% of the lowest values (cut everything below 1 sigma from the mean).
4. Run probabilistic Hough algorithm on the cleaned score maps, obtain a set of streak segments.
5. Check streak alignment and connect parts of the same streak. (Alignment checks can be improved)
6. Calculate masks for all streaks and combine them into a single mask. Compare with the true mask.
7. (Optional) Repeat 1-6 for multiple images with different parameters to get a ROC curve.

## How to use. Example.

## Description of the scripts, inputs and outputs
roc_launcher.zsh - the script for submitting the jobs to the DESY cluster. Runs roc_single.zsh with a set of given parameters.

roc_single.zsh - the script that is run on a cluster. Creates images in the scratch directory, analyzes them and puts the output where requested.
```
roc_single $1 $2 $3 $4 $5
```
where $1 - star number, $2 - streak number, $3 - star SNR, $4 streak SNR, $5 - seed.

make_streak.py - makes images with and without streaks
```
python make_streak.py $1 $2 $3 $4 $5
```
where $1 - star number, $2 - streak number, $3 - star SNR, $4 streak SNR, $5 - seed.

image_splitter.py - splits images for faster processing
```
python image_splitter.py dataset xparts yparts
```

find_streak.py - looks for streaks in images
```
python find_streak.py dataset xparts yparts
```
Here 'dataset' stands for a part of file names:
```
%s_bkg.fit % dataset
%s_new.fit % dataset
```
Dataset name is generated automatically if the script is ran from roc_launcher.zsh, otherwise must be specified by the user.

overall input - sets of star numbers, streak numbers, star signal-to-noise ratios (SNR), streak SNRs, seeds provided in roc_launcher.zsh

overall output - output parameters sorted by folders star number -> star_number/streak_number/star_SNR/streak_SNR/$seed_out.txt and $seed_in.txt



Output is analysed using streaks_stats.ipynb. Sample output is included.

In case of manual analysis two images are required

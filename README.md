# last_streaks
A set of tools to make ROC curves of the streak detection algorithm

roc_launcher.zsh - the script for submitting the jobs to the DESY cluster. Runs roc_single.zsh with a set of given parameters.
roc_single.zsh - the script that is run on a cluster. Creates images in the scratch directory, analyzes them and puts the output where requested.

make_streak.py - makes images with and without streaks

image_splitter.py - splits images for faster processing

find_streak.py - looks for streaks in images

overall input - sets of star numbers, streak numbers, star signal-to-noise ratios (SNR), streak SNRs, seeds provided in roc_launcher.zsh

overall output - output parameters sorted by folders star number -> star_number/streak_number/star_SNR/streak_SNR/$seed_out.txt and $seed_in.txt

Output is analysed using streaks_stats.ipynb. Sample output is included.

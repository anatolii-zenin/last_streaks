#!/bin/zsh
#
#(otherwise the default shell would be used)
#$ -S /bin/zsh
#(the running time for this job)
#$ -l h_rt=02:00:00
#
#(the maximum memory usage of this job)
#$ -l h_rss=6G
#$ -l tmpdir_size=3G
#
#(stderr and stdout are merged together to stdout)
#$ -j y
#
#(send mail on job's end and abort)
#$ -m ae

# change to scratch directory
export OUTDIR=/afs/ifh.de/user/a/azenin/roc_output
export HOME=/afs/ifh.de/user/a/azenin
source /afs/ifh.de/user/a/azenin/.zshrc
source /afs/ifh.de/user/a/azenin/.zprofile
conda activate last
cd $TMPDIR
# copy scripts to scratch directory
cp /lustre/fs23/group/icecube/azenin/last/roc_jobs/*.py $TMPDIR
cp /lustre/fs23/group/icecube/azenin/last/roc_jobs/makesmap.m $TMPDIR
echo 'copied'
python make_streak.py $1 $2 $3 $4 $5
python -- image_splitter.py out 5 4
echo 'python'
/opt/matlab/R2019a/bin/matlab -nodesktop -nosplash -r "makesmap out 5 4"
python find_streak.py out 5 4
# copy output where needed
mkdir $1
mkdir $1/$2
mkdir $1/$2/$3
mkdir $1/$2/$3/$4
mv std* $1/$2/$3/$4/
mv out_eff.txt $1/$2/$3/$4/"$5"_out.txt
mv input_params.txt $1/$2/$3/$4/"$5"_in.txt
cp -r $1 $OUTDIR

for i in 10000 100000 1000000 20000000 #stars
do
    for j in 1 5 10 #streaks
    do
        for k in 100 1000 10000 #star snr
        do
            for l in 0.5 1 10 100 1000 10000 #streak snr
            do
                qsub roc_single.zsh $i $j $k $l $1
            done
        done
        sleep 700 #give time to process
    done
done

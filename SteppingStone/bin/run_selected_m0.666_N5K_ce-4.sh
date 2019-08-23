for i in $(seq 0.0298 -0.0001 0.0001)
do
    echo $i
    mkdir $i
    python /home/booker/FishersWaveFst/simulations/bin/runSelectedRuns.py --direc /home/booker/FishersWaveFst/simulations/neutral_m0.666_N5K_formatted/ -s $i -m 0.666 -c 0.0001 -o $i/ --batch --reps 500 --nproc 45
done


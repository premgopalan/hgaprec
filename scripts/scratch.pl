
my @cmds2 = (
    "/scratch/pgopalan/gaprec/src/gaprec -dir $datadir/data/movielens-1m -n 6040 -m 3900 -k $K -rfreq 10 -label O",
    "/scratch/pgopalan/gaprec/src/gaprec -dir $datadir/data/movielens-1m -n 6040 -m 3900 -k $K -rfreq 10 -bias -label O",
    "/scratch/pgopalan/gaprec/src/gaprec -dir $datadir/data/movielens-1m -n 6040 -m 3900 -k $K -rfreq 10 -binary-data -label O",
    "/scratch/pgopalan/gaprec/src/gaprec -dir $datadir/data/movielens-1m -n 6040 -m 3900 -k $K -rfreq 10 -bias -binary-data -label O");


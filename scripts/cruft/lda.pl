#!/usr/bin/perl

my $datadir = "/scratch/pgopalan/hgaprec/analysis";
my $K=10;
my $SEED = 123456;
my @cmds = (
    #"/scratch/pgopalan/hgaprec/src/hgaprec -dir $datadir/data/movielens-1m -n 6040 -m 3900 -k $K -rfreq 10",
    "/scratch/pgopalan/hgaprec/src/hgaprec -dir $datadir/data/movielens-1m -n 6040 -m 3900 -k $K -rfreq 10 -bias",
    #"/scratch/pgopalan/hgaprec/src/hgaprec -dir $datadir/data/movielens-1m -n 6040 -m 3900 -k $K -rfreq 10 -hier",
    #"/scratch/pgopalan/hgaprec/src/hgaprec -dir $datadir/data/movielens-1m -n 6040 -m 3900 -k $K -rfreq 10 -binary-data",
    "/scratch/pgopalan/hgaprec/src/hgaprec -dir $datadir/data/movielens-1m -n 6040 -m 3900 -k $K -rfreq 10 -bias -binary-data",
    #"/scratch/pgopalan/hgaprec/src/hgaprec -dir $datadir/data/movielens-1m -n 6040 -m 3900 -k $K -rfreq 10 -hier -binary-data"
    );

my @echonest = (
    "/scratch/pgopalan/hgaprec/src/hgaprec -dir $datadir/data/echonest -n 1019318 -m 384546 -k $K -rfreq 1 -seed $SEED",
    "/scratch/pgopalan/hgaprec/src/hgaprec -dir $datadir/data/echonest -n 1019318 -m 384546 -k $K -rfreq 1 -bias -seed $SEED",
    "/scratch/pgopalan/hgaprec/src/hgaprec -dir $datadir/data/echonest -n 1019318 -m 384546 -k $K -rfreq 1 -hier -bias -seed $SEED",
    "/scratch/pgopalan/hgaprec/src/hgaprec -dir $datadir/data/echonest -n 1019318 -m 384546 -k $K -rfreq 1 -hier -seed $SEED",
    "/scratch/pgopalan/hgaprec/src/hgaprec -dir $datadir/data/echonest -n 1019318 -m 384546 -k $K -rfreq 1 -binary-data -seed $SEED",
    "/scratch/pgopalan/hgaprec/src/hgaprec -dir $datadir/data/echonest -n 1019318 -m 384546 -k $K -rfreq 1 -binary-data -hier -seed $SEED",
    "/scratch/pgopalan/hgaprec/src/hgaprec -dir $datadir/data/echonest -n 1019318 -m 384546 -k $K -rfreq 1 -binary-data -bias -seed $SEED",
    "/scratch/pgopalan/hgaprec/src/hgaprec -dir $datadir/data/echonest -n 1019318 -m 384546 -k $K -rfreq 1 -binary-data -bias -hier -seed $SEED",
    );

my @netflix = (
    "/scratch/pgopalan/hgaprec/src/hgaprec -dir $datadir/data/netflix -n 480189 -m 17770 -k $K -rfreq 1 -seed $SEED",
    "/scratch/pgopalan/hgaprec/src/hgaprec -dir $datadir/data/netflix -n 480189 -m 17770 -k $K -rfreq 1 -bias -seed $SEED",
    "/scratch/pgopalan/hgaprec/src/hgaprec -dir $datadir/data/netflix -n 480189 -m 17770 -k $K -rfreq 1 -hier -bias -seed $SEED",
    "/scratch/pgopalan/hgaprec/src/hgaprec -dir $datadir/data/netflix -n 480189 -m 17770 -k $K -rfreq 1 -hier -seed $SEED",
    "/scratch/pgopalan/hgaprec/src/hgaprec -dir $datadir/data/netflix -n 480189 -m 17770 -k $K -rfreq 1 -binary-data -seed $SEED",
    "/scratch/pgopalan/hgaprec/src/hgaprec -dir $datadir/data/netflix -n 480189 -m 17770 -k $K -rfreq 1 -binary-data -hier -seed $SEED",
    "/scratch/pgopalan/hgaprec/src/hgaprec -dir $datadir/data/netflix -n 480189 -m 17770 -k $K -rfreq 1 -binary-data -bias -seed $SEED",
    "/scratch/pgopalan/hgaprec/src/hgaprec -dir $datadir/data/netflix -n 480189 -m 17770 -k $K -rfreq 1 -binary-data -bias -hier -seed $SEED",
    );


my @netflixmsr1 = (
    "/scratch/pgopalan/hgaprec/src/hgaprec -dir $datadir/data/netflixmsr -n 478614 -m 17770 -k 20 -rfreq 10 -seed $SEED",
    "/scratch/pgopalan/hgaprec/src/hgaprec -dir $datadir/data/netflixmsr -n 478614 -m 17770 -k 20 -rfreq 10 -bias -seed $SEED",
    "/scratch/pgopalan/hgaprec/src/hgaprec -dir $datadir/data/netflixmsr -n 478614 -m 17770 -k 20 -rfreq 10 -hier -bias -seed $SEED",
    "/scratch/pgopalan/hgaprec/src/hgaprec -dir $datadir/data/netflixmsr -n 478614 -m 17770 -k 20 -rfreq 10 -hier -seed $SEED",
    "/scratch/pgopalan/hgaprec/src/hgaprec -dir $datadir/data/netflixmsr -n 478614 -m 17770 -k 20 -rfreq 10 -seed $SEED -binary-data",
    "/scratch/pgopalan/hgaprec/src/hgaprec -dir $datadir/data/netflixmsr -n 478614 -m 17770 -k 20 -rfreq 10 -bias -seed $SEED -binary-data",
    "/scratch/pgopalan/hgaprec/src/hgaprec -dir $datadir/data/netflixmsr -n 478614 -m 17770 -k 20 -rfreq 10 -hier -bias -seed $SEED -binary-data",
    "/scratch/pgopalan/hgaprec/src/hgaprec -dir $datadir/data/netflixmsr -n 478614 -m 17770 -k 20 -rfreq 10 -hier -seed $SEED -binary-data",
    "/scratch/pgopalan/hgaprec/src/hgaprec -dir $datadir/data/netflixmsr -n 478614 -m 17770 -k 100 -rfreq 10 -seed $SEED",
    "/scratch/pgopalan/hgaprec/src/hgaprec -dir $datadir/data/netflixmsr -n 478614 -m 17770 -k 100 -rfreq 10 -bias -seed $SEED",
    "/scratch/pgopalan/hgaprec/src/hgaprec -dir $datadir/data/netflixmsr -n 478614 -m 17770 -k 100 -rfreq 10 -hier -bias -seed $SEED",
    "/scratch/pgopalan/hgaprec/src/hgaprec -dir $datadir/data/netflixmsr -n 478614 -m 17770 -k 100 -rfreq 10 -hier -seed $SEED",
    );

my @netflixmsr2 = (
    "/scratch/pgopalan/hgaprec/src/hgaprec -dir $datadir/data/netflixmsr2 -n 478614 -m 17770 -k 20 -rfreq 10 -seed $SEED",
    "/scratch/pgopalan/hgaprec/src/hgaprec -dir $datadir/data/netflixmsr2 -n 478614 -m 17770 -k 20 -rfreq 10 -bias -seed $SEED",
    "/scratch/pgopalan/hgaprec/src/hgaprec -dir $datadir/data/netflixmsr2 -n 478614 -m 17770 -k 20 -rfreq 10 -hier -bias -seed $SEED",
    "/scratch/pgopalan/hgaprec/src/hgaprec -dir $datadir/data/netflixmsr2 -n 478614 -m 17770 -k 20 -rfreq 10 -hier -seed $SEED",
    "/scratch/pgopalan/hgaprec/src/hgaprec -dir $datadir/data/netflixmsr2 -n 478614 -m 17770 -k 100 -rfreq 10 -seed $SEED",
    "/scratch/pgopalan/hgaprec/src/hgaprec -dir $datadir/data/netflixmsr2 -n 478614 -m 17770 -k 100 -rfreq 10 -bias -seed $SEED",
    "/scratch/pgopalan/hgaprec/src/hgaprec -dir $datadir/data/netflixmsr2 -n 478614 -m 17770 -k 100 -rfreq 10 -hier -bias -seed $SEED",
    "/scratch/pgopalan/hgaprec/src/hgaprec -dir $datadir/data/netflixmsr2 -n 478614 -m 17770 -k 100 -rfreq 10 -hier -seed $SEED",
    );


my @mendeley = (
    "/scratch/pgopalan/hgaprec/src/hgaprec -dir $datadir/data/mendeley -n 80278 -m 261248 -k $K -rfreq 1 -seed $SEED",
    "/scratch/pgopalan/hgaprec/src/hgaprec -dir $datadir/data/mendeley -n 80278 -m 261248 -k $K -rfreq 1 -bias -seed $SEED",
    "/scratch/pgopalan/hgaprec/src/hgaprec -dir $datadir/data/mendeley -n 80278 -m 261248 -k $K -rfreq 1 -hier -bias -seed $SEED",
    "/scratch/pgopalan/hgaprec/src/hgaprec -dir $datadir/data/mendeley -n 80278 -m 261248 -k $K -rfreq 1 -hier -seed $SEED",
    "/scratch/pgopalan/hgaprec/src/hgaprec -dir $datadir/data/mendeley -n 80278 -m 261248 -k $K -rfreq 1 -binary-data -seed $SEED",
    "/scratch/pgopalan/hgaprec/src/hgaprec -dir $datadir/data/mendeley -n 80278 -m 261248 -k $K -rfreq 1 -binary-data -hier -seed $SEED",
    "/scratch/pgopalan/hgaprec/src/hgaprec -dir $datadir/data/mendeley -n 80278 -m 261248 -k $K -rfreq 1 -binary-data -bias -seed $SEED",
    "/scratch/pgopalan/hgaprec/src/hgaprec -dir $datadir/data/mendeley -n 80278 -m 261248 -k $K -rfreq 1 -binary-data -bias -hier -seed $SEED",
    );

my @movielens1m = (
    "/scratch/pgopalan/hgaprec/src/hgaprec -dir $datadir/data/movielens-1m -n 6040 -m 3900 -k $K -rfreq 1 -seed $SEED",
    "/scratch/pgopalan/hgaprec/src/hgaprec -dir $datadir/data/movielens-1m -n 6040 -m 3900 -k $K -rfreq 1 -bias -seed $SEED",
    "/scratch/pgopalan/hgaprec/src/hgaprec -dir $datadir/data/movielens-1m -n 6040 -m 3900 -k $K -rfreq 1 -hier -bias -seed $SEED",
    "/scratch/pgopalan/hgaprec/src/hgaprec -dir $datadir/data/movielens-1m -n 6040 -m 3900 -k $K -rfreq 1 -hier -seed $SEED",
    "/scratch/pgopalan/hgaprec/src/hgaprec -dir $datadir/data/movielens-1m -n 6040 -m 3900 -k $K -rfreq 1 -binary-data -seed $SEED",
    "/scratch/pgopalan/hgaprec/src/hgaprec -dir $datadir/data/movielens-1m -n 6040 -m 3900 -k $K -rfreq 1 -binary-data -hier -seed $SEED",
    "/scratch/pgopalan/hgaprec/src/hgaprec -dir $datadir/data/movielens-1m -n 6040 -m 3900 -k $K -rfreq 1 -binary-data -bias -seed $SEED",
    "/scratch/pgopalan/hgaprec/src/hgaprec -dir $datadir/data/movielens-1m -n 6040 -m 3900 -k $K -rfreq 1 -binary-data -bias -hier -seed $SEED",
    );

# SMALL

my @mendeleysmall = (
    "/scratch/pgopalan/hgaprec/src/hgaprec -dir $datadir/data/mendeley -n 10000 -m 40000 -k $K -rfreq 1 -seed $SEED",
    "/scratch/pgopalan/hgaprec/src/hgaprec -dir $datadir/data/mendeley -n 10000 -m 40000 -k $K -rfreq 1 -bias -seed $SEED",
    "/scratch/pgopalan/hgaprec/src/hgaprec -dir $datadir/data/mendeley -n 10000 -m 40000 -k $K -rfreq 1 -hier -bias -seed $SEED",
    "/scratch/pgopalan/hgaprec/src/hgaprec -dir $datadir/data/mendeley -n 10000 -m 40000 -k $K -rfreq 1 -hier -seed $SEED",
    "/scratch/pgopalan/hgaprec/src/hgaprec -dir $datadir/data/mendeley -n 10000 -m 40000 -k $K -rfreq 1 -binary-data -seed $SEED",
    "/scratch/pgopalan/hgaprec/src/hgaprec -dir $datadir/data/mendeley -n 10000 -m 40000 -k $K -rfreq 1 -binary-data -hier -seed $SEED",
    "/scratch/pgopalan/hgaprec/src/hgaprec -dir $datadir/data/mendeley -n 10000 -m 40000 -k $K -rfreq 1 -binary-data -bias -seed $SEED",
    "/scratch/pgopalan/hgaprec/src/hgaprec -dir $datadir/data/mendeley -n 10000 -m 40000 -k $K -rfreq 1 -binary-data -bias -hier -seed $SEED",
    );

foreach my $cmd (@netflixmsr2) {
    system("$cmd 2>&1 > /dev/null &");
}


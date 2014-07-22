#!/usr/bin/perl

my $datadir = "/scratch/pgopalan/hgaprec/analysis";
my $K=100;
my $N=69878;
my $M=10677;
my @cmds = (
    "/scratch/pgopalan/hgaprec/src/hgaprec -dir $datadir/data/movielens-10m -n $N -m $M -k $K -rfreq 1",
    "/scratch/pgopalan/hgaprec/src/hgaprec -dir $datadir/data/movielens-10m -n $N -m $M -k $K -rfreq 1 -bias",
    "/scratch/pgopalan/hgaprec/src/hgaprec -dir $datadir/data/movielens-10m -n $N -m $M -k $K -rfreq 1 -hier -bias",
    "/scratch/pgopalan/hgaprec/src/hgaprec -dir $datadir/data/movielens-10m -n $N -m $M -k $K -rfreq 1 -hier",
    "/scratch/pgopalan/hgaprec/src/hgaprec -dir $datadir/data/movielens-10m -n $N -m $M -k $K -rfreq 1 -binary-data",
    "/scratch/pgopalan/hgaprec/src/hgaprec -dir $datadir/data/movielens-10m -n $N -m $M -k $K -rfreq 1 -bias -binary-data",
    "/scratch/pgopalan/hgaprec/src/hgaprec -dir $datadir/data/movielens-10m -n $N -m $M -k $K -rfreq 1 -hier -binary-data",
    "/scratch/pgopalan/hgaprec/src/hgaprec -dir $datadir/data/movielens-10m -n $N -m $M -k $K -rfreq 1 -hier -bias -binary-data"
    );

my @cmds2 = (
    "/scratch/pgopalan/gaprec/src/gaprec -dir $datadir/data/movielens-10m -n $N -m $M -k $K -rfreq 10 -label O",
    "/scratch/pgopalan/gaprec/src/gaprec -dir $datadir/data/movielens-10m -n $N -m $M -k $K -rfreq 10 -bias -label O",
    "/scratch/pgopalan/gaprec/src/gaprec -dir $datadir/data/movielens-10m -n $N -m $M -k $K -rfreq 10 -binary-data -label O",
    "/scratch/pgopalan/gaprec/src/gaprec -dir $datadir/data/movielens-10m -n $N -m $M -k $K -rfreq 10 -bias -binary-data -label O");

foreach my $cmd (@cmds) {
    system("$cmd 2>&1 > /dev/null &");
}


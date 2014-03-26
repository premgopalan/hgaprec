#!/usr/bin/perl

use strict;
use warnings;


my $f = $ARGV[0];
print "using file $f\n";
open F, "<$f";

my $sum = 0;
while (<F>) {
    if ($_ = /(\d+)\s+(\d+)\s+(\d+).*/) {
	my $user = $1 + 0;
	my $item = $2 + 0;
	my $c = $3 + 0;
	$sum += $c;
    }
}

close F;

print "sum = $sum\n";
	
    

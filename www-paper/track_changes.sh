#!/bin/bash
#
# script to show word-style "track changes" from a previous revision
#

if [ $# -lt 1 ]
then
    echo 'please supply a revision hash argument (see "git log")'
    exit 1
fi

# get the revision to compare to
rev=$1
[ -d $rev ] || mkdir $rev

# pull old versions of .tex files
sections=`grep -i input bpf-www.tex | awk -F'[{}]' '{print $2}'`
for section in bpf-www.tex $sections
do
    git show $rev:www-paper/$section > $rev/$section
done

cd $rev

# symlink supporting files
for f in sig-alternate.cls figures newfigs bib.bib
do
    ln -s ../$f .
done

# generate diff with latexdiff, clean up, and produce pdf
latexdiff --flatten bpf-www.tex ../bpf-www.tex > changes.tex
dos2unix changes.tex
grep -v '^$' changes.tex > tmp
mv tmp changes.tex
pdflatex changes.tex

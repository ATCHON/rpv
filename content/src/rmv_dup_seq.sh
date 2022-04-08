#!/bin/bash

echo Concatenate sequences ...
cat *.fsa > $1.fsa
nb_sequence=$(cat $1.fsa | grep '>' | wc -l)
echo $1.fsa contains $nb_sequence sequences 
echo Remove potential duplicate sequences...
awk '/^>/{f=!d[$1];d[$1]=1}f' $1.fsa > $2.fsa
nb_sequence2=$(cat $2.fsa | grep '>' | wc -l)
echo $2.fsa contains $nb_sequence2 sequences
rm $1.fsa 
echo Finish
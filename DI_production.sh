#!/usr/bin/env bash


awk '{print $1,$2,$4}' $1 > $2.DI
awk '{if($2-$1>4)print}' $2.DI > $2_nonlocal.DI
sort -g -k 3 -r $2_nonlocal.DI > $2_ranked.DI
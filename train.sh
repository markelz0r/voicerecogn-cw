#!/bin/bash
declare -a names=( "Adrian" "Ali" "Andrew" "Andy" "Ce" "Chaorong" "Jeremy" "Ke" "Liam" "Mateusz" "Minghong" "Martino" "Nicole" "Nicholas" "Oliver" "Sarah" "Shaun" "Travis" "Vincent" "Vinny" "Silent")

for i in "${names[@]}"
do
	HInit -S lists/list.txt -l $i -L labels/train -M hmms -o $i -T 1 lib/proto10States.txt;
done

#!/bin/bash
combineCards.py -S central*incb.txt > ${1}.txt
combine -M Asymptotic ${1}.txt

#!/bin/bash
combineCards.py -S *incb.txt > ${1}.txt
combine -M Asymptotic ${1}.txt

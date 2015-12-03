#!/bin/bash
combineCards.py -S XXINPUTXX > ${1}.txt
combine -M Asymptotic ${1}.txt

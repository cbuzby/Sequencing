#!/usr/bin/env bash

for i in I II III IV V VI VII VIII IX X XI XII XIII XIV XV XVI; do python mp_inference.py -n 1000 C_A_$i.txt U_A_$i.txt -m contrast -np -o Aout_$i.txt; done


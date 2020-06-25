#!/bin/bash

estados=(
	SC
    MG
    AL
    SP
)

for INST in "${estados[@]}"; do
	for BETA in 10 20 30; do
		./PCTFGenetico "$INST" 30 "$BETA"
	done
done
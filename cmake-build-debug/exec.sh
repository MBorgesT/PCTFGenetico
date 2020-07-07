#!/bin/bash

estados=(
    MG
    SP
)

for INST in "${estados[@]}"; do
	./PCTFGenetico "$INST" 30 30
done
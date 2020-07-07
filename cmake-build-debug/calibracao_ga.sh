#!/bin/bash


for PRC_REM in 5 20; do
	for QTD_REM_GUL in 2 4; do
		for PRC_MUT in 10 80; do
			./PCTFGenetico SP 30 30 3 150.0 100 "$PRC_REM" "$QTD_REM_GUL" "$PRC_MUT" 10
		done
	done
done
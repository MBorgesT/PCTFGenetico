#!/bin/bash

for TAM_POP in 100 300; do
	for PRC_REM in 30 65; do
		for QTD_REM_GUL in 5 20; do
			for PRC_MUT in 25 80; do
				for PRC_TAM_CRO in 8 25; do
					./PCTFGenetico SP 30 30 5 240.0 "$TAM_POP" "$PRC_REM" "$QTD_REM_GUL" "$PRC_MUT" "$PRC_TAM_CRO"
				done
			done
		done
	done
done
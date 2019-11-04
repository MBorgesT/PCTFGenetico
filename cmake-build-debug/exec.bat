REM [executavel	inst	alfa	beta	num_exe	tem_ini	tem_con	tax_res	num_ite	tempo tam_pop tam_cem prc_mut]

FOR /L %%A IN (400,100,600) DO (
	FOR /L %%B IN (5,10,25) DO (
		FOR /L %%C IN (3,4,11) DO (
			CALL PCTFGenetico.exe	AC	30	10	5	5	%%A	%%B	%%C
		)
	)
)


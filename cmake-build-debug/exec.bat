REM [executavel	inst	alfa	beta	num_exe	tem_ini	tem_con	tax_res	num_ite	tempo tam_pop tam_cem prc_mut]

FOR /L %%A IN (100,300,700) DO (
	FOR /L %%B IN (5,15,35) DO (
		FOR /L %%C IN (5,10,25) DO (
			FOR /L %%D IN (10,45,100) DO (
				FOR /L %%E IN (10,45,100) DO (
					CALL PCTFGenetico.exe	AC	30	10	5	4	%%A	%%B	%%C	%%D	%%E
				)
			)
		)
	)
)


REM [executavel	inst	alfa	beta	num_exe	tem_ini	tem_con	tax_res	num_ite	tempo]

FOR /L %%A IN (60,40,300) DO (
	FOR /L %%B IN (20,10,60) DO (
		FOR /L %%C IN (15,5,45) DO (
			FOR /L %%D IN (5,5,30) DO (
				CALL PCTFGenetico.exe	AC	30	10	3	1	0.01	0.975	1	20	%%A	%%B	%%C	%%D
			)
		)
	)
)

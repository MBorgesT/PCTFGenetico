arq = open('calibracao_ga.txt', 'r')
linhas = arq.readlines()
arq.close()

resultados = []

for l in linhas:
	aux = l.split('\t')
	aux3 = []
	for i in aux:
		aux2 = i.split()
		aux3.append(aux2)

	resultados.append(aux3)

resultados.sort(key=lambda x:x[0][2], reverse=True)
for r in resultados:
	print(str(r) + '\n')
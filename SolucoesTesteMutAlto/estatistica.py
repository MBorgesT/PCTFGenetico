import statistics as st

arq = open('resumo.txt', 'r')
linhas = arq.readlines()
for i in range(len(linhas)):
    linhas[i] = linhas[i].split()
   
   
def aux(l):
    return(l[13])

linhas.sort(reverse=True, key=aux)
for i in range(len(linhas)):
    print(linhas[i])
    
aux1 = []
aux2 = []
for i in linhas:
    aux1.append(float(i[13]))
    aux2.append(float(i[14]))
print('media fo:', st.mean(aux1))
print('media tempo:', st.mean(aux2))
    

    
    
#include "pctfga.h"

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <omp.h>

#define DBG // habilita o modo DEBUG (exibe na tela as melhoras na FO)
#define RND // habilita m?ltiplas (NUM_EXE) execu??es com sementes aleat?rias para a inst?ncia
//#define DBG_TEMPO_GERACAO // habilita medidor de tempo de cada parte da gera??o da popula??o
//#define DBG_TEMPO_GA // habilita medidor de tempo de cada parte do algoritmo genetico
//#define DBG_TEMPO_SA
#define ESCREVER_SOL

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

//============================== DADOS DE ENTRADA ==============================
char INST[50] = "MG";  // arquivo com a inst?ncia (estado)
int ALFA = 30;    // limite de contadores instalados (% do total de arestas)
int BETA = 30;    // limite de faixas cobertas (% do total de faixas)
int NUM_EXE = 3;     // n?mero de execu??es do metodo
double MAX_TIME = 1200.0;  // tempo m?ximo de execu??o (segundos)
int MAX_ITE = 1;     // n?mero m?ximo de itera??es (MAX_ITE = MAX_ITE * n?mero de arestas)
// ------------------------------- GENETICO ---------------------------------
int TAM_POP = 100;
int PRC_REM = 65; // percentual de remo??o de contadores (usado na gera??o da popula??o inicial)
int QTD_REM_GUL = 3; // numero dentre quantas arestas serao testadas para decidir qual a melhor inserir na geracao de populacao
int PRC_MUT = 80; // chance de se ocorrer uma mutacao (geracao de vizinho)
int PRC_TAM_CRO = 10;
// ---------------------------------- SA ------------------------------------
int PRM_SAN = 20;
float TMP_INI = 3000;
float TMP_CNG = 0.01;
float TAX_RSF = 0.985;
int SAN_MAX = 10;
// ----------------------------- aux calibracao -----------------------------
int vetFOs[5];
//================================== PRINCIPAL =================================
int main(int argc, char *argv[]) {
    char arq[50];
#ifdef RND
    srand(time(NULL));
#else
    NUM_EXE = 1;
#endif
    if (argc > 1) {
        strcpy(INST, argv[1]);
        ALFA = atoi(argv[2]);
        BETA = atoi(argv[3]);
        NUM_EXE = atoi(argv[4]);
        MAX_TIME = atof(argv[5]);
        TAM_POP = atoi(argv[6]);
        PRC_REM = atoi(argv[7]);
        QTD_REM_GUL = atoi(argv[8]);
        PRC_MUT = atoi(argv[9]);
        PRC_TAM_CRO = atoi(argv[10]);
    }

    numCro_ = TAM_POP * ((float)PRC_TAM_CRO / 100);

    sprintf(arq, "../Instancias/%s.txt", INST);
    lerInstancia(arq);

    montarRede();


    #ifdef ESCREVER_SOL
    FILE *f = fopen("calibracao_ga.txt", "a");
    float mediaFOs = 0;
    int melhorFO = -1;

	printf("INST: %s\tALFA: %i\tBETA: %i\tNUM_EXE: %i\tMAX_TIME: %.1f\nTAM_POP: %i\tPRC_REM: %i\tQTD_REM_GUL: %i\tPRC_MUT: %i\tPRC_TAM_CRO: %i\n",
			INST, ALFA, BETA, NUM_EXE, MAX_TIME, TAM_POP, PRC_REM, QTD_REM_GUL, PRC_MUT, PRC_TAM_CRO);
	#endif

    for (int r = 0; r < NUM_EXE; r++) {
        execGA(r);
		#ifdef ESCREVER_SOL
        mediaFOs += vetFOs[r];
        if (vetFOs[r] > melhorFO)
        	melhorFO = vetFOs[r];
		#endif
    }
    mediaFOs /= NUM_EXE;
    float desvio = ((melhorFO - mediaFOs)/melhorFO) * 100;

	#ifdef ESCREVER_SOL
    fprintf(f, "MEDIA FOS: %.2f\tDESVIO: %f\tNUM_EXE: %i\tMAX_TIME: %.1f\tTAM_POP: %i\tPRC_REM: %i\tQTD_REM_GUL: %i\tPRC_MUT: %i\tPRC_TAM_CRO: %i\n",
			mediaFOs, desvio, NUM_EXE, MAX_TIME, TAM_POP, PRC_REM, QTD_REM_GUL, PRC_MUT, PRC_TAM_CRO);

    fclose(f);

    printf("MEDIA FOS: %.2f\tDESVIO: %f\n\n", mediaFOs, desvio);
	#endif

    return 0;
}
//==============================================================================
//==================================== GA ======================================

void execGA(const int r) {
	#ifdef DBG
	/*
	printf("INST:%s\tALFA:%i\tBETA:%i\tNUM_EXE:%i\tMAX_TIME:%.1f\tTAM_POP:%i\tPER_REM:%i\tQTD_REM_GUL:%i\tPRC_MUT:%i\n",
		   INST, ALFA, BETA, NUM_EXE, MAX_TIME, TAM_POP, PRC_REM, QTD_REM_GUL, PRC_MUT);
		   */
	printf("\nr: %i\n", r);
	#endif

    double h1, h2, h3;
    double tempoGastoEpi = 0, auxEpi1, auxEpi2;

    #ifdef DBG_TEMPO_GA
    double aux1, aux2, aux3;
    #endif

    int flag, qtdGen = 0, qtdEpi = 0;
    limRemocao_ = floor((PRC_REM / 100.0) * (numAre_ - maxContReal_));
    h1 = omp_get_wtime();
    gerarPopulacao();

    #ifdef DBG
    h3 = omp_get_wtime();
    printf("tempo geracao: %.4f\n", h3 - h1);
    #endif

    ordenarPopulacao(1);

    #ifdef DBG
	printf("NOVO: bstSol: %i\n", populacao[0].numParCob);
	#endif

    h2 = omp_get_wtime();
    while (h2 - h1 < MAX_TIME) {
        int p1, p2;

        #ifdef DBG_TEMPO_GA
        aux1 = omp_get_wtime();
        #endif

        for (int i = 0; i < numCro_; i++) {
            p1 = rand() % TAM_POP;
            do
                p2 = rand() % TAM_POP;
            while (p2 == p1);
            gerarFilho(p1, p2, TAM_POP + i);
            if (populacao[TAM_POP + i].numParCob > populacao[0].numParCob) {
                bstTime_ = omp_get_wtime() - h1;
                #ifdef DBG
                printf("NOVO: bstSol: %i\tbstTempo: %.4f\n", populacao[TAM_POP + i].numParCob, bstTime_);
                #endif
            }
        }

        #ifdef DBG_TEMPO_GA
        aux2 = omp_get_wtime();
        #endif

        ordenarPopulacao(TAM_POP);

        #ifdef DBG_TEMPO_GA
        aux3 = omp_get_wtime();
        #endif

        flag = 1;
        for (int i = 1; i < TAM_POP; i++) {
            if (populacao[i].numParCob != populacao[0].numParCob) {
                flag = 0;
                break;
            }
        }
        if (flag == 1) {
            #ifdef DBG
            auxEpi1 = omp_get_wtime();
            printf("epidemia\n");
            #endif
            epidemia();
            #ifdef DBG
            auxEpi2 = omp_get_wtime();
            tempoGastoEpi += auxEpi2 - auxEpi1;
            qtdEpi++;
            #endif
        }

        qtdGen++;

        h2 = omp_get_wtime();

        #ifdef DBG_TEMPO_GA
        if (qtdGen == 1) printf("geracao filho: %.4f\tordenacao: %.4f\tresto: %.4f\n", ((aux2-aux1)/(h2-aux1)), ((aux3-aux2)/(h2-aux1)), ((h2-aux3)/(h2-aux1)));
        #endif
    }



	#ifdef ESCREVER_SOL
	vetFOs[r] = populacao[0].numParCob;
    char path[100];
    sprintf(path, "../solucoes/%s/ALFA:%i-BETA:%i.txt", INST, ALFA, BETA);
    escreverResultado(populacao[0], path);
	#endif

    #ifdef DBG
    printf("FINAL: bstSol: %i\tbstTempo: %.4f\tqtd gen: %i\tqtd epi: %i\ttempo gasto em epi: %.3f\n\n\n",
           populacao[0].numParCob, bstTime_, qtdGen, qtdEpi, tempoGastoEpi);
    #endif
}

void gerarFilho(const int &p1, const int &p2, int f) {
	int aux = 1 + rand() % (numAre_ - 2);
	memset(populacao[f].vetAre, 0, sizeof(populacao[f].vetAre));
	populacao[f].numConIns = populacao[f].numFaiCob = populacao[f].numAreCom = 0;
	for (int i = 0; i < aux; i++) {
		populacao[f].vetAre[i] = populacao[p1].vetAre[i];
		populacao[f].numConIns += (populacao[p1].vetAre[i] * vetArestas_[i].nCon);
		populacao[f].numFaiCob += (populacao[p1].vetAre[i] * vetArestas_[i].nFai);
		populacao[f].numAreCom += (populacao[p1].vetAre[i]);
	}
	for (int i = aux; i < numAre_; i++) {
		populacao[f].vetAre[i] = populacao[p2].vetAre[i];
		populacao[f].numConIns += (populacao[p2].vetAre[i] * vetArestas_[i].nCon);
		populacao[f].numFaiCob += (populacao[p2].vetAre[i] * vetArestas_[i].nFai);
		populacao[f].numAreCom += (populacao[p2].vetAre[i]);
	}

	/*
	int pos;
	if ((populacao[f].numConIns > maxContReal_) || (populacao[f].numFaiCob > maxFaix_)) {
		while ((populacao[f].numConIns > maxContReal_) || (populacao[f].numFaiCob > maxFaix_)) {
			do
				pos = rand() % numAre_;
			while (populacao[f].vetAre[pos] == 0);
			populacao[f].vetAre[pos] = 0;
			populacao[f].numConIns -= vetArestas_[pos].nCon;
			populacao[f].numFaiCob -= vetArestas_[pos].nFai;
			populacao[f].numAreCom--;
		}
	} else {
		while ((populacao[f].numConIns <= maxContReal_) && (populacao[f].numFaiCob <= maxFaix_)) {
			do
				pos = rand() % numAre_;
			while (populacao[f].vetAre[pos] == 1);
			populacao[f].vetAre[pos] = 1;
			populacao[f].numConIns += vetArestas_[pos].nCon;
			populacao[f].numFaiCob += vetArestas_[pos].nFai;
			populacao[f].numAreCom++;
		}
		populacao[f].vetAre[pos] = 0;
		populacao[f].numConIns -= vetArestas_[pos].nCon;
		populacao[f].numFaiCob -= vetArestas_[pos].nFai;
		populacao[f].numAreCom--;
	}

    if (rand() % 100 < PRC_MUT){
    	gerVizinho(populacao[f]);
    }

    calcParCob(populacao[f]);

    */

	if (rand() % 100 < PRC_MUT){
		gerVizinho(populacao[f]);
	}

	calcParCob(populacao[f]);

	gulosidade(populacao[f]);

    if (rand() % 1000 < PRM_SAN) {
        execSA(f);
    }

}

void gerarPopulacao() {
    int foAntes, bstFO, bstPos, aux, limite, auxRandGul;

#ifdef DBG_TEMPO_GERACAO
    double aux1, aux2, aux3, aux4, aux5;
#endif

    for (int p = 0; p < TAM_POP + numCro_; p++) {
        memset(populacao[p].vetAre, 0, sizeof(populacao[p].vetAre));
        populacao[p].numConIns = populacao[p].numFaiCob = 0;
        populacao[p].numAreCom = numAre_;

#ifdef DBG_TEMPO_GERACAO
        aux1 = omp_get_wtime();
#endif
        // insere contador em todas as arestas
        for (int i = 0; i < numAre_; i++) {
            populacao[p].vetAre[i] = 1;
            populacao[p].numConIns += vetArestas_[i].nCon;
            populacao[p].numFaiCob += vetArestas_[i].nFai;
        }
        populacao[p].numParCob = maxPar_;

#ifdef DBG_TEMPO_GERACAO
        aux2 = omp_get_wtime();
#endif
        // remove alguns (at? o limite) contadores aleatoriamente
        limite = rand() % limRemocao_;
        for (int i = 0; i < limite; i++) {
            do
                aux = rand() % numAre_;
            while (populacao[p].vetAre[aux] == 0);
            populacao[p].vetAre[aux] = 0;
            populacao[p].numConIns -= vetArestas_[aux].nCon;
            populacao[p].numFaiCob -= vetArestas_[aux].nFai;
            populacao[p].numAreCom--;
        }

#ifdef DBG_TEMPO_GERACAO
        aux3 = omp_get_wtime();
#endif
        // remove contadores de forma GULOSA, at? encontrar uma solu??o vi?vel
        while ((populacao[p].numConIns > maxContReal_) || (populacao[p].numFaiCob > maxFaix_)) {
            bstFO = bstPos = -1;
            for (int i = 0; i < QTD_REM_GUL; i++) {
                do
                    auxRandGul = rand() % numAre_;
                while (populacao[p].vetAre[auxRandGul] == 0);
                foAntes = populacao[p].numParCob;
                populacao[p].vetAre[auxRandGul] = 0;
                calcParCob(populacao[p]);
                if (populacao[p].numParCob > bstFO) {
                    bstFO = populacao[p].numParCob;
                    bstPos = auxRandGul;
                }
                populacao[p].vetAre[auxRandGul] = 1;
                populacao[p].numParCob = foAntes;
            }
            populacao[p].vetAre[bstPos] = 0;
            populacao[p].numConIns -= vetArestas_[bstPos].nCon;
            populacao[p].numFaiCob -= vetArestas_[bstPos].nFai;
            populacao[p].numAreCom--;
            populacao[p].numParCob = bstFO;
        }

#ifdef DBG_TEMPO_GERACAO
        aux4 = omp_get_wtime();
#endif

        calcParCob(populacao[p]);

#ifdef DBG_TEMPO_GERACAO
        aux5 = omp_get_wtime();
#endif
    }

#ifdef DBG_TEMPO_GERACAO
    printf("insercao contador: %.4f\nremocao aleatoria: %.4f\nremocao gulosa: %.4f\ncalculo FO: %.4f\n", ((aux2-aux1)/(aux5-aux1)), ((aux3-aux2)/(aux5-aux1)), ((aux4-aux3)/(aux5-aux1)), ((aux5-aux4)/(aux5-aux1)));
#endif
}

void gulosidade(Solucao &s) {
    int melhorFO, melhorId, flag = 0, aux2, foAnterior, randId;
    while ((s.numConIns > maxContReal_) || (s.numFaiCob > maxFaix_)) {
        flag = 1;
        melhorFO = melhorId = -1;
        for (int i = 0; i < QTD_REM_GUL; i++) {
        	do
        		randId = rand() % numAre_;
        	while (s.vetAre[randId] == 0);
			aux2 = s.numParCob;
			s.vetAre[randId] = 0;
			calcParCob(s);
			if (s.numParCob > melhorFO) {
				melhorFO = s.numParCob;
				melhorId = randId;
			}
			s.vetAre[randId] = 1;
			s.numParCob = aux2;
        }
        s.vetAre[melhorId] = 0;
        s.numConIns -= vetArestas_[melhorId].nCon;
        s.numFaiCob -= vetArestas_[melhorId].nFai;
        s.numAreCom--;
        s.numParCob = melhorFO;
    }
    if (flag == 0 && s.numConIns != maxContReal_ && s.numFaiCob != maxFaix_) {
        while ((s.numConIns < maxContReal_) && (s.numFaiCob < maxFaix_)) {
            melhorFO = melhorId = -1;
            foAnterior = s.numParCob;
            for (int i = 0; i < QTD_REM_GUL; i++) {
            	do
            		randId = rand() % numAre_;
            	while (s.vetAre[randId] == 1);
				aux2 = s.numParCob;
				s.vetAre[randId] = 1;
				calcParCob(s);
				if (s.numParCob > melhorFO) {
					melhorFO = s.numParCob;
					melhorId = randId;
				}
				s.vetAre[randId] = 0;
				s.numParCob = aux2;
            }
            s.vetAre[melhorId] = 1;
            s.numConIns += vetArestas_[melhorId].nCon;
            s.numFaiCob += vetArestas_[melhorId].nFai;
            s.numAreCom++;
            s.numParCob = melhorFO;
        }
        if (!(((s.numConIns == maxContReal_) && (s.numFaiCob <= maxFaix_)) ||
              ((s.numConIns <= maxContReal_) && (s.numFaiCob == maxFaix_)))) {
            s.vetAre[melhorId] = 0;
            s.numConIns -= vetArestas_[melhorId].nCon;
            s.numFaiCob -= vetArestas_[melhorId].nFai;
            s.numAreCom--;
            s.numParCob = foAnterior;
        }
    } else {
        calcParCob(s);
    }
}

void epidemia() {
    Solucao melhor;
    memcpy(&melhor, &populacao[0], sizeof(populacao[0]));
    gerarPopulacao();
    memcpy(&populacao[0], &melhor, sizeof(melhor));
    ordenarPopulacao(1);
}

void ordenarPopulacao(const int &inicio) {
    Solucao escolhido;
    int j;
    for (int i = inicio; i < TAM_POP + numCro_; i++) {
        memcpy(&escolhido, &populacao[i], sizeof(populacao[i]));
        j = i - 1;
        while ((j >= 0) && (populacao[j].numParCob < escolhido.numParCob)) {
            memcpy(&populacao[j + 1], &populacao[j], sizeof(populacao[i]));
            j--;
        }
        memcpy(&populacao[j + 1], &escolhido, sizeof(populacao[i]));
    }

}

//=================================== SA =======================================

void execSA(const int p){
	Solucao melhorSolucao, sAux;
	int nIter = 0, delta;
	float temperaturaAtual = TMP_INI, x;
	memcpy(&melhorSolucao, &populacao[p], sizeof(populacao[p]));

	#ifdef DBG_TEMPO_SA
	double h1;
	h1 = omp_get_wtime();
	#endif

	while (temperaturaAtual > TMP_CNG){
		while (nIter < SAN_MAX){
			nIter++;
			memcpy(&sAux, &melhorSolucao, sizeof(melhorSolucao));
			gerVizinho(sAux);
			calcParCob(sAux);
			delta = populacao[p].numParCob - sAux.numParCob;
			if (delta < 0){
				memcpy(&populacao[p], &sAux, sizeof(melhorSolucao));
				if(sAux.numParCob > melhorSolucao.numParCob){
					memcpy(&melhorSolucao, &sAux, sizeof(melhorSolucao));
				}
			}else{
				x = (float(rand() % 101)) / 100;
				if (x < exp(float(-delta) / temperaturaAtual)){
					memcpy(&populacao[p], &sAux, sizeof(populacao[p]));
				}
			}
		}
		temperaturaAtual = TAX_RSF * temperaturaAtual;
		nIter = 0;
	}

	#ifdef DBG_TEMPO_SA
	printf("tempo SA: %.4f\n", omp_get_wtime() - h1);
	#endif

	memcpy(&populacao[p], &melhorSolucao, sizeof(melhorSolucao));
}

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
void gerVizinho(Solucao &s) {
    int pos;
	do
		pos = rand() % numAre_;
	while (s.vetAre[pos] == 1);
	s.vetAre[pos] = 1;
	s.numConIns += vetArestas_[pos].nCon;
	s.numFaiCob += vetArestas_[pos].nFai;
    //-----------------------
    // remover contadores at? que a solu??o fique vi?vel
    while ((s.numConIns > maxContReal_) || (s.numFaiCob > maxFaix_)) {
        do
            pos = rand() % numAre_;
        while (s.vetAre[pos] == 0);
        s.vetAre[pos] = 0;
        s.numConIns -= vetArestas_[pos].nCon;
        s.numFaiCob -= vetArestas_[pos].nFai;
    }
}
//------------------------------------------------------------------------------

//==============================================================================


//================================== SOLU??O ===================================

//------------------------------------------------------------------------------
void criarSolucao(Solucao &s) {
    int foAntes, bstFO, bstPos;
    memset(s.vetAre, 0, sizeof(s.vetAre));
    s.numConIns = s.numFaiCob = 0;
    s.numAreCom = numAre_;
    for (int i = 0; i < numAre_; i++) {
        s.vetAre[i] = 1;
        s.numConIns += vetArestas_[i].nCon;
        s.numFaiCob += vetArestas_[i].nFai;
    }
    s.numParCob = maxPar_;
    while ((s.numConIns > maxContReal_) || (s.numFaiCob > maxFaix_)) {
        bstFO = bstPos = -1;
        for (int i = 0; i < numAre_; i++)
            if (s.vetAre[i] == 1) {
                foAntes = s.numParCob;
                s.vetAre[i] = 0;
                calcParCob(s);
                if (s.numParCob > bstFO) {
                    bstFO = s.numParCob;
                    bstPos = i;
                }
                s.vetAre[i] = 1;
                s.numParCob = foAntes;
            }
        s.vetAre[bstPos] = 0;
        s.numConIns -= vetArestas_[bstPos].nCon;
        s.numFaiCob -= vetArestas_[bstPos].nFai;
        s.numAreCom--;
        s.numParCob = bstFO;
    }
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
void calcParCob(Solucao &s) {
    int res = 0;
#pragma omp parallel for reduction(+:res)
    for (int n = 0; n < numPar_ - 1; n++)
        res += calcPar(s, n);
    s.numParCob = maxPar_ - res;
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
int calcPar(Solucao &s, const int &n1) {
    int qtd;
    int vetBFSVis[MAX_NOS];
    int vetBFSFila[MAX_NOS];
    int res, u, v, iniFila, tamFila;
    res = 0;
    qtd = numPar_ - (n1 + 1);
    memset(vetBFSVis, 0, sizeof(vetBFSVis));
    vetBFSVis[n1] = 1;
    vetBFSFila[0] = n1;
    iniFila = 0;
    tamFila = 1;
    while (tamFila > 0) {
        u = vetBFSFila[iniFila];
        iniFila++;
        tamFila--;
        for (int i = 0; i < vetQtdAdj_[u]; i++)
            if (s.vetAre[matAreAdj_[u][i]] == 0) {
                v = matNosAdj_[u][i];
                if (!vetBFSVis[v]) {
                    if ((v > n1) && (v < numPar_)) {
                        res++;
                        qtd--;
                        if (qtd == 0)
                            return res;
                    }
                    vetBFSVis[v] = 1;
                    vetBFSFila[iniFila +
                               tamFila] = v; // circular: vetBFSFila_[(iniFila+tamFila)%MAX_NOS] = v; - n?o foi necess?rio
                    tamFila++;
                }
            }
    }
    return res;
}
//------------------------------------------------------------------------------

//==============================================================================


//============================== ENTRADA E SA?DA ===============================

//------------------------------------------------------------------------------
void lerInstancia(char *arq) {
    bool flag;
    int qtdNos, no1, no2;
    FILE *f = fopen(arq, "r");
    fscanf(f, "%d %d %d %d %d %d\n", &numNos_, &numAre_, &numPar_, &numAreVir_, &numTotAre_, &numTotFai_);
    vetIdPar_ = new int[numPar_];
    vetIdNos_ = new int[numNos_];
    for (int i = 0; i < numPar_; i++) {
        fscanf(f, "%d\n", &vetIdPar_[i]);
        vetIdNos_[i] = vetIdPar_[i];
    }
    qtdNos = numPar_;
    vetArestas_ = new Aresta[numAre_];
    for (int i = 0; i < numAre_; i++) {
        fscanf(f, "%d %d %d %d %d\n", &vetArestas_[i].id, &no1, &no2, &vetArestas_[i].nCon, &vetArestas_[i].nFai);
        flag = 1;
        for (int j = 0; j < qtdNos; j++)
            if (no1 == vetIdNos_[j]) {
                vetArestas_[i].no1 = j;
                flag = 0;
                break;
            }
        if (flag) {
            vetIdNos_[qtdNos] = no1;
            vetArestas_[i].no1 = qtdNos;
            qtdNos++;
        }
        flag = 1;
        for (int j = 0; j < qtdNos; j++)
            if (no2 == vetIdNos_[j]) {
                vetArestas_[i].no2 = j;
                flag = 0;
                break;
            }
        if (flag) {
            vetIdNos_[qtdNos] = no2;
            vetArestas_[i].no2 = qtdNos;
            qtdNos++;
        }
    }
    fclose(f);
    maxPar_ = (numPar_ * (numPar_ - 1)) / 2;
    maxCont_ = floor(numTotAre_ * ALFA / 100.0);
    maxFaix_ = floor(numTotFai_ * BETA / 100.0);
    maxContReal_ = MIN(maxCont_, maxFaix_ / 2);
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
void printarSolucao(Solucao &s) {
    printf("\n< ---------------------------- SOLUCAO ----------------------------- >\n");
    printf("Numero de postos de contagem instalados.............: %d (max = %d)\n", s.numConIns, maxCont_);
    printf("Numero de faixas cobertas...........................: %d (max = %d)\n", s.numFaiCob, maxFaix_);
    printf("Numero de pares OD cobertos (FO)....................: %d\n", s.numParCob);
    printf("\nVetor solucao: [ ");
    for (int i = 0; i < numAre_; i++)
        printf("%d ", s.vetAre[i]);
    printf("]\n\n");
    printf("-----------------------\n>>> ARESTAS PARA INSTALACAO DOS POSTOS DE CONTAGEM:\n");
    for (int i = 0; i < numAre_; i++)
        if (s.vetAre[i] == 1)
            printf("%d (%d)\n", vetArestas_[i].id, vetArestas_[i].nCon);
}

void escreverSolucao(Solucao &s, FILE *f) {
    if (f == NULL)
        f = stdout;
    fprintf(f, "\n< ---------------------------- SOLUCAO ----------------------------- >\n");
    fprintf(f, "Numero de postos de contagem instalados.............: %d (max = %d)\n", s.numConIns, maxCont_);
    fprintf(f, "Numero de faixas cobertas...........................: %d (max = %d)\n", s.numFaiCob, maxFaix_);
    fprintf(f, "Numero de pares OD cobertos (FO)....................: %d\n", s.numParCob);
    fprintf(f, "\nVetor solucao: [ ");
    for (int i = 0; i < numAre_; i++)
        fprintf(f, "%d ", s.vetAre[i]);
    fprintf(f, "]\n\n");
    fprintf(f, "-----------------------\n>>> ARESTAS PARA INSTALACAO DOS POSTOS DE CONTAGEM:\n");
    for (int i = 0; i < numAre_; i++)
        if (s.vetAre[i] == 1)
            fprintf(f, "%d (%d)\n", vetArestas_[i].id, vetArestas_[i].nCon);
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
void escreverResultado(Solucao &s, char *path) {
    FILE *f = fopen(path, "w");

	fprintf(f, "INST:%s\tALFA:%i\tBETA:%i\tNUM_EXE:%i\tMAX_TIME:%.1f\tTAM_POP:%i\tPER_REM:%i\tQTD_REM_GUL:%i\tPRC_MUT:%i\n",
			INST, ALFA, BETA, NUM_EXE, MAX_TIME, TAM_POP, PRC_REM, QTD_REM_GUL, PRC_MUT);

    fprintf(f, "< ---------------------------- PROBLEMA ---------------------------- >\n");
    fprintf(f, "Instancia...........................................: %s\n", INST);
    fprintf(f, "Reduzida............................................: NAO\n");
    fprintf(f, "Limite de contadores instalados.....................: %d (%d%% arestas sem reducao)\n", maxCont_, ALFA);
    fprintf(f, "Limite de faixas cobertas...........................: %d (%d%% faixas sem reducao)\n", maxFaix_, BETA);
    fprintf(f, "Numero de nos.......................................: %d\n", numNos_);
    fprintf(f, "Numero de arestas...................................: %d\n", numAre_);
    fprintf(f, "Numero de arestas virtuais..........................: %d\n", numAreVir_);
    fprintf(f, "Numero total de arestas (sem redu??o)...............: %d\n", numTotAre_);
    fprintf(f, "Numero total de faixas (sem redu??o)................: %d\n", numTotFai_);
    fprintf(f, "Numero de nos que definem os pares..................: %d\n", numPar_);
    fprintf(f, "Numero total de pares...............................: %d\n", maxPar_);
    fprintf(f, "\n< --------------------------- PARAMETROS --------------------------- >\n");
    fprintf(f, "Tamanho da popula??o................................: %i\n", TAM_POP);
    //fprintf(f, "Quantidade de cruzamentos por gera??o...............: %i\n", tamCem);
    fprintf(f, "Porcentagem da chance de muta??o....................: %i\n", PRC_MUT);
    fprintf(f, "N. maximo de iteracoes..............................: %d\n", MAX_ITE);
    fprintf(f, "Tempo maximo de processamento.......................: %.5f segundos!\n", MAX_TIME);
    fprintf(f, "\n< --------------------------- RESULTADOS --------------------------- >\n");
    fprintf(f, "Tempo para encontrar a melhor solucao (BST).........: %.5f segundos!\n", bstTime_);
    fprintf(f, "Tempo de execucao do metodo.........................: %.5f segundos!\n", excTime_);
    escreverSolucao(s, f);
    fclose(f);
}

void escreverResumo(int *fos, double *tempos, char *path) {
    FILE *f = fopen(path, "a");
    //fprintf(f, "%i\t%i\t%i\t%i\t%i\t", TAM_POP, PRC_CEM, PRC_MUT, PRC_ARE, PRC_GUL);
    for (int i = 0; i < NUM_EXE; i++) fprintf(f, "%i\t", fos[i]);
    for (int i = 0; i < NUM_EXE; i++) fprintf(f, "%.3f\t", tempos[i]);

    float medValores = 0, medTempos = 0;
    for (int i = 0; i < NUM_EXE; i++) {
        medValores += fos[i];
        medTempos += tempos[i];
    }
    medValores /= 5;
    medTempos /= 5;
    fprintf(f, "%.2f\t%.2f", medValores, medTempos);

    fprintf(f, "\n");

    fclose(f);
}
//------------------------------------------------------------------------------

//==============================================================================


//================================= AUXILIARES =================================

//------------------------------------------------------------------------------
void montarRede() {
    int aux;
    int vetNos[MAX_NOS];
    int vetAre[MAX_NOS];
    matNosAdj_ = new int *[numNos_];
    matAreAdj_ = new int *[numNos_];
    vetQtdAdj_ = new int[numNos_];
    for (int i = 0; i < numNos_; i++) {
        vetQtdAdj_[i] = 0;
        for (int j = 0; j < numAre_; j++) {
            aux = -1;
            if (vetArestas_[j].no1 == i)
                aux = vetArestas_[j].no2;
            else if (vetArestas_[j].no2 == i)
                aux = vetArestas_[j].no1;
            if (aux != -1) {
                vetNos[vetQtdAdj_[i]] = aux;
                vetAre[vetQtdAdj_[i]] = j;
                vetQtdAdj_[i]++;
            }
        }
        matNosAdj_[i] = new int[vetQtdAdj_[i]];
        matAreAdj_[i] = new int[vetQtdAdj_[i]];
        for (int j = 0; j < vetQtdAdj_[i]; j++) {
            matNosAdj_[i][j] = vetNos[j];
            matAreAdj_[i][j] = vetAre[j];
        }
    }
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
void limparMemoria() {
    delete[] vetIdPar_;
    delete[] vetIdNos_;
    delete[] vetArestas_;
    delete[] vetQtdAdj_;
    for (int i = 0; i < numNos_; i++)
        delete[] matNosAdj_[i];
    delete[] matNosAdj_;
    for (int i = 0; i < numNos_; i++)
        delete[] matAreAdj_[i];
    delete[] matAreAdj_;
}
//------------------------------------------------------------------------------

//==============================================================================

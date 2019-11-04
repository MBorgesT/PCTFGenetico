#include "pctfga.h"

#include <conio.h>
#include <time.h>
#include <omp.h>
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <io.h>

#define DBG // habilita o modo DEBUG (exibe na tela as melhoras na FO)
#define RND // habilita múltiplas (NUM_EXE) execuções com sementes aleatórias para a instância

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))


//============================== DADOS DE ENTRADA ==============================
char INST[50] = "AC";  // arquivo com a instância (estado)
int ALFA = 30;    // limite de contadores instalados (% do total de arestas)
int BETA = 10;    // limite de faixas cobertas (% do total de faixas)
int NUM_EXE = 5;     // número de execuções do método
double MAX_TIME = 5.0;  // tempo máximo de execução (segundos)
int MAX_ITE = 1;     // número máximo de iterações (MAX_ITE = MAX_ITE * número de arestas)
double INI_TMP = 2;     // temperatura inicial (INI_TMP = INI_TMP * fo da solução inicial)
double FRZ_TMP = 0.01;  // temperatura de congelamento
double COO_RTE = 0.975; // taxa de resfriamento
// ------------------------------- GENETICO ---------------------------------
int TAM_POP = 600;
int PRC_CEM = 5;
int PRC_MUT = 3;
int PRC_ELT = 25;
//==============================================================================


//================================== PRINCIPAL =================================
int main(int argc, char *argv[]) {
    FILE *f;
    Solucao sol;
    char arq[50], dir[50], slv[50];
    int fos[NUM_EXE];
    double tempos[NUM_EXE];
#ifdef RND
    srand(time(NULL));
#else
    NUM_EXE = 1;
#endif
    //-----------------------
    // parâmetros
    if (argc > 1) {
        strcpy(INST, argv[1]);
        ALFA = atoi(argv[2]);
        BETA = atoi(argv[3]);
        NUM_EXE = atoi(argv[4]);
        MAX_TIME = atoi(argv[5]);
        TAM_POP = atoi(argv[6]);
        PRC_CEM = atoi(argv[7]);
        PRC_MUT = atoi(argv[8]);
    }

    tamCem = (int) (TAM_POP * ((float) PRC_CEM / 100));
    tamElt = (int) (TAM_POP * ((float) PRC_ELT / 100));

    populacao = new Solucao[TAM_POP];
    sprintf(arq, "..\\Instancias\\%s.txt", INST);
    lerInstancia(arq);

    strcpy(dir, "Solucoes");
    sprintf(slv, "..\\SolucoesTesteGen\\tamPop-%i prcCem-%i prcMut-%i", TAM_POP, PRC_CEM, PRC_MUT);
    mkdir(slv);

    montarRede();

    for (int r = 1; r <= NUM_EXE; r++) {
        printf("\n\n>>> Resolvendo a instancia %s ALFA = %d e BETA = %d TAM_POP = %i PRC_CEM = %i PRC_MUT = %i - rodada %d\n\n",
               INST, ALFA, BETA, TAM_POP, PRC_CEM, PRC_MUT, r);

        execGA();

        copiarSolucao(sol, populacao[0]);
        fos[r-1] = sol.numParCob;
        tempos[r-1] = bstTime_;
    }
    //-----------------------
    // imprimir resumo
    sprintf(arq, "..\\SolucoesTesteGen\\resumo.txt");
    escreverResumo(fos, tempos, arq);
    return 0;

    /*
    //-----------------------
    // dados de entrada
    tamCem = (int) (TAM_POP * ((float) PRC_CEM / 100));
    populacao = new Solucao[TAM_POP];

    sprintf(arq, "..\\Instancias\\%s.txt", INST);
    strcpy(dir, "Solucoes");
    sprintf(slv, "..\\SolucoesTesteGen\\tamPop-%i prcCem-%i prcMut-%i prcElt-%i", TAM_POP, PRC_CEM, PRC_MUT);
    mkdir(slv);
    lerInstancia(arq);
    montarRede();
    //testEntRede();
    //-----------------------
    // executar o SA
    bstSol = 0;
    solMed = medSolAva = temMed = desvio = 0.0;
    for (int r = 1; r <= NUM_EXE; r++) {
        printf("\n\n>>> Resolvendo a instancia %s ALFA = %d e BETA = %d TAM_POP = %i PRC_CEM = %i PRC_MUT = %i - rodada %d\n\n",
               INST, ALFA, BETA, TAM_POP, PRC_CEM, PRC_MUT, r);
        execGA();
        sol = populacao[0];
        fos[r-1] = sol.numParCob;
        tempos[r-1] = bstTime_;
        sprintf(arq, "%s\\sol-%s-%da-%db-%d.txt", slv, INST, ALFA, BETA, r);
        escreverResultado(sol, arq);
        if (sol.numParCob > bstSol) {
            sprintf(arq, "%s\\bstSol-%s-%da-%db.txt", slv, INST, ALFA, BETA);
            escreverResultado(sol, arq);
            bstSol = sol.numParCob;
        }
        solMed += sol.numParCob;
        temMed += bstTime_;
        f = fopen("saida-full.txt", "at");
        fprintf(f, "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t\t%d\t%d\t%d\t%.2f\t\t%.2f\n",
                INST, ALFA, BETA, numNos_, numAre_, numTotAre_, numTotFai_, numPar_, solAva_, sol.numConIns,
                sol.numFaiCob, sol.numParCob, bstTime_, excTime_);
        fclose(f);
    }
    //-----------------------
    // calcular as médias
    solMed = solMed / (double) NUM_EXE;
    medSolAva = medSolAva / (double) NUM_EXE;
    temMed = temMed / (double) NUM_EXE;
    desvio = ((bstSol - solMed) / bstSol) * 100;
    f = fopen("saida.txt", "at");
    fprintf(f, "%s\t%d\t%d\t%d\t%.2f\t%.2f\t%.2f\t%.2f\n", INST, ALFA, BETA, bstSol, solMed, desvio, temMed, medSolAva);
    fclose(f);
    //-----------------------
    // imprimir resumo
    sprintf(arq, "..\\SolucoesTesteGen\\resumo.txt");
    escreverResumo(fos, tempos, arq);
    //-----------------------
    // limpar memória
    limparMemoria();
    if (argc == 1) {
        printf("\n\n>>> Pressione ENTER para encerrar: ");
        _getch();
    }
    return 0;
    */
}
//==============================================================================
//==================================== GA ======================================

void execGA() {
    int melhorFO, flag, aux, qtdGen = 0, qtdEpi = 0;

    gerarPopulacao();
    //---- ordena a população toda, o que não é necessário depois
    Solucao escolhido;
    int j;

    for (int i = 1; i < TAM_POP; i++) {
        copiarSolucao(escolhido, populacao[i]);
        j = i - 1;

        while ((j >= 0) && (populacao[j].numParCob < escolhido.numParCob)) {
            copiarSolucao(populacao[j + 1], populacao[j]);
            j--;
        }

        copiarSolucao(populacao[j + 1], escolhido);
    }
    //----

    melhorFO = populacao[0].numParCob;

    clock_t h1, h2;
    h1 = clock();
    h2 = clock();

    //---- teste de tempo
    /*
    clock_t h3;
    double cross, ord, epi, total;

    h3 = clock();
    crossover();
    cross = (clock() - h3) / (float) CLOCKS_PER_SEC;


    for (int i = 0; i < TAM_POP; i++) printf("%i: %i\n", i, populacao[i].numParCob);

    h3 = clock();
    ordenarPopulacao();
    ord = (clock() - h3) / (float) CLOCKS_PER_SEC;

    h3 = clock();
    aux = populacao[0].numParCob;
    flag = 1;
    for (int i = 1; i < TAM_POP - tamCem; i++) {
        if (populacao[i].numParCob != aux) {
            flag = 0;
            break;
        }
    }
    if (flag == 1) {
        epidemia();
        qtdEpi++;
    }
    epi = (clock() - h3) / (float) CLOCKS_PER_SEC;

    h3 = clock();
    total = (h3 - h1) / (float) CLOCKS_PER_SEC;

    printf("cross: %.6f\tord: %.6f\tepi: %.6f\n", cross/total, ord/total, epi/total);
    for (int i = 0; i < TAM_POP; i++) printf("%i: %i\n", i, populacao[i].numParCob);

    if (populacao[0].numParCob > melhorFO) {
        melhorFO = populacao[0].numParCob;
        bstTime_ = (clock() - h1) / (float) CLOCKS_PER_SEC;
    }

    h2 = clock();
    qtdGen++;
    */
    //----

    while (((double) ((h2 - h1)) / CLOCKS_PER_SEC) < MAX_TIME) {
        crossover();
        ordenarPopulacao();

        aux = populacao[0].numParCob;
        flag = 1;
        for (int i = 1; i < TAM_POP - tamCem; i++) {
            if (populacao[i].numParCob != aux) {
                flag = 0;
                break;
            }
        }
        if (flag == 1) {
            epidemia();
            qtdEpi++;
        }

        if (populacao[0].numParCob > melhorFO) {
            melhorFO = populacao[0].numParCob;
            bstTime_ = (clock() - h1) / (float) CLOCKS_PER_SEC;
        }

        h2 = clock();
        qtdGen++;
    }

    printf("bstSol: %i\tbstTempo: %.4f\tqtd gen: %i\tqtd epi: %i\tdiv: %.5f\n", populacao[0].numParCob, bstTime_,
           qtdGen, qtdEpi, (float) qtdEpi / (float) qtdGen);
}

void gerarPopulacao() {
    for (int i = 0; i < TAM_POP; i++) {
        heuAleGA(populacao[i]);
    }
}

void heuAleGA(Solucao &s) {
    int aux, i = 0;
    memset(s.vetAre, 0, sizeof(s.vetAre));
    s.numConIns = s.numFaiCob = 0;
    s.numAreCom = 0;
    do {
        aux = rand() % numTotAre_;
        s.vetAre[aux] = 1;
        s.numConIns += vetArestas_[aux].nCon;
        s.numFaiCob += vetArestas_[aux].nFai;
        s.vetAreIds[i] = aux;
        i++;
    } while ((s.numConIns < maxContReal_) && (s.numFaiCob < maxFaix_));

    s.vetAre[aux] = 0;
    s.numConIns -= vetArestas_[aux].nCon;
    s.numFaiCob -= vetArestas_[aux].nFai;
    s.numAreCom = i;

    calcParCob(s);
}

void crossover() {
    int a, b;
    for (int i = TAM_POP - tamCem; i < TAM_POP; i++) {
        a = rand() % (TAM_POP - tamCem);
        b = rand() % (TAM_POP - tamCem);
        gerarFilho(populacao[i], populacao[a], populacao[b]);
    }
}

void gerarFilho(Solucao &filho, Solucao &pai, Solucao &mae) {
    int aux = (rand() % ((pai.numAreCom + mae.numAreCom) / 2) - 1) +
              1;  // possivel problema com array out of limits na hora de puxar

    memset(filho.vetAre, 0, sizeof(filho.vetAre));
    filho.numConIns = filho.numFaiCob = 0;
    filho.numAreCom = 0;

    for (int i = 0; i < aux; i++) {
        filho.vetAre[pai.vetAreIds[i]] = 1;
        filho.numConIns += vetArestas_[pai.vetAreIds[i]].nCon;
        filho.numFaiCob += vetArestas_[pai.vetAreIds[i]].nFai;
        filho.vetAreIds[i] = pai.vetAreIds[i];
        filho.numAreCom++;
    }
    for (int i = aux; i < mae.numAreCom; i++) {
        filho.vetAre[mae.vetAreIds[i]] = 1;
        filho.numConIns += vetArestas_[mae.vetAreIds[i]].nCon;
        filho.numFaiCob += vetArestas_[mae.vetAreIds[i]].nFai;
        filho.vetAreIds[i] = mae.vetAreIds[i];
        filho.numAreCom++;
    }

    if (rand() % 100 < PRC_MUT) {
        gerVizinho(filho);
    }

    int pos;
    while ((filho.numConIns > maxContReal_) || (filho.numFaiCob > maxFaix_)) {
        do
            pos = rand() % numAre_;
        while (filho.vetAre[pos] == 0);
        filho.vetAre[pos] = 0;
        filho.numConIns -= vetArestas_[pos].nCon;
        filho.numFaiCob -= vetArestas_[pos].nFai;
    }

    calcParCob(filho);
}

void ordenarPopulacao() {
    Solucao escolhido;
    int j;

    for (int i = TAM_POP - tamCem; i < TAM_POP; i++) {
        copiarSolucao(escolhido, populacao[i]);
        j = i - 1;

        while ((j >= 0) && (populacao[j].numParCob < escolhido.numParCob)) {
            copiarSolucao(populacao[j + 1], populacao[j]);
            j--;
        }

        copiarSolucao(populacao[j + 1], escolhido);
    }
}

void copiarSolucao(Solucao &destino, Solucao &origem) {
    memcpy(&destino, &origem, sizeof(origem));
}

void epidemia() {
    Solucao melhor;
    copiarSolucao(melhor, populacao[0]);

    gerarPopulacao();
    copiarSolucao(populacao[0], melhor);
    //---- ordena a população toda, o que não é necessário depois
    int j;
    Solucao escolhido;
    for (int i = 1; i < TAM_POP; i++) {
        copiarSolucao(escolhido, populacao[i]);
        j = i - 1;

        while ((j >= 0) && (populacao[j].numParCob < escolhido.numParCob)) {
            copiarSolucao(populacao[j + 1], populacao[j]);
            j--;
        }

        copiarSolucao(populacao[j + 1], escolhido);
    }
    //----
}

//==================================== CS ======================================

//------------------------------------------------------------------------------
void execSA(Solucao &s) {
    clock_t hI, hF;
    Solucao sA, sV;
    double x, iniTmp, tmp;
    int iter, itMax, foAux;
    iniTmp = MAX(10, INI_TMP * s.numParCob);
    itMax = MAX_ITE * numAre_;
    hI = clock();
    criarSolucao(s);
    hF = clock();
    bstTime_ = ((double) hF - hI) / CLOCKS_PER_SEC;
    memcpy(&sA, &s, sizeof(s));
    foIni_ = s.numParCob;
#ifdef DBG
    printf("\n>>> Sol. Ini.: %d\tTempo: %.2f\n\n", foIni_, bstTime_);
#endif
    solAva_ = 1;
    excTime_ = 0.0;
    while (excTime_ < MAX_TIME) {
        tmp = iniTmp;
        while (tmp > FRZ_TMP) {
            for (int i = 0; i < itMax; i++) {
                solAva_++;
                memcpy(&sV, &sA, sizeof(sA));
                gerVizinho(sV);
                if (sV.numParCob > sA.numParCob) {
                    memcpy(&sA, &sV, sizeof(sV));
                    if (sV.numParCob > s.numParCob) {
                        hF = clock();
                        bstTime_ = ((double) hF - hI) / CLOCKS_PER_SEC;
                        memcpy(&s, &sV, sizeof(sV));
#ifdef DBG
                        printf(">>> Best Sol: %d\tTempo: %.2f\n", s.numParCob, bstTime_);
#endif
                    }
                } else {
                    x = rand() % 1001;
                    x = x / 1000.0;
                    if (x < expl(-(sA.numParCob - sV.numParCob) / tmp))
                        memcpy(&sA, &sV, sizeof(sV));
                }
                hF = clock();
                excTime_ = ((double) hF - hI) / CLOCKS_PER_SEC;
                if (excTime_ >= MAX_TIME)
                    goto FIM;
            }
            tmp = COO_RTE * tmp;
        }
    }
    hF = clock();
    excTime_ = ((double) hF - hI) / CLOCKS_PER_SEC;
    FIM :;
    foFin_ = s.numParCob;
#ifdef DBG
    printf("\n>>> Sol. Fin.: %d\tTempo: %.2f\tSol. Ava.: %d\n", foFin_, excTime_, solAva_);
#endif
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
void gerVizinho(Solucao &s) {
    int pos;
    pos = rand() % numAre_;
    if (s.vetAre[pos] == 1) {
        s.vetAre[pos] = 0;
        s.numConIns -= vetArestas_[pos].nCon;
        s.numFaiCob -= vetArestas_[pos].nFai;
    } else {
        s.vetAre[pos] = 1;
        s.numConIns += vetArestas_[pos].nCon;
        s.numFaiCob += vetArestas_[pos].nFai;
        //-----------------------
        // remover contadores até que a solução fique viável
        while ((s.numConIns > maxContReal_) || (s.numFaiCob > maxFaix_)) {
            do
                pos = rand() % numAre_;
            while (s.vetAre[pos] == 0);
            s.vetAre[pos] = 0;
            s.numConIns -= vetArestas_[pos].nCon;
            s.numFaiCob -= vetArestas_[pos].nFai;
        }
    }
}
//------------------------------------------------------------------------------

//==============================================================================


//================================== SOLUÇÃO ===================================

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
                               tamFila] = v; // circular: vetBFSFila_[(iniFila+tamFila)%MAX_NOS] = v; - não foi necessário
                    tamFila++;
                }
            }
    }
    return res;
}
//------------------------------------------------------------------------------

//==============================================================================


//============================== ENTRADA E SAÍDA ===============================

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
    fprintf(f, "< ---------------------------- PROBLEMA ---------------------------- >\n");
    fprintf(f, "Instancia...........................................: %s\n", INST);
    fprintf(f, "Reduzida............................................: NAO\n");
    fprintf(f, "Limite de contadores instalados.....................: %d (%d%% arestas sem reducao)\n", maxCont_, ALFA);
    fprintf(f, "Limite de faixas cobertas...........................: %d (%d%% faixas sem reducao)\n", maxFaix_, BETA);
    fprintf(f, "Numero de nos.......................................: %d\n", numNos_);
    fprintf(f, "Numero de arestas...................................: %d\n", numAre_);
    fprintf(f, "Numero de arestas virtuais..........................: %d\n", numAreVir_);
    fprintf(f, "Numero total de arestas (sem redução)...............: %d\n", numTotAre_);
    fprintf(f, "Numero total de faixas (sem redução)................: %d\n", numTotFai_);
    fprintf(f, "Numero de nos que definem os pares..................: %d\n", numPar_);
    fprintf(f, "Numero total de pares...............................: %d\n", maxPar_);
    fprintf(f, "\n< --------------------------- PARAMETROS --------------------------- >\n");
    fprintf(f, "Tamanho da população................................: %i\n", TAM_POP);
    fprintf(f, "Quantidade de cruzamentos por geração...............: %i\n", tamCem);
    fprintf(f, "Porcentagem da chance de mutação....................: %i\n", PRC_MUT);
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
    fprintf(f, "%i\t%i\t%i\t", TAM_POP, PRC_CEM, PRC_MUT);
    for (int i = 0; i < NUM_EXE; i++) fprintf(f, "%i\t", fos[i]);
    for (int i = 0; i < NUM_EXE; i++) fprintf(f, "%.3f\t", tempos[i]);

    float medValores = 0, medTempos = 0;
    for (int i = 0; i < NUM_EXE; i++){
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
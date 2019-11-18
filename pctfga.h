#ifndef uPCTFH
#define uPCTFH

#include <iostream>

#define MAX_NOS 2500 // número máximo de nós - baseado no tamanho das instâncias
#define MAX_ARE 2600 // número máximo de arestas - baseado no tamanho das instâncias
#define MAX_PAR 900  // número máximo de pares O-D - baseado no tamanho das instâncias
// ----- Genetico

// ------------------------------- ESTRUTURAS -------------------------------
typedef struct tAresta {
    int id;   // identificador real da aresta
    int no1;  // id (posição no vetor de nós) do primeiro nó da aresta
    int no2;  // id (posição no vetor de nós) do segundo nó da aresta
    int nCon; // número de contadores da aresta
    int nFai; // número de faixas da aresta
} Aresta;

typedef struct tSolucao {
    int vetAre[MAX_ARE]; // vetor binário de arestas (1 - contador instalado; 0 - caso contrário)
    int vetAreIds[50]; // para usar na geração de filhos do GA; obs: alterar tamanho do vetor
    int numAreCom; // número de arestas COM contador
    int numConIns; // número de contadores instalados
    int numFaiCob; // número de faixas cobertas
    int numParCob; // número de pares cobertos
} Solucao;
//---------------------------------------------------------------------------


//---------------------------- VARIÁVEIS GLOBAIS ----------------------------
// ----- Genetico
Solucao *populacao;
int tamCem;
int tamGul;
int tamElt;
// ----- Dados de entrada
int numPar_; // número de pares O-D 
int numNos_; // número de nós na rede
int numAre_; // número de arestas na rede
int numAreVir_; // número de arestas virtuais
int numTotAre_; // número total de arestas (instância completa)
int numTotFai_; // número total de faixas das arestas (instância completa)
int *vetIdPar_; // vetor com o id real do nós usados como pares O-D
int *vetIdNos_; // vetor com o id real dos nós
Aresta *vetArestas_; // vetor de arestas
// ----- Rede
int *vetQtdAdj_;  // vetor com a quantidade de nós (e arestas) adjacentes a cada nó
int **matNosAdj_; // matriz com os nós adjacentes a cada nó
int **matAreAdj_; // matriz com as arestas adjacentes a cada nó
// ----- Resultados
int foIni_, foFin_;        // função objetivo das soluções inicial e final
double bstTime_, excTime_; // tempo para encontrar a melhor solução e tempo total
int solAva_;              // número de soluções avaliadas pelo CS
// ----- Variáveis auxiliares
int maxPar_;      // número máximo de pares cobertos
int maxCont_;     // número máximo de contadores a serem instalados  
int maxFaix_;     // número máximo de faixas a serem cobertas
int maxContReal_; // número máximo REAL de contadores (definido com base no limite de contadores e faixas)
//---------------------------------------------------------------------------


//--------------------------------- MÉTODOS ---------------------------------
// ------------ Genetico
void execGA();

void gerarPopulacao();

void heuAleGA(Solucao &s);

void crossover();

void gerarFilho(Solucao &filho, Solucao &pai, Solucao &mae);

void gulosidade(Solucao &s);

void ordenarPopulacao();

void copiarSolucao(Solucao &destino, Solucao &origem);

void epidemia();

// ------------ SA
void execSA(Solucao &s);

void gerVizinho(Solucao &s);

// ------------ Solução
void criarSolucao(Solucao &s);

void calcParCob(Solucao &s);

int calcPar(Solucao &s, const int &n1);

// ------------ Entrada e Saída
void lerInstancia(char *arq);

void printarSolucao(Solucao &s);

void escreverSolucao(Solucao &s, FILE *f);

void escreverResultado(Solucao &s, char *path);

void escreverResumo(int* fos, double* tempos, char* path);

// ------------ Auxiliares
void montarRede();

void limparMemoria();

#endif
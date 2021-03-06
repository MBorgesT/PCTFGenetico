#ifndef uPCTFH
#define uPCTFH

#include <iostream>

#define MAX_NOS 2500 // n�mero m�ximo de n�s - baseado no tamanho das inst�ncias
#define MAX_ARE 2600 // n�mero m�ximo de arestas - baseado no tamanho das inst�ncias
#define MAX_PAR 900  // n�mero m�ximo de pares O-D - baseado no tamanho das inst�ncias
// ----- Genetico
#define MAX_POP 2000 // geraldo

// ------------------------------- ESTRUTURAS -------------------------------
typedef struct tAresta {
    int id;   // identificador real da aresta
    int no1;  // id (posi��o no vetor de n�s) do primeiro n� da aresta
    int no2;  // id (posi��o no vetor de n�s) do segundo n� da aresta
    int nCon; // n�mero de contadores da aresta
    int nFai; // n�mero de faixas da aresta
} Aresta;

typedef struct tSolucao {
    int vetAre[MAX_ARE]; // vetor bin�rio de arestas (1 - contador instalado; 0 - caso contr�rio)
    //int vetAreIds[50]; // para usar na gera��o de filhos do GA; obs: alterar tamanho do vetor
    int numAreCom; // n�mero de arestas COM contador
    int numConIns; // n�mero de contadores instalados
    int numFaiCob; // n�mero de faixas cobertas
    int numParCob; // n�mero de pares cobertos
} Solucao;
//---------------------------------------------------------------------------


//---------------------------- VARI�VEIS GLOBAIS ----------------------------
// ----- Genetico
Solucao populacao[MAX_POP]; // geraldo
int limRemocao_; // geraldo
int tamElt_;
int limSup_;
// ----- Dados de entrada
int numPar_; // n�mero de pares O-D 
int numNos_; // n�mero de n�s na rede
int numAre_; // n�mero de arestas na rede
int numAreVir_; // n�mero de arestas virtuais
int numTotAre_; // n�mero total de arestas (inst�ncia completa)
int numTotFai_; // n�mero total de faixas das arestas (inst�ncia completa)
int *vetIdPar_; // vetor com o id real do n�s usados como pares O-D
int *vetIdNos_; // vetor com o id real dos n�s
Aresta *vetArestas_; // vetor de arestas
// ----- Rede
int *vetQtdAdj_;  // vetor com a quantidade de n�s (e arestas) adjacentes a cada n�
int **matNosAdj_; // matriz com os n�s adjacentes a cada n�
int **matAreAdj_; // matriz com as arestas adjacentes a cada n�
// ----- Resultados
int bstFO_;
double bstTime_, excTime_; // tempo para encontrar a melhor solu��o e tempo total
// ----- Vari�veis auxiliares
int maxPar_;      // n�mero m�ximo de pares cobertos
int maxCont_;     // n�mero m�ximo de contadores a serem instalados  
int maxFaix_;     // n�mero m�ximo de faixas a serem cobertas
int maxContReal_; // n�mero m�ximo REAL de contadores (definido com base no limite de contadores e faixas)
int numCro_;
//---------------------------------------------------------------------------


//--------------------------------- M�TODOS ---------------------------------
// ------------ Genetico
// geraldo
void execGA(const int r);
void gerarFilho(const int &p1,const int &p2,int f);
void gerarPopulacao();
void gulosidade(Solucao &s);
void epidemia();
void ordenarPopulacao(const int &inicio);
void gerVizinhoGul(const int p);

// ------------ SA
void execSA(const int p);

void gerVizinho(Solucao &s);

// ------------ Solu��o
void criarSolucao(Solucao &s);

void calcParCob(Solucao &s);

int calcPar(Solucao &s, const int &n1);

// ------------ Entrada e Sa�da
void lerInstancia(char *arq);

void printarSolucao(Solucao &s);

void escreverSolucao(Solucao &s, FILE *f);

void escreverResultado(Solucao &s, char *path);

void escreverResumo(int* fos, double* tempos, char* path);

// ------------ Auxiliares
void montarRede();

void limparMemoria();

#endif

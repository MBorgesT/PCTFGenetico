#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <time.h>

#include "main.h"

#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))

int main() {
    srand(time(NULL));

    lerDados("entradas.txt");
    testarDados("teste.txt");

    for (int i = 0; i < 10; i++) {
        algGenetico();
    }

    return 0;
}


void lerDados(char *arq) {
    int aux;

    FILE *f = fopen(arq, "r");
    fscanf(f, "%i %i", &qtdCalLin, &qtdCalCol);
    fscanf(f, "%i %i %i", &qtdProf, &qtdHorLin, &qtdHorCol);

    for (int i = 0; i < qtdCalLin; i++) {
        for (int j = 0; j < qtdCalCol; j++) {
            fscanf(f, "%d", &aux);
            calendario.matriz[i][j] = aux;
        }
    }

    for (int i = 0; i < qtdProf; i++) {
        for (int j = 0; j < qtdHorLin; j++) {
            for (int k = 0; k < qtdHorCol; k++) {
                fscanf(f, "%d", &aux);
                listaProfessores[i].horario[j][k] = aux;
            }
        }
    }

    for (int i = 0; i < qtdProf; i++) {
        for (int j = 0; j < 2; j++) {
            for (int k = 0; k < 5; k++) {
                fscanf(f, "%d", &aux);
                listaProfessores[i].tercTunro[j][k] = aux;
            }
        }
    }

    for (int i = 0; i < qtdProf; i++){
        for (int j = 0; j < 6; j++){
            for (int k = 0; k < 5; k++){
                horarioProfSoma[j][k] += listaProfessores[i].horario[j][k];
            }
        }

        for (int j = 0; j < 2; j++){
            for (int k = 0; k < 5; k++){
                horarioTercTurnoSoma[j][k] += listaProfessores[i].tercTunro[j][k];
            }
        }
    }

    fclose(f);
}

void testarDados(char *arq) {
    FILE *f;
    if (arq == "")
        f = stdout;
    else
        f = fopen(arq, "w");

    fprintf(f, "%i %i\n", qtdCalLin, qtdCalCol);
    fprintf(f, "%i %i %i\n\n", qtdProf, qtdHorLin, qtdHorCol);

    for (int i = 0; i < qtdCalLin; i++) {
        for (int j = 0; j < qtdCalCol; j++) {
            fprintf(f, "%i ", calendario.matriz[i][j]);
        }
        fprintf(f, "\n");
    }
    fprintf(f, "\n");

    for (int i = 0; i < qtdProf; i++) {
        for (int j = 0; j < qtdHorLin; j++) {
            for (int k = 0; k < qtdHorCol; k++) {
                fprintf(f, "%i ", listaProfessores[i].horario[j][k]);
            }
            fprintf(f, "\n");
        }
        fprintf(f, "\n");
    }

    for (int i = 0; i < qtdProf; i++) {
        for (int j = 0; j < 2; j++) {
            for (int k = 0; k < 5; k++) {
                fprintf(f, "%i ", listaProfessores[i].tercTunro[j][k]);
            }
            fprintf(f, "\n");
        }
        fprintf(f, "\n");
    }


    if (arq == "")
        fclose(f);
}

void zerarTotais(){
    for (int i = 0; i < qtdProf; i++){
        listaProfessores[i].qtdAulasPerdidas = 0;
        listaProfessores[i].qtdReuTercTurno = 0;
        listaProfessores[i].qtdAulasNoiteAntesReuniao = 0;
    }
    qtdAulasPerdidasTotal = 0;
    qtdReuTerTurnoTotal = 0;
    qtdAulasNoiteAntesReuniTotal = 0;
}

void heuCon(Solucao &s){
    int vet[9];
    for (int i = 0; i < 9; i++){
        vet[i] = rand() % 126;
        if ((vet[i] % 7 == 0) || (vet[i] % 7 == 6))
            i--;
    }

    int escolhido, j;

    for (int i = 1; i < 9; i++) {
        escolhido = vet[i];
        j = i - 1;

        while ((j >= 0) && (vet[j] > escolhido)) {
            vet[j + 1] = vet[j];
            j--;
        }

        vet[j + 1] = escolhido;
    }

    for (int i = 0; i < 9; i++){
        s.reunioes[i].data[0] = (int)(vet[i]/7);
        s.reunioes[i].data[1] = vet[i] % 7;
        s.reunioes[i].horario = rand() % 4;
    }

    calFO(s);
}

void calTotais(Solucao &s) {
    for (int i = 0; i < 9; i++) { //reunioes
        for (int j = 0; j < qtdProf; j++) { //professores
            //Aulas perdidas
            if (listaProfessores[j].horario[s.reunioes[i].horario][s.reunioes[i].data[1] - 1]) { //-1 porque a matriz de calendário conta sábado e domingo, sendo que a de horario de professores não
                listaProfessores[j].qtdAulasPerdidas++;
                qtdAulasPerdidasTotal++;
            }

            //Aulas de noite antes da reuniao
            if (s.reunioes[i].data[1] != 1) {
                if ((listaProfessores[j].horario[4][s.reunioes[i].data[1] - 2] ||
                     listaProfessores[j].horario[5][s.reunioes[i].data[1] - 2]) && (s.reunioes[i].horario == 0 || s.reunioes[i].horario == 1)) { //-2 porque a matriz de calendário conta sábado e domingo, sendo que a de horario de professores não
                    listaProfessores[j].qtdAulasNoiteAntesReuniao++;
                    qtdAulasNoiteAntesReuniTotal++;
                }
            }

            //Terceiro turno
            if (listaProfessores[j].tercTunro[(int)(s.reunioes[i].horario/2)][s.reunioes[i].data[1] - 1]) {
                listaProfessores[j].qtdReuTercTurno++;
                qtdReuTerTurnoTotal++;
            }
        }
    }
}

float calRestricao2(){
    float MAP = qtdAulasPerdidasTotal / qtdProf;
    float resultado = 0;
    for (int i = 0; i < qtdProf; i++){
        resultado += abs(listaProfessores[i].qtdAulasPerdidas - MAP);
    }
    return resultado;
}

float calRestricao5(){
    float resultado = 0;
    for (int i = 0; i < qtdProf; i++){
        resultado += MAX(0, listaProfessores[i].qtdAulasPerdidas - 3);
    }
    return resultado;
}

float calRestricao6(Solucao &s){
    float resultado = 0;
    for (int i = 0; i < 4; i++){
        resultado += abs(calQtdReuMes(s, i) - 2);
    }
    resultado += abs(calQtdReuMes(s, 4) - 1);
    return resultado;
}

float calRestricao7(Solucao &s){
    float resultado = 0;
    for (int i = 1; i < 9; i++){
        int dia1 = (s.reunioes[i-1].data[0] * 7) + s.reunioes[i-1].data[1];
        int dia2 = (s.reunioes[i].data[0] * 7) + s.reunioes[i].data[1];
        resultado += MAX(0, 10 - dia2 + dia1);
    }
    return resultado;
}

int calQtdReuMes(Solucao &s, int op){
    int qtd = 0;
    switch(op){
        case 0:
            for (int i = 0; i < 9; i++){
                int dia = (s.reunioes[i].data[0] * 7) + s.reunioes[i].data[1] + 1;
                if(dia <= 22){
                    qtd++;
                }
            }
            break;
        case 1:
            for (int i = 0; i < 9; i++){
                int dia = (s.reunioes[i].data[0] * 7) + s.reunioes[i].data[1] + 1;
                if((dia > 22 ) && (dia <= 52)){
                    qtd++;
                }
            }
            break;
        case 2:
            for (int i = 0; i < 9; i++){
                int dia = (s.reunioes[i].data[0] * 7) + s.reunioes[i].data[1] + 1;
                if((dia > 52 ) && (dia <= 83)){
                    qtd++;
                }
            }
            break;
        case 3:
            for (int i = 0; i < 9; i++){
                int dia = (s.reunioes[i].data[0] * 7) + s.reunioes[i].data[1] + 1;
                if((dia > 83) && (dia <= 113)){
                    qtd++;
                }
            }
            break;
        case 4:
            for (int i = 0; i < 9; i++){
                int dia = (s.reunioes[i].data[0] * 7) + s.reunioes[i].data[1] + 1;
                if((dia > 113 ) && (dia <= 126)){
                    qtd++;
                }
            }
            break;
        default:
            printf("Parametro incorreto na funcao calQtdReuMes");
            break;
    }
    return qtd;
}

void calFO(Solucao &s){
    zerarTotais();
    calTotais(s);

    s.fo = 0;
    s.fo = (pesos[0] * qtdAulasPerdidasTotal) +
            (pesos[1] * calRestricao2()) +
            (pesos[2] * qtdAulasNoiteAntesReuniTotal) +
            (pesos[3] * qtdReuTerTurnoTotal) +
            (pesos[4] * calRestricao5()) +
            (pesos[5] * calRestricao6(s)) +
            (pesos[6] * calRestricao7(s));

    //é adicionado 100000 à FO caso alguma reunião esteja marcada para um dia não letivo
    for (int i = 0; i < 9; i++){
        if (calendario.matriz[s.reunioes[i].data[0]][s.reunioes[i].data[1]] == 0)
            s.fo += 100000;
    }

}

void escreverResultado(Solucao &s){
    zerarTotais();
    calTotais(s);

    printf("\n---------------------------------------------\n\nFO: %i\n\n", s.fo);
    printf("Reunioes: \n");
    for (int i = 0; i < 9; i++){
        int dia = (s.reunioes[i].data[0] * 7) + s.reunioes[i].data[1] + 1;
        if(dia <= 22){
            printf("%i de marco\t", dia + 10);
        }else if((dia > 22 ) && (dia <= 52)){

            printf("%i de abril\t", dia - 22);

        }else if((dia > 52 ) && (dia <= 83)){
            printf("%i de maio\t", dia - 52);

        }else if((dia > 83 ) && (dia <= 113)){
            printf("%i de junho\t", dia - 83);

        }else if((dia > 113 ) && (dia <= 126)){
            printf("%i de julho\t", dia - 113);

        }

        switch (s.reunioes[i].data[1]){
            case 1:
                printf("Segunda feira\t");
                break;
            case 2:
                printf("Terca feira\t\t");
                break;
            case 3:
                printf("Quarta feira\t");
                break;
            case 4:
                printf("Quinta feira\t");
                break;
            case 5:
                printf("Sexta feira\t\t");
                break;
        }

        printf("Horario: %i\n", s.reunioes[i].horario);
    }
    printf("\n");

    for (int i = 0; i < qtdProf; i++){
        printf("Prof %i\t\tAulas Perdidas %i\tReu Ter Turno: %i\tReuni depois de aula a noite: %i\n", i, listaProfessores[i].qtdAulasPerdidas, listaProfessores[i].qtdReuTercTurno, listaProfessores[i].qtdAulasNoiteAntesReuniao);
    }

    printf("\nAulas medias Perdias: %.2f\tAulas perdidas total: %i\tReuni Ter Turno total: %i\tReuni depois de aula a noite total: %i\n",((float)qtdAulasPerdidasTotal)/16,  qtdAulasPerdidasTotal, qtdReuTerTurnoTotal, qtdAulasNoiteAntesReuniTotal);

    printf("\n---------------------------\n");
}

void escreverResultadoSimp(Solucao &s){
    printf("\n-------------------------------\n");
    printf("FO: %i\n", s.fo);

    for (int i = 0; i < 9; i++){
        printf("%i: dia %i\t horario %i\n", i, ((s.reunioes[i].data[0] * 7) + s.reunioes[i].data[1]), s.reunioes[i].horario);
    }
}

void copiarSolucao(Solucao &destino, Solucao &origem){
    memcpy(&destino, &origem, sizeof(origem));
}

void algGenetico(){
    int aux, flag, melhorFO;

    gerarPop();
    ordenarPop();
    melhorFO = pop[0].fo;

    clock_t h1, h2, hMaxFO;
    h1 = clock();
    h2 = clock();
    hMaxFO = clock();

    for (int i = 0; ((double)((h2 - h1))/CLOCKS_PER_SEC) < 60; i++){
        crossover();
        ordenarPop();
        //printf("%i: %i\n", i, pop[0].fo);

        aux = pop[0].fo;
        flag = 1;
        for (int i = 1; i < TAM_POP - QTD_CRUZ; i++){
            if (pop[i].fo != aux){
                flag = 0;
                break;
            }
        }
        if (flag == 1) {
            epidemia();
        }

        if (pop[0].fo < melhorFO){
            melhorFO = pop[0].fo;
            hMaxFO = clock();
        }

        h2 = clock();
    }

    printf("Melhor fo: %i\nTempo: %.3f\n---------------------\n", melhorFO, (((double)hMaxFO) - ((double)h1))/ CLOCKS_PER_SEC);
}

void crossover(){
    int a, b;
    for (int i = TAM_POP - QTD_CRUZ; i < TAM_POP; i++){
        do{
            a = rand() % (TAM_ELT);
            b = rand() % (TAM_POP - QTD_CRUZ);
        }while(a == b);
        if (rand() % 2){
            gerarFilho(i, a, b);
        }else{
            gerarFilho(i, b, a);
        }
    }
}

void gerarFilho(int i, int a, int b){
    int aux = (rand() % 8) + 1;
    for (int j = 0; j < aux; j++){
        pop[i].reunioes[j].data[0] = pop[a].reunioes[j].data[0];
        pop[i].reunioes[j].data[1] = pop[a].reunioes[j].data[1];
        pop[i].reunioes[j].horario = pop[a].reunioes[j].horario;
    }
    for (int j = aux; j < 9; j++){
        pop[i].reunioes[j].data[0] = pop[b].reunioes[j].data[0];
        pop[i].reunioes[j].data[1] = pop[b].reunioes[j].data[1];
        pop[i].reunioes[j].horario = pop[b].reunioes[j].horario;
    }

    if (rand() % 100 < MUTACAO){
        mutacao(i);
    }

    calFO(pop[i]);
}

void gerarPop(){
    for (int i = 0; i < TAM_POP; i++){
        heuCon(pop[i]);
    }
}

void ordenarPop(){
    Solucao escolhido;
    int j;

    for (int i = 1; i < TAM_POP; i++) {
        copiarSolucao(escolhido, pop[i]);
        j = i - 1;

        while ((j >= 0) && (pop[j].fo > escolhido.fo)) {
            copiarSolucao(pop[j+1], pop[j]);
            j--;
        }

        copiarSolucao(pop[j+1], escolhido);
    }

}

void epidemia(){
    Solucao melhor;
    copiarSolucao(melhor, pop[0]);

    gerarPop();
    copiarSolucao(pop[0], melhor);
    ordenarPop();
}

void mutacao(int i){
    int aux = rand() % 9;
    if (rand() % 2){
        int dia = (pop[i].reunioes[aux].data[0] * 7) + pop[i].reunioes[aux].data[1];
        dia++;
        pop[i].reunioes[aux].data[0] = (int)(dia/7);
        pop[i].reunioes[aux].data[1] = dia % 7;
    }else{
        pop[i].reunioes[aux].horario = (pop[i].reunioes[aux].horario + 1) % 4;
    }
}

void mutacaoVizinhos(int i){
    int aux = rand() % 9;

    int melhorHorario = 0, valorAux, valorAux2;
    valorAux = (10 * horarioProfSoma[0][pop[i].reunioes[aux].data[1]]) + horarioTercTurnoSoma[0][pop[i].reunioes[aux].data[1]];
    for (int i = 1; i < 4; i++){
        valorAux2 = (10 * horarioProfSoma[1][pop[i].reunioes[aux].data[1]]) + horarioTercTurnoSoma[i%2][pop[i].reunioes[aux].data[1]];
        if (valorAux2 < valorAux){
            melhorHorario = i;
            valorAux = valorAux2;
        }
    }

    pop[i].reunioes[aux].horario = melhorHorario;
}

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

//voy a empezar generando numeros aleatorios
//Me he quedado con la matriz las condiicones de contorno y la energia
//Lo que no se es como hacer lo de los 100 pasos de montecarlo y que es lo que tengo que haver despues

#define N 100//Dimension de la matriz

//PASO 0: Doy un T, genero la matriz inicial.
#define T 3.0
#define PMC 500

void GenerarMatriz (int s[][N+2])
{
int i,j;
for (i=1;i<=N;i++){
    for (j=1;j<=N;j++){
    //Genero un numero aleatorio en cada posicion, 0 o 1
    s[i][j]=rand() %2;
    //Transformo los valores obtenidos en 1 y -1
    s[i][j] = 2 *s[i][j]-1;

     }

}
}
//Comprobada en un fichero aparte que funciona
 
int main() {
    int s[N+2][N+2];
    int i, j;
    double suma;
    int E;
    double p;
    double deltaE;
    int numero1, numero2;
    int aleatorio;
    double epsilon;
    int n;

    // Semilla para la función rand() basada en el tiempo actual
    srand(time(NULL));
    

    //Genero mi matriz desordenada de espines mediante una funcion que llamo a continuacion segun la dimension
    // Llenar la matriz con valores aleatorios de 1 y -1
    GenerarMatriz(s);

    // Imprimir la matriz
    //printf("Matriz generada:\n");
    //for (i =1; i <= N; i++) {
       // for (j=1; j <=N; j++) {
       //     printf("%d ", s[i][j]);
        //}
        //printf("\n");
 //   }

    //Condiciones de contorno
    //¿No podría meter todo esto en un solo ciclo for?
    for (j = 0; j <= N+1; j++){

     s[0][j]=s[N][j];
     s[N+1][j]=s[1][j];

    }

     for (j = 0; j <= N+1; j++){

     s[j][0]=s[j][N];
     s[j][N+1]=s[j][1];
     
    }
    
    //Abro el fichero donde voy a escribir mis matrices
    FILE *fichs;
    fichs = fopen("ising_datos.txt", "w");
    n=0;
    while (n<(PMC*pow(N,2)))
     {
    
    //PASO 1: ELEGIR UN VALOR DE LA RED

    //Ahora tengo que elegir un valor de la red de manera aleatoria
    //Para ello voy a generar dos numeros aleatorios y los introduzco como coeficientes de la matriz

    numero1 = rand() % N; //Numero aleatorio entero entre [0,N-1]
    numero2 = rand() % N;
    numero1+=1; //Hago esta correccion del valor obtenido para que se quede en la submatriz central
    numero2+=1;


    //PASO 2: CALCULAR LA ENERGÍA Y CALCULAR P 
    //Calculo la energia
    deltaE=2*s[numero1][numero2]*(s[numero1 +1][numero2]+s[numero1-1][numero2]+s[numero1][numero2+1]+s[numero1][numero2-1]) ;
    
    //Hago un if para quedarme con el minimo de los dos valores
    if(1<=exp(-deltaE/T)){
        p=1;
       
    }
    else{
        p=exp(-deltaE/T);

    }
        //PASO 3: CALCULAR EPSILON ALEATORIO Y DECIDIR EL CAMBIO DE SIGNO DEL SPIN
        //Genero el epsilon en [0,1] (no entero)
        epsilon=((double) rand()) /RAND_MAX; 
        if(epsilon < p){
        s[numero1][numero2]=-s[numero1][numero2];
        //Hay que reajustar las condiciones de contorno
        if(numero1==1) s[N+1][numero2]=s[numero1][numero2];
        if(numero1==N) s[0][numero2]=s[numero1][numero2];
        if(numero2==1) s[numero1][N+1]=s[numero1][numero2];
        if(numero2==N) s[numero1][0]=s[numero1][numero2];
    }
    
    
    //Compruebo que se este haciendo bien
    //printf("%d\t%d",-deltaE,p);


//Escribo la matriz en el fichero, solo para cada paso de montecarlo
 if (n % (N*N) == 0)
{
for(i=1;i<=N;i++)
{
for(j=1;j<N;j++) // no incluyo el ultimo elemento de cada fila para que no lleve coma
{
fprintf(fichs, "%d,\t", s[i][j]);
}
fprintf(fichs, "%d\n", s[i][N]); // para escribir el N-esimo elemento de cada fila sin , en el fichero y saltar linea
}
fprintf(fichs, "\n"); // para separar los bloques de tiempo en el fichero
}

    n++;
    
    }
    //Cierro el fichero
    fclose(fichs);

    //Aqui se termina el while y hemos obtenido un paso de montecarlo, necesito un millon

    return 0;
}





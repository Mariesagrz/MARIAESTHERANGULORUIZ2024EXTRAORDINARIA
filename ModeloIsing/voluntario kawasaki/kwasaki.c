#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define N 16 // Dimensión de la matriz
#define T 1.0
#define PMC 100 

void GenerarMatriz(int s[][N+2]) {
    int i, j;
    for (i = 1; i <= N; i++) {
        for (j = 1; j <= N; j++) {
            s[i][j] = rand() % 2;
            s[i][j] = 2 * s[i][j] - 1;
        }
    }
}

void ActualizarCondicionesDeContorno(int s[][N+2]) {
    int j;
    for (j = 0; j <= N + 1; j++) {
        s[0][j] = s[N][j];
        s[N + 1][j] = s[1][j];
        s[j][0] = s[j][N];
        s[j][N + 1] = s[j][1];
    }
}

double CalcularMagnetizacionTotal(int s[][N+2]) {
    int i, j;
    double magnetizacion = 0.0;
    for (i = 1; i <= N; i++) {
        for (j = 1; j <= N; j++) {
            magnetizacion += s[i][j];
        }
    }
    return magnetizacion ;
}

double CalcularEnergia(int s[][N+2]) {
    int i, j;
    double energia = 0.0;
    for (i = 1; i <= N; i++) {
        for (j = 1; j <= N; j++) {
            energia -= s[i][j] * (s[i+1][j] + s[i-1][j] + s[i][j+1] + s[i][j-1]);
        }
    }
    return energia / 2.0; // Para evitar contar dos veces cada par de espines
}

//Para la magnetizacion tambien voy a hacer una funcion
void MagnetizacionPorDominios(int s[][N+2], double *magSuperior, double *magInferior) {
    int i, j;
    *magSuperior = 0.0;
    *magInferior = 0.0;

    for (i = 1; i <= N / 2; i++) {
        for (j = 1; j <= N; j++) {
            *magSuperior += s[i][j];
        }
    }
    for (i = N / 2 + 1; i <= N; i++) {
        for (j = 1; j <= N; j++) {
            *magInferior += s[i][j];
        }
    }
    *magSuperior /= (N * N / 2);
    *magInferior /= (N * N / 2);
}


//Para la densidad de particulas en el eje y hago un promedio de los espines de cada fila y luego hago un promedio de los promedios
double Densidad (int s[][N+2]) {
    int i, j;
    double densidad, PromDensidad;
    for (i = 1; i <= N; i++) {
        densidad=0.0;
        for (j = 1; j <= N; j++) {
            densidad+= s[i][j];
        }
        densidad=densidad/N;
        PromDensidad+= densidad;
    }
    PromDensidad=PromDensidad/N;
    return PromDensidad;
}




int main() {
    int s[N+2][N+2];
    int i, j, n;
    double deltaE, energiaAntes, energiaDespues, magSuperior, magInferior, densidad, magnetizacion;
    int x1, y1, x2, y2;
    double p, epsilon, EnergiaMedia, EnergiaMediaParticula;

    srand(time(NULL));
    GenerarMatriz(s);
    ActualizarCondicionesDeContorno(s);

    FILE *fichs;
    fichs = fopen("kawasaki_datos.txt", "w");
    n = 0;
    FILE *fichm;
    fichm = fopen("magnetizacion.txt", "w");
    FILE *fichMG;
    fichMG = fopen("magnetizacion_total.txt", "w");
    FILE *fichd;
    fichd = fopen("densidad.txt", "w");

    magnetizacion=0.0;

    while (n < (PMC * pow(N, 2))) {
        // Seleccionar un espín inicial aleatorio
        x1 = (rand() % N) + 1;
        y1 = (rand() % N) + 1;

        // Seleccionar un vecino aleatorio
        int direction = rand() % 4; // Generar un número aleatorio entre 0 y 3

        switch (direction) {
            case 0: // Vecino a la derecha
                x2 = x1;
                y2 = y1 + 1;
                break;
            case 1: // Vecino a la izquierda
                x2 = x1;
                y2 = y1 - 1;
                break;
            case 2: // Vecino arriba
                x2 = x1 - 1;
                y2 = y1;
                break;
            case 3: // Vecino abajo
                x2 = x1 + 1;
                y2 = y1;
                break;
        }

        // Manejar condiciones de contorno
        if (x2 == 0) x2 = N;
        if (x2 == N + 1) x2 = 1;
        if (y2 == 0) y2 = N;
        if (y2 == N + 1) y2 = 1;

        // Si los espines son diferentes, considerar el intercambio
        if (s[x1][y1] != s[x2][y2]) {
            // Calcular la energía antes del intercambio
            energiaAntes = CalcularEnergia(s);

            // Realizar el intercambio
            int temp = s[x1][y1];
            s[x1][y1] = s[x2][y2];
            s[x2][y2] = temp;

            // Actualizar condiciones de contorno
            ActualizarCondicionesDeContorno(s);

            // Calcular la energía después del intercambio
            energiaDespues = CalcularEnergia(s);

            // Calcular deltaE
            deltaE = energiaDespues - energiaAntes;

            // Calcular la probabilidad de aceptación
            if (deltaE <= 0) {
                p = 1;
            } else {
                p = exp(-deltaE / T);
            }

            // Generar un número aleatorio epsilon en [0, 1)
            epsilon = ((double) rand()) / RAND_MAX;

            // Decidir si aceptamos el intercambio
            if (epsilon >= p) {
                // Si no aceptamos, revertir el intercambio
                temp = s[x1][y1];
                s[x1][y1] = s[x2][y2];
                s[x2][y2] = temp;

                // Actualizar condiciones de contorno
                ActualizarCondicionesDeContorno(s);
                magnetizacion= CalcularMagnetizacionTotal(s);
            }
        }
        MagnetizacionPorDominios(s, &magSuperior, &magInferior);
        
         fprintf(fichm, "%d\t%lf\t%lf\n", n, magSuperior, magInferior);
         densidad=Densidad(s);
         fprintf(fichd, "%d\t%lf\n", n, densidad);
         

        
         // Vamos a ir calculando la energía media
        EnergiaMedia+=CalcularEnergia(s);

        // Escribir la matriz en el fichero cada N*N pasos de Monte Carlo
        if (n % (N * N) == 0) {
            for (i = 1; i <= N; i++) {
                for (j = 1; j < N; j++) {
                    fprintf(fichs, "%d,\t", s[i][j]);
                }
                fprintf(fichs, "%d\n", s[i][N]);
            }
            fprintf(fichs, "\n");
        }
        
        n++;
    }
    EnergiaMedia=EnergiaMedia/(PMC*pow(N, 2));
    EnergiaMediaParticula=EnergiaMedia/pow(N, 2);
    //magnetizacion = magnetizacion / (PMC * pow(N, 2));
    fprintf(fichMG, "%d\t%lf\n", n, magnetizacion);

 fclose(fichs);
 fclose(fichm);
 fclose(fichd);
 fclose(fichMG);
    return 0;
}




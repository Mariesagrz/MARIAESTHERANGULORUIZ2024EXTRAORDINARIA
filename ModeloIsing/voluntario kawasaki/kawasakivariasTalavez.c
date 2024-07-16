#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define N  48// Dimensión de la matriz
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

double Densidad(int s[][N+2]) {
    int i, j;
    double densidad = 0.0;
    double PromDensidad = 0.0;
    for (i = 1; i <= N; i++) {
        densidad = 0.0;
        for (j = 1; j <= N; j++) {
            densidad += s[i][j];
        }
        densidad = densidad / N;
        PromDensidad += densidad;
    }
    PromDensidad = PromDensidad / N;
    return PromDensidad;
}

int main() {
    int s[N+2][N+2];
    int i, j, n;
    double deltaE, energiaAntes, energiaDespues, magSuperior, magInferior, densidad;
    int x1, y1, x2, y2;
    double p, epsilon, EnergiaMedia, EnergiaMediaParticula;
    double T;

    srand(time(NULL));
    GenerarMatriz(s);
    ActualizarCondicionesDeContorno(s);

    FILE *fichs;
    fichs = fopen("energia_media.txt", "w"); // Almacena energía media por partícula
    FILE *fichm;
    fichm = fopen("magnetizacion.txt", "w"); // Almacena magnetización por dominios
    FILE *fichd;
    fichd = fopen("densidad.txt", "w"); // Almacena densidad media

    for (T = 1.0; T <4.0; T += 1.0) { // Iterar sobre diferentes temperaturas
        n = 0;
        EnergiaMedia = 0.0;
        double SumaDensidad = 0.0; // Suma de densidades para promedio
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
                    p = 1.0;
                } else {
                    p = exp(-deltaE / T);
                }

                epsilon = ((double) rand() / RAND_MAX);
                if (epsilon >= p) {
                    // Deshacer el intercambio si no se acepta
                    temp = s[x1][y1];
                    s[x1][y1] = s[x2][y2];
                    s[x2][y2] = temp;

                    // Actualizar condiciones de contorno
                    ActualizarCondicionesDeContorno(s);
                } else {
                    n++;
                    EnergiaMedia += energiaDespues;
                    SumaDensidad += Densidad(s); // Acumular densidad
                }
            }
        }
        EnergiaMediaParticula = EnergiaMedia / (PMC * pow(N, 2));
        MagnetizacionPorDominios(s, &magSuperior, &magInferior);
        densidad = SumaDensidad / (PMC * pow(N, 2)); // Promedio de densidad

        // Guardar resultados en los ficheros correspondientes
        fprintf(fichs, "%lf\t%lf\n", T, EnergiaMediaParticula);
        fprintf(fichm, "%lf\t%lf\t%lf\n", T, magSuperior, magInferior);
        fprintf(fichd, "%lf\t%lf\n", T, densidad);
    }

    fclose(fichs);
    fclose(fichm);
    fclose(fichd);

    return 0;
}

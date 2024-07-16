#include <stdio.h>
#include <math.h>
#include <stdlib.h>

// Defino todas las variables globales que hay en este problema

#define G 6.67e-11 // Constante gravitacional 
#define PI 3.1415
#define M_T 5.9736e24 // Masa de la Tierra en kg
#define M_L 0.07349e24 // Masa de la Luna en kg
#define d_TL 3.844e8 // Distancia Tierra-Luna en metros
#define w 2.6617e-6 // Velocidad angular de la Luna en rad/segundo
#define R_T 6.37816e6 // Radio de la Tierra en metros
#define R_L 1.7374e6 // Radio de la Luna en metros
#define TSIM 10 *24 *60*60 // numero de dias de la simulacion en segundos


#define Delta (G * M_T / pow(d_TL, 3))
#define mu (M_L / M_T)

int main() {
    int i, j, contador;
    double t, h, theta, v, m_cohete, pos[3][2], r, phi, p_r, p_phi, k[5][5];
    double T, V, H, Hprima;
    FILE *fichpos, *fichk, *fichH, *fichT, *fichV, *fichwpphi, *fichHprima;

    // Datos
    t = 10; // tiempo en días
    h = 0.8*60*5; // paso del tiempo en segundos
    m_cohete = 1; // masa del cohete en kg
    phi = 27; // latitud en grados
    theta = 0; // ángulo de lanzamiento en grados
    v = 11176; // velocidad de salida en metros por segundo

    // Convertir unidades
    v = v / d_TL;
    phi = phi * PI / 180.0; // radianes
    theta = theta * PI / 180.0; // radianes

    // Inicializar posiciones
    pos[0][0] = 0.0; // eje x Tierra
    pos[0][1] = 0.0; // eje y Tierra
    pos[1][0] = 1.0; // eje x Luna
    pos[1][1] = 0.0; // eje y Luna
    pos[2][1] = 0.0; // eje y cohete
    pos[2][0] = R_T / d_TL; // eje x cohete

    // Inicializar momentos y r
    r = R_T / d_TL;
    p_r = v;
    p_phi = 0.0;
    
    //printf("%lf\t%lf\t%lf\t%lf\n", r, phi, p_r, p_phi);
    //printf("\n");

    // Abrir archivos
    fichpos = fopen("posiciones.txt", "w");
    fichk = fopen("k.txt", "w");
    fichT = fopen("T.txt", "w");
    fichV = fopen("V.txt", "w");
    fichH = fopen("hamiltoniano.txt", "w");
    fichwpphi = fopen("wpphi.txt", "w"); 
    fichHprima = fopen("hamiltonianoprima.txt", "w");

    // Escribir posiciones iniciales
    fprintf(fichpos, "%lf, %lf\n", pos[0][0], pos[0][1]);
    fprintf(fichpos, "%lf, %lf\n", pos[1][0], pos[1][1]);
    fprintf(fichpos, "%lf, %lf\n", pos[2][0], pos[2][1]);
    fprintf(fichpos, "\n");

    // Algoritmo de Runge-Kutta
    t = 0.0;
    contador = 0;
    double rprima;
    while (t <= TSIM) {
        rprima = pow(1 + pow(r, 2) - 2 * r * cos(phi - w * t), 1.5);

        // Calcular k1
        k[1][1] = h * p_r;
        k[1][2] = h * p_phi / pow(r, 2);
        k[1][3] = h * (p_phi * p_phi / pow(r, 3) - Delta * (1.0 / pow(r, 2) + mu * (r - cos(phi - w * t)) / rprima));
        k[1][4] = h * (-Delta * mu * r * sin(phi - w * t) / rprima);

        // Calcular k2
        k[2][1] = h * (p_r + k[1][3] / 2.0);
        k[2][2] = h * (p_phi + k[1][4] / 2.0) / pow(r + k[1][1] / 2.0, 2);
        k[2][3] = h * ((p_phi + k[1][4] / 2.0) * (p_phi + k[1][4] / 2.0) / pow(r + k[1][1] / 2.0, 3) - Delta * (1.0 / pow(r + k[1][1] / 2.0, 2) + mu * ((r + k[1][1] / 2.0) - cos((phi + k[1][2] / 2.0) - w * (t + h / 2.0))) / pow(1 + pow((r + k[1][1] / 2.0), 2) - 2 * (r + k[1][1] / 2.0) * cos((phi + k[1][2] / 2.0) - w * (t + h / 2.0)), 1.5)));
        k[2][4] = h * (-Delta * mu * (r + k[1][1] / 2.0) * sin((phi + k[1][2] / 2.0) - w * (t + h / 2.0)) / pow(1 + pow((r + k[1][1] / 2.0), 2) - 2 * (r + k[1][1] / 2.0) * cos((phi + k[1][2] / 2.0) - w * (t + h / 2.0)), 1.5));

        // Calcular k3
        k[3][1] = h * (p_r + k[2][3] / 2.0);
        k[3][2] = h * (p_phi + k[2][4] / 2.0) / pow(r + k[2][1] / 2.0, 2);
        k[3][3] = h * ((p_phi + k[2][4] / 2.0) * (p_phi + k[2][4] / 2.0) / pow(r + k[2][1] / 2.0, 3) - Delta * (1.0 / pow(r + k[2][1] / 2.0, 2) + mu * ((r + k[2][1] / 2.0) - cos((phi + k[2][2] / 2.0) - w * (t + h / 2.0))) / pow(1 + pow((r + k[2][1] / 2.0), 2) - 2 * (r + k[2][1] / 2.0) * cos((phi + k[2][2] / 2.0) - w * (t + h / 2.0)), 1.5)));
        k[3][4] = h * (-Delta * mu * (r + k[2][1] / 2.0) * sin((phi + k[2][2] / 2.0) - w * (t + h / 2.0)) / pow(1 + pow((r + k[2][1] / 2.0), 2) - 2 * (r + k[2][1] / 2.0) * cos((phi + k[2][2] / 2.0) - w * (t + h / 2.0)), 1.5));

        // Calcular k4
        k[4][1] = h * (p_r + k[3][3]);
        k[4][2] = h * (p_phi + k[3][4]) / pow(r + k[3][1], 2);
        k[4][3] = h * ((p_phi + k[3][4]) * (p_phi + k[3][4]) / pow(r + k[3][1], 3) - Delta * (1.0 / pow(r + k[3][1], 2) + mu * ((r + k[3][1]) - cos((phi + k[3][2]) - w * (t + h))) / pow(1 + pow((r + k[3][1]), 2) - 2 * (r + k[3][1]) * cos((phi + k[3][2]) - w * (t + h)), 1.5)));
        k[4][4] = h * (-Delta * mu * (r + k[3][1]) * sin((phi + k[3][2]) - w * (t + h)) / pow(1 + pow((r + k[3][1]), 2) - 2 * (r + k[3][1]) * cos((phi + k[3][2]) - w * (t + h)), 1.5));

        // Actualizar valores
        r = r + (k[1][1] + 2.0 * k[2][1] + 2.0 * k[3][1] + k[4][1]) / 6.0;
        phi = phi + (k[1][2] + 2.0 * k[2][2] + 2.0 * k[3][2] + k[4][2]) / 6.0;
        p_r = p_r + (k[1][3] + 2.0 * k[2][3] + 2.0 * k[3][3] + k[4][3]) / 6.0;
        p_phi = p_phi + (k[1][4] + 2.0 * k[2][4] + 2.0 * k[3][4] + k[4][4]) / 6.0;

        //printf("%lf\t%lf\t%lf\t%lf\n", r, phi, p_r, p_phi);
        //printf("\n");

        // Escribir resultados
        if (contador % 1 == 0) {
            pos[0][0] = 0.0;
            pos[0][1] = 0.0;
            pos[1][0] = cos(w*t);
            pos[1][1] = sin(w*t);
            pos[2][0] = r * cos(phi);
            pos[2][1] = r * sin(phi);
            fprintf(fichpos, "%lf, %lf\n", pos[0][0], pos[0][1]);
            fprintf(fichpos, "%lf, %lf\n", pos[1][0], pos[1][1]);
            fprintf(fichpos, "%lf, %lf\n", pos[2][0], pos[2][1]);
            fprintf(fichpos, "\n");
        }
        contador++;

        // Calcular energía
        T = 0.5 * (p_r * p_r + p_phi * p_phi / (r * r));
        V = -Delta * (1.0 / r + mu * (r * cos(phi - w * t) - 1.0) / pow(1.0 + r * r - 2.0 * r * cos(phi - w * t), 0.5));
        H = T + V;
        Hprima = (H * pow(d_TL / 60.0, 2) / M_T);

        // Escribir energía en archivos
        fprintf(fichT, "%lf\t%lf\n", t, T);
        fprintf(fichV, "%lf\t%lf\n", t, V);
        fprintf(fichH, "%lf\t%lf\n", t, H);
        fprintf(fichHprima, "%lf\t%lf\n", t, Hprima);

        // Incrementar tiempo
        t += h;
    }

    // Cerrar archivos
    fclose(fichpos);
    fclose(fichk);
    fclose(fichT);
    fclose(fichV);
    fclose(fichH);
    fclose(fichwpphi);
    fclose(fichHprima);

    return 0;
}

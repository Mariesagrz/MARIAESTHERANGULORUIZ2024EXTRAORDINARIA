#include <stdio.h>
#include <math.h>
#include <stdlib.h>

// Defino todas las variables globales que hay en este problema

#define G 9.81 // Aceleración gravitacional
#define PI 3.1415
#define TSIM 10 // número de segundos de la simulación

// Parámetros del péndulo doble
#define m1 1.0
#define m2 1.0
#define l1 1.0
#define l2 1.0

void derivadas(double *derivs, double *y) {
    double theta1 = y[0];
    double theta2 = y[1];
    double p_theta1 = y[2];
    double p_theta2 = y[3];

    double delta = theta2 - theta1;

    double denom1 = (m1 + m2) * l1 * l1 - m2 * l1 * l1 * cos(delta) * cos(delta);
    double denom2 = (l2 / l1) * denom1;

    derivs[0] = (p_theta1 * l2 - p_theta2 * l1 * cos(delta)) / (l1 * l1 * l2 * denom1);
    derivs[1] = (p_theta2 * (m1 + m2) * l1 - p_theta1 * m2 * l2 * cos(delta)) / (m2 * l1 * l2 * l2 * denom2);
    derivs[2] = -((m1 + m2) * G * l1 * sin(theta1) + derivs[0] * derivs[1] * m2 * l1 * l2 * sin(delta)) / l1;
    derivs[3] = -m2 * G * l2 * sin(theta2) + derivs[0] * derivs[1] * m2 * l1 * l2 * sin(delta) / l2;
}

void runge_kutta(double *y, double h, int n) {
    double k1[4], k2[4], k3[4], k4[4], ytemp[4];
    double derivs[4];

    derivadas(derivs, y);
    for (int i = 0; i < 4; i++) {
        k1[i] = h * derivs[i];
        ytemp[i] = y[i] + k1[i] / 2.0;
    }

    derivadas(derivs, ytemp);
    for (int i = 0; i < 4; i++) {
        k2[i] = h * derivs[i];
        ytemp[i] = y[i] + k2[i] / 2.0;
    }

    derivadas(derivs, ytemp);
    for (int i = 0; i < 4; i++) {
        k3[i] = h * derivs[i];
        ytemp[i] = y[i] + k3[i];
    }

    derivadas(derivs, ytemp);
    for (int i = 0; i < 4; i++) {
        k4[i] = h * derivs[i];
    }

    for (int i = 0; i < 4; i++) {
        y[i] += (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0;
    }
}

int main() {
    double t, h;
    double y[4];
    FILE *fichpos;

    // Datos
    t = 0.0; // tiempo inicial en segundos
    h = 0.01; // paso del tiempo en segundos

    // Condiciones iniciales
    y[0] = PI / 2; // theta1
    y[1] = PI / 2; // theta2
    y[2] = 0.0; // p_theta1
    y[3] = 0.0; // p_theta2

    // Abrir archivo
    fichpos = fopen("pendulo_doble_posiciones.txt", "w");

    // Algoritmo de Runge-Kutta
    while (t <= TSIM) {
        runge_kutta(y, h, 4);
        t += h;

        // Calcular las coordenadas x e y de las masas
        double x1 = l1 * sin(y[0]);
        double y1 = -l1 * cos(y[0]);
        double x2 = x1 + l2 * sin(y[1]);
        double y2 = y1 - l2 * cos(y[1]);

        // Escribir resultados
        fprintf(fichpos, "%lf, %lf\n", x1, y1);
        fprintf(fichpos, "%lf, %lf\n",x2, y2);
        fprintf(fichpos, "\n");
    }

    // Cerrar archivo
    fclose(fichpos);

    return 0;
}

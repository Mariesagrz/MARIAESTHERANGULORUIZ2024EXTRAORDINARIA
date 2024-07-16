#include <stdio.h>
#include <math.h>

// Declaramos las constantes globales
#define G 9.81  // Aceleración debida a la gravedad
#define L1 1.0  // Longitud del primer péndulo
#define L2 1.0  // Longitud del segundo péndulo
#define M1 1.0  // Masa del primer péndulo
#define M2 1.0  // Masa del segundo péndulo
#define DT 0.01  // Paso de tiempo
#define TSIM 1000 // Número de pasos de tiempo

// Función que calcula las derivadas según las ecuaciones del péndulo doble dadas
void derivadas(double t, double y[], double dydt[]) {
    double phi = y[0];
    double phidot = y[1];
    double psi = y[2];
    double psidot = y[3];

    double delta = psi - phi;
    double den1 = (M1 + M2) * L1 - M2 * L1 * cos(delta) * cos(delta);
    double den2 = (L2 / L1) * den1;

    dydt[0] = phidot;
    dydt[1] = (M2 * L1 * phidot * phidot * sin(delta) * cos(delta) +
              M2 * G * sin(psi) * cos(delta) +
              M2 * L2 * psidot * psidot * sin(delta) -
              (M1 + M2) * G * sin(phi)) / den1;
    dydt[2] = psidot;
    dydt[3] = (-M2 * L2 * psidot * psidot * sin(delta) * cos(delta) +
              (M1 + M2) * G * sin(phi) * cos(delta) -
              (M1 + M2) * L1 * phidot * phidot * sin(delta) -
              (M1 + M2) * G * sin(psi)) / den2;
}

// Función que calcula las constantes de Runge Kutta Clásico
void rungeKutta(double y[], double t, double h, int n) {
    int i;
    double k1[4], k2[4], k3[4], k4[4], yt[4];

    derivadas(t, y, k1);
    for (i = 0; i < n; i++) yt[i] = y[i] + h * k1[i] / 2;
    derivadas(t + h / 2, yt, k2);
    for (i = 0; i < n; i++) yt[i] = y[i] + h * k2[i] / 2;
    derivadas(t + h / 2, yt, k3);
    for (i = 0; i < n; i++) yt[i] = y[i] + h * k3[i];
    derivadas(t + h, yt, k4);
    for (i = 0; i < n; i++) y[i] += h * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) / 6;
}

int main() {
    double t = 0;
    double y[4];  // Array para almacenar phi, phidot, psi, psidot
    double y_pert[4];  // Array para almacenar las condiciones perturbadas
    double d0 = 1e-8;  // Perturbación inicial
    double energia;
    FILE *f = fopen("lyapunov.txt", "w");

    // Condiciones iniciales
    y[0] = M_PI / 4;  // Ángulo inicial de phi
    y[1] = 0.5;         // Velocidad angular inicial de phi
    y[2] = M_PI / 4;  // Ángulo inicial de psi
    y[3] = sqrt(2 * 3.0); // Velocidad inicial ajustada para la energía deseada

    // Condiciones perturbadas
    for (int i = 0; i < 4; i++) y_pert[i] = y[i];
    y_pert[0] += d0;  // Perturbación en el ángulo phi

    // Simulación
    for (int i = 0; i < TSIM; i++) {
        // Calcular la distancia entre las trayectorias
        double dist = sqrt(pow(y_pert[0] - y[0], 2) + pow(y_pert[1] - y[1], 2) + 
                           pow(y_pert[2] - y[2], 2) + pow(y_pert[3] - y[3], 2));
        fprintf(f, "%lf, %lf\n", t, dist);

        // Actualizar los valores utilizando Runge-Kutta
        rungeKutta(y, t, DT, 4);
        rungeKutta(y_pert, t, DT, 4);
        t += DT;
    }

    fclose(f);
    return 0;
}

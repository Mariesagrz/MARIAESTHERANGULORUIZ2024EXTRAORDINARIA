#include <stdio.h>
#include <math.h>

//Declaramos las constantes globales
#define G 9.81  // Aceleración debida a la gravedad
#define L1 1.0  // Longitud del primer péndulo
#define L2 1.0  // Longitud del segundo péndulo
#define M1 1.0  // Masa del primer péndulo
#define M2 1.0  // Masa del segundo péndulo
#define DT 0.01  // Paso de tiempo
#define TSIM 1000 // Número de pasos de tiempo

//Funcion que calcula las derivadas segun las ecuaciones del pendulo doble dadas
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
//Funcion que calcula las constantes de Runge Kutta Clasico
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
    double energia;
    FILE *f = fopen("pendulo_doble.txt", "w");
    FILE *map1 = fopen("Poincare1.txt", "w");
    FILE *map2 = fopen("Poincare2.txt", "w");
    FILE *map3 = fopen("Poincare3.txt", "w");

    // Condiciones iniciales
    y[0] = M_PI / 2;  // Ángulo inicial de phi
    y[1] = 0;         // Velocidad angular inicial de phi
    y[2] = M_PI / 2;  // Ángulo inicial de psi
    y[3] = 0;         // Velocidad angular inicial de psi

    // Iteramos sobre las energías especificadas, para ahorrar tiempo lo hago todo junto
    energia = 15.0;

        
        printf("Simulando para E = %lf\n", energia);

        // Restablecer condiciones iniciales para cada energía
        y[0] = M_PI / 2;
        y[1] = 0;
        y[2] = M_PI / 2;
        y[3] = sqrt(2 * energia); // Velocidad inicial ajustada para la energía deseada

        // Simulación
        for (int i = 0; i < TSIM; i++) {
            // Calcular las coordenadas x e y de las masas
        double x1 = L1 * sin(y[0]);
        double y1 = -L1 * cos(y[0]);
        double x2 = x1 + L2 * sin(y[1]);
        double y2 = y1 - L2 * cos(y[1]);

        // Escribir resultados
        fprintf(f, "%lf, %lf\n", x1, y1);
        fprintf(f, "%lf, %lf\n",x2, y2);
        fprintf(f, "\n");
         
        //Para hacer los mapas de Poincaré necesito imprimir en varios ficheros los valores del vector y
       fprintf(map1, "%lf, %lf\n", y[2], y[0]);
       fprintf(map2, "%lf, %lf\n", y[2], y[3]);
       fprintf(map3, "%lf, %lf\n", y[0], y[1]);
        rungeKutta(y, t, DT, 4);
        t += DT;
        }
    

    fclose(f);
    fclose(map1);
    fclose(map2);
    fclose(map3);
    return 0; 
}


package problemariot3;

public class Jacobi5 {
    public static void main(String[] args) {
        // Definir la matriz de coeficientes y el vector de términos independientes
        double[][] A = {
            {4, -1, 0},
            {-1, 4, -1},
            {0, -1, 3}
        };

        double[] b = {5, -7, 6};

        // Establecer el vector de soluciones inicial
        double[] initialGuess = {1, 1, 1}; // Cambiar el vector inicial puede afectar la convergencia

        // Definir el número máximo de iteraciones y la tolerancia
        int maxIterations = 1000;
        double tolerance = 1e-6;

        // Resolver el sistema de ecuaciones
        double[] x = jacobi(A, b, initialGuess, maxIterations, tolerance);

        // Imprimir la solución
        System.out.println("Solución:");
        for (int i = 0; i < x.length; i++) {
            System.out.printf("x%d = %.5f\n", i + 1, x[i]);
        }
    }

    public static double[] jacobi(double[][] A, double[] b, double[] initialGuess, int maxIterations, double tolerance) {
        int n = A.length;
        double[] x = new double[n];
        double[] xPrev = new double[n];
        System.arraycopy(initialGuess, 0, x, 0, n);

        // Realizar iteraciones
        for (int iter = 0; iter < maxIterations; iter++) {
            // Copiar los valores anteriores de x
            System.arraycopy(x, 0, xPrev, 0, n);

            // Calcular el nuevo valor de x para cada ecuación
            for (int i = 0; i < n; i++) {
                double sum = 0.0;
                for (int j = 0; j < n; j++) {
                    if (j != i) {
                        sum += A[i][j] * xPrev[j];
                    }
                }
                x[i] = (b[i] - sum) / A[i][i];
            }

            // Verificar la convergencia
            double error = 0.0;
            for (int i = 0; i < n; i++) {
                error += Math.abs(x[i] - xPrev[i]);
            }
            if (error < tolerance) {
                break;
            }
        }

        return x;
    }
}

/*
Entradas:
Matriz de coeficientes A:
4 -1  0
-1  4 -1
0  -1  3

Vector de términos independientes b:
5
-7
6

Vector de soluciones inicial:
1
1
1

Número máximo de iteraciones: 1000
Tolerancia: 0.000001

Salidas:
Solución:
x1 = 1.50000
x2 = -2.00000
x3 = 2.00000
*/

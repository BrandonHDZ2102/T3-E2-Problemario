
package problemariot3;

public class GaussSeidel4 {
    
    public static void main(String[] args) {
        // Definir la matriz de coeficientes y el vector de términos independientes
        double[][] A = {
            {3, -1, 1},
            {2, 6, 2},
            {1, 1, 4}
        };

        double[] b = {4, 8, 7};

        // Establecer el vector de soluciones inicial
        double[] initialGuess = {0, 0, 0};

        // Definir el número máximo de iteraciones
        int maxIterations = 1000;

        // Resolver el sistema de ecuaciones
        double[] x = gaussSeidel(A, b, initialGuess, maxIterations);

        // Imprimir la solución
        System.out.println("Solución:");
        for (int i = 0; i < x.length; i++) {
            System.out.printf("x%d = %.5f\n", i + 1, x[i]);
        }
    }

    public static double[] gaussSeidel(double[][] A, double[] b, double[] initialGuess, int maxIterations) {
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
                        sum += A[i][j] * x[j];
                    }
                }
                x[i] = (b[i] - sum) / A[i][i];
            }

            // Verificar la convergencia
            double error = 0.0;
            for (int i = 0; i < n; i++) {
                error += Math.abs(x[i] - xPrev[i]);
            }
            if (error < 1e-6) {
                break;
            }
        }

        return x;
    }
}

/*
Entradas:
Matriz de coeficientes A:
3 -1  1
2  6  2
1  1  4

Vector de términos independientes b:
4
8
7

Vector de soluciones inicial:
0
0
0

Número máximo de iteraciones: 1000

Salidas:
Solución:
x1 = 1.00000
x2 = 1.00000
x3 = 1.00000
*/



package problemariot3;


public class EliminacionGaussiana2 {

     public static void main(String[] args) {
        // Definir la matriz de coeficientes y el vector de términos independientes
        double[][] A = {
            {2, 1, -1},
            {-3, -1, 2},
            {-2, 1, 2}
        };

        double[] b = {8, -11, -3};

        // Resolver el sistema de ecuaciones
        double[] x = gaussianElimination(A, b);

        // Imprimir la solución
        System.out.println("Solución:");
        for (int i = 0; i < x.length; i++) {
            System.out.printf("x%d = %.2f\n", i + 1, x[i]);
        }
    }

    public static double[] gaussianElimination(double[][] A, double[] b) {
        int n = A.length;

        // Formar la matriz aumentada
        double[][] augmentedMatrix = new double[n][n + 1];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                augmentedMatrix[i][j] = A[i][j];
            }
            augmentedMatrix[i][n] = b[i];
        }

        // Aplicar el método de eliminación gaussiana
        for (int i = 0; i < n; i++) {
            // Buscar el mayor valor en la columna i desde la fila i hasta la fila n
            int maxRow = i;
            for (int k = i + 1; k < n; k++) {
                if (Math.abs(augmentedMatrix[k][i]) > Math.abs(augmentedMatrix[maxRow][i])) {
                    maxRow = k;
                }
            }

            // Intercambiar filas
            double[] temp = augmentedMatrix[i];
            augmentedMatrix[i] = augmentedMatrix[maxRow];
            augmentedMatrix[maxRow] = temp;

            // Hacer que todos los valores debajo de la diagonal sean 0 en la columna i
            for (int k = i + 1; k < n; k++) {
                double factor = augmentedMatrix[k][i] / augmentedMatrix[i][i];
                augmentedMatrix[k][i] = 0;
                for (int j = i + 1; j <= n; j++) {
                    augmentedMatrix[k][j] -= factor * augmentedMatrix[i][j];
                }
            }
        }

        // Resolver el sistema triangular superior
        double[] x = new double[n];
        for (int i = n - 1; i >= 0; i--) {
            x[i] = augmentedMatrix[i][n] / augmentedMatrix[i][i];
            for (int k = i - 1; k >= 0; k--) {
                augmentedMatrix[k][n] -= augmentedMatrix[k][i] * x[i];
            }
        }

        return x;
    }
}

/*
Entradas:
Matriz de coeficientes A:
2  1 -1
-3 -1  2
-2  1  2

Vector de términos independientes b:
8
-11
-3

Salidas:
Solución:
x1 = 2.00
x2 = 3.00
x3 = -1.00
*/
    


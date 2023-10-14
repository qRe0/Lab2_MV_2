import java.util.Arrays;

public class LinearSystemSolver {

    public static void main(String[] args) {
        int N = 3; // Размерность системы
        double eps = 1e-6; // Порог сходимости
        int kmax = 10; // Максимальное число итераций
        double w = 1.0; // Параметр для метода релаксации

        // Инициализация массивов для векторов X и B
        double[] x = new double[N]; // Вектор приближенных решений
        double[] xNext = new double[N]; // Вектор для хранения следующей итерации при релаксации
        double[] b = new double[N]; // Вектор B

        // Инициализация матрицы A
        double[][] A = new double[N][N];
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                if (i == j) {
                    A[i][j] = 15 * Math.pow(i + 1, 2); // Диагональные элементы
                } else {
                    A[i][j] = -1.0 / ((i + 1) * (j + 1)); // Внедиагональные элементы
                }
            }
            x[i] = 0.0; // Начальное приближение
        }

        // Инициализация вектора b
        for (int i = 0; i < N; i++) {
            b[i] = 2.0 + i;
        }

        // Метод Якоби
        int iteration = 0;
        while (iteration < kmax) {
            for (int i = 0; i < N; i++) {
                double sum = 0.0;
                for (int j = 0; j < N; j++) {
                    if (j != i) {
                        sum += A[i][j] * x[j];
                    }
                }
                xNext[i] = (1.0 / A[i][i]) * (b[i] - sum);
            }

            double maxDiff = 0.0;
            for (int i = 0; i < N; i++) {
                double diff = Math.abs(xNext[i] - x[i]);
                if (diff > maxDiff) {
                    maxDiff = diff;
                }
                x[i] = xNext[i];
            }

            iteration++;

            if (maxDiff < eps) {
                System.out.println("Jacobi method converged after " + iteration + " iterations.");
                break;
            }
        }

        if (iteration >= kmax) {
            System.out.println("Jacobi method did not converge within " + kmax + " iterations.");
        }

        // Вывод приближенного решения
        System.out.println("Approximate solution (Jacobi method):");
        for (int i = 0; i < N; i++) {
            System.out.println("x[" + i + "] = " + x[i]);
        }

        // Проверка правильности решения и модификация матрицы для пунктов 2 и 3
        double[] exactSolution = new double[N];
        for (int i = 0; i < N; i++) {
            exactSolution[i] = 2.0 + i;
        }

        double normDifference = 0.0;
        double normExact = 0.0;
        for (int i = 0; i < N; i++) {
            double diff = Math.abs(exactSolution[i] - x[i]);
            normDifference += diff * diff;
            normExact += exactSolution[i] * exactSolution[i];
        }

        double relativeError = Math.sqrt(normDifference) / Math.sqrt(normExact);
        System.out.println("Cubic norm of the difference: " + Math.sqrt(normDifference));
        System.out.println("Relative error: " + relativeError);

        // Модификация матрицы (пункт 3)
        for (int i = 1; i < N; i++) {
            for (int j = 0; j < N; j++) {
                A[i][j] = Math.abs(A[i][j]);
            }
        }

        // Метод Гаусса-Зейделя
        iteration = 0;
        // Сброс начального приближения
        Arrays.fill(x, 0.0);

        while (iteration < kmax) {
            for (int i = 0; i < N; i++) {
                double sum1 = 0.0;
                double sum2 = 0.0;

                for (int j = 0; j < i; j++) {
                    sum1 += A[i][j] * x[j];
                }

                for (int j = i + 1; j < N; j++) {
                    sum2 += A[i][j] * x[j];
                }

                xNext[i] = (1.0 / A[i][i]) * (b[i] - sum1 - sum2);
            }

            double maxDiff = 0.0;
            for (int i = 0; i < N; i++) {
                double diff = Math.abs(xNext[i] - x[i]);
                if (diff > maxDiff) {
                    maxDiff = diff;
                }
                x[i] = xNext[i];
            }

            iteration++;

            if (maxDiff < eps) {
                System.out.println("Gauss-Seidel method converged after " + iteration + " iterations.");
                break;
            }
        }

        if (iteration >= kmax) {
            System.out.println("Gauss-Seidel method did not converge within " + kmax + " iterations.");
        }

        // Вывод приближенного решения
        System.out.println("Approximate solution (Gauss-Seidel method):");
        for (int i = 0; i < N; i++) {
            System.out.println("x[" + i + "] = " + x[i]);
        }

        // Метод релаксации
        iteration = 0;
        // Сброс начального приближения
        Arrays.fill(x, 0.0);

        while (iteration < kmax) {
            for (int i = 0; i < N; i++) {
                double sum1 = 0.0;
                double sum2 = 0.0;

                for (int j = 0; j < i; j++) {
                    sum1 += A[i][j] * x[j];
                }

                for (int j = i + 1; j < N; j++) {
                    sum2 += A[i][j] * x[j];
                }

                xNext[i] = (1 - w) * x[i] + (w / A[i][i]) * (b[i] - sum1 - sum2);
            }

            double maxDiff = 0.0;
            for (int i = 0; i < N; i++) {
                double diff = Math.abs(xNext[i] - x[i]);
                if (diff > maxDiff) {
                    maxDiff = diff;
                }
                x[i] = xNext[i];
            }

            iteration++;

            if (maxDiff < eps) {
                System.out.println("Relaxation method with w=" + w + " converged after " + iteration + " iterations.");
                break;
            }
        }

        if (iteration >= kmax) {
            System.out.println("Relaxation method with w=" + w + " did not converge within " + kmax + " iterations.");
        }

        // Вывод приближенного решения
        System.out.println("Approximate solution (Relaxation method with w=" + w + "):");
        for (int i = 0; i < N; i++) {
            System.out.println("x[" + i + "] = " + x[i]);
        }

//        //Тестовый вывод
//        // Вывод матрицы A
//        System.out.println("Matrix A:");
//        for (int i = 0; i < N; i++) {
//            for (int j = 0; j < N; j++) {
//                System.out.printf("%.2f\t", A[i][j]);
//            }
//            System.out.println();
//        }
//
//        // Вывод вектора B
//        System.out.println("Vector B:");
//        for (int i = 0; i < N; i++) {
//            System.out.printf("%.2f\t", b[i]);
//        }
//        System.out.println();
    }
}

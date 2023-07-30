#pragma once
#ifndef QR_Decomposition
#define QR_Decomposition

#include <iostream> 
#include <iomanip>

class SquareMatrix
{
    size_t n;                                                             // Размер матрицы
    double** A;                                                           // Начальная матрица

    void AllocMatrixMemory(double**& matrix)                              // Выделение памяти
    {
        matrix = new double* [n];
        for (size_t i = 0; i < n; i++)
        {
            matrix[i] = new double[n];
        }
    }
    void FreeMatrixMemory(double**& matrix)                               // Освободжение памяти
    {
        for (size_t i = 0; i < n; i++)
        {
            delete[] matrix[i];
        }
        delete[] matrix;
    }
    void IdentityMatrix(double**& matrix)                                 // Единичная матрица
    {
        for (size_t i = 0; i < n; i++)
        {
            for (size_t j = 0; j < n; j++)
            {
                i == j ? matrix[i][j] = 1 : matrix[i][j] = 0;
            }
        }
    }
    void MartixTranspose(double**& matrix)                                // Транспонирование
    {
        double a(0);
        for (size_t j = 0; j < n - 1; j++)
        {
            for (size_t i = j + 1; i < n; i++)
            {
                a = matrix[i][j];
                matrix[i][j] = matrix[j][i];
                matrix[j][i] = a;
            }
        }
    }
    void MatrixProd(double**& result, double** matrix1, double** matrix2) // Произведение матриц
    {
        double temp(0);
        for (size_t i = 0; i < n; i++)
        {
            for (size_t j = 0; j < n; j++)
            {
                for (size_t k = 0; k < n; k++)
                    temp += matrix1[k][j] * matrix2[i][k];
                result[i][j] = temp;
                temp = 0;
            }
        }
    }
    void MatrixPrint(double** matrix)                                     // Вывод матрицы
    {
        std::cout << "\n";
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < n; j++)
                std::cout << std::setw(15) << matrix[i][j] << "   ";
            std::cout << "\n";
        }
    }
    void RandomMatrix(double**& matrix)                                   // Заполнение матрицы произвольными значениями
    {
        for (size_t i = 0; i < n; i++)
        {
            for (size_t j = 0; j < n; j++)
            {
                matrix[i][j] = rand() % 10;
            }
        }
    }

public:
    SquareMatrix(size_t size)                                             // Конструктор
    {
        n = size;
        AllocMatrixMemory(A);
        RandomMatrix(A);
    }
    ~SquareMatrix()                                                       // Деструктор
    {
        FreeMatrixMemory(A);
    }

    void QR()                                                             // QR разложение
    {
        double s, c, temp;                                                // Синус, косинус

        double** T;                                                       // Создание матрицы вращений
        AllocMatrixMemory(T);
        IdentityMatrix(T);

        std::cout << "Initial matrix\n";                       // Вывод начальной матрицы
        MatrixPrint(A);

        double rho(0);
        for (size_t j = 0; j < n - 1; j++)                                // Проход по всем элементам под диагональю (Вращение Гивенса)
            for (size_t i = j + 1; i < n; i++)
            {
                rho = sqrt(pow(A[j][j], 2) + pow(A[i][j], 2));

                c = A[j][j] / rho;                                        // ШАГ1 - вычисление /cos/phi
                s = A[i][j] / rho;                                        //      - вычисление /sin/phi

                for (size_t k = j; k < n; k++)                            // ШАГ2 - вращение j-ой и i-ой строки
                {
                    temp = A[j][k];
                    A[j][k] = A[j][k] * c + A[i][k] * s;
                    A[i][k] = A[i][k] * c - temp * s;
                }

                for (size_t k = 0; k < n; k++)                            // ШАГ3 - формирование матрицы Т
                {
                    temp = T[j][k];
                    T[j][k] = T[j][k] * c + T[i][k] * s;
                    T[i][k] = T[i][k] * c - temp * s;
                }
            }

        std::cout << "\nMatrix R\n";                           // Вывод матрицы R
        MatrixPrint(A);

        MartixTranspose(T);                                               // Получение матрицы Q
        std::cout << "\nMatrix Q:\n";                          // Вывод матрицы Q
        MatrixPrint(T);

        std::cout << "\nMatrix Q*R:\n";                        // Вывод матрицы Q*R
        double** QR;
        AllocMatrixMemory(QR);
        MatrixProd(QR, A, T);
        MatrixPrint(QR);

        FreeMatrixMemory(T);
        FreeMatrixMemory(QR);
    }
};

#endif
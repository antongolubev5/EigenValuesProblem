
// Lab2slau.cpp: ���������� ����� ����� ��� ����������� ����������.
//

#include "stdafx.h"
#include <fstream>
#include <iostream>
#include <malloc.h>
#include <iomanip>
#include <cmath>
using namespace std;

double** readfromfile(string filename, int &n)
{
	char buff[500];
	double element;
	int count = 0;
	ifstream fin(filename); // ������� ���� ��� ������		

	while (!fin.eof())//������������ ����� �����
	{
		fin.getline(buff, 500);
		count++;
	}
	//�������� ������ ��� ������
	int i;
	double **extMatrix;
	extMatrix = new double*[count];
	for (i = 0; i<count; i++)
		//a[i], a[i] �������� � ��������� ���� double
		extMatrix[i] = new double[count + 1];

	fin.close(); // ��������� ����
	fin.open(filename); // ������� ���� ��� ������
	int elcnt = 0;
	while (!fin.eof())
	{
		fin >> element; // ������� ��������� �������
		div_t ij = div(elcnt, count + 1);
		extMatrix[ij.quot][ij.rem] = element;
		elcnt++;
	}
	n = count;
	fin.close(); // ��������� ����
	return extMatrix;
}

double** readfromfilesquare(string filename, int &n)
{
	char buff[500];
	double element;
	int count = 0;
	ifstream fin(filename); // ������� ���� ��� ������		

	while (!fin.eof())//������������ ����� �����
	{
		fin.getline(buff, 500);
		count++;
	}
	//�������� ������ ��� ������
	int i;
	double **Matrix;
	Matrix = new double*[count];
	for (i = 0; i<count; i++)
		//a[i], a[i] �������� � ��������� ���� double
		Matrix[i] = new double[count];

	fin.close(); // ��������� ����
	fin.open(filename); // ������� ���� ��� ������
	int elcnt = 0;

	while (!fin.eof())
	{
		fin >> element; // ������� ��������� �������
		div_t ij = div(elcnt, count);
		Matrix[ij.quot][ij.rem] = element;
		elcnt++;
	}
	n = count;
	fin.close(); // ��������� ����
	return Matrix;
}

void matrix_destroyer(double** ary, int n)
{
	if (ary != nullptr)
	{
		for (int i = 0; i < n; i++) {
			delete[] ary[i];
		}
		delete[] ary;
		ary = nullptr;
	}
}

void vector_destroyer(double* vec, int n)
{
	if (vec != nullptr)
	{

		delete[] vec;
		vec = nullptr;
	}
}

double** readfromscreen(int &k)
{
	int i, j, raz;
	cout << "������� ����������� �������: ";
	cin >> raz;
	k = raz;
	double **extMatrix;
	extMatrix = new double*[raz];

	for (i = 0; i<raz; i++)
		extMatrix[i] = new double[raz + 1];

	for (i = 0; i < raz; i++)
	{
		for (j = 0; j < raz + 1; j++)
		{
			cin >> extMatrix[i][j];
		}
	}
	return extMatrix;
}

double** readfromscreensquare(int &k)
{
	int i, j, raz;
	cout << "������� ����������� ������� A: ";
	cin >> raz;
	k = raz;
	double **extMatrix;
	extMatrix = new double*[raz];

	for (i = 0; i<raz; i++)
		//a[i], a[i] �������� � ��������� ���� double
		extMatrix[i] = new double[raz];

	for (i = 0; i < raz; i++)
	{
		for (j = 0; j < raz; j++)
		{
			cin >> extMatrix[i][j];
		}
	}
	return extMatrix;
}

void printMatrix(double** extMatrix, int k, int m, bool extended)
{
	cout << endl;
	for (int i = 0; i < k; i++)
	{
		for (int j = 0; j < m; j++)
		{
			char simb = ' ';
			//if ((j == k - 1) && extended) { simb = '='; }
			//else { simb = ' '; }
			//cout <<  setprecision(2) << fixed<< extMatrix[i][j] << simb;
			cout  << extMatrix[i][j] << simb;
		}
		cout << endl;
	}

}

void printVector(double* Vector, int k)
{
	cout << endl;
	for (int i = 0; i < k; i++)
	{
		cout << Vector[i] /*<< " " << i+1 */ << endl;
	}

}

double** multiplyMatrix(double** Matrix1, double** Matrix2, int n)
{
	double **result;
	result = new double*[n];
	double s = 0;
	for (int ii = 0; ii<n; ii++)
		result[ii] = new double[n];

	for (int i = 0; i < n; i++)
	{
		for (int l = 0; l < n; l++)
		{
			s = 0;

			for (int j = 0; j < n; j++)
			{
				s += Matrix1[i][j] * Matrix2[j][l];
			}
			result[i][l] = s;
		}
	}
	return result;
}

double* multiplyMatrixVector(double** Matrix, double* Vector, int n)
{
	double *result;
	result = new double[n];
	double s = 0;

	for (int i = 0; i < n; i++)
	{
		s = 0;

		for (int j = 0; j < n; j++)
		{
			s += Matrix[i][j] * Vector[j];
		}
		result[i] = s;

	}
	return result;
}

double** multiplyMatrixNumber(double** Matrix, double Number, int n)
{
	double **result;
	result = new double*[n];
	for (int ii = 0; ii<n; ii++)
		result[ii] = new double[n];

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			result[i][j] = Matrix[i][j] * Number;
		}
	}
	return result;
}

double* multiplyVectorNumber(double* Vector, double Number, int n)
{
	double *result;
	result = new double[n];

	for (int i = 0; i < n; i++)
	{
		result[i] = Vector[i] * Number;
	}
	return result;
}

double* substractVector(double* Vector1, double* Vector2, int n)
{
	double *result;
	result = new double[n];

	for (int i = 0; i < n; i++)
	{
		result[i] = Vector1[i] - Vector2[i];
	}

	return result;
}

double** substractMatrix(double** Matrix1, double** Matrix2, int n)
{
	double **result;
	result = new double*[n];

	for (int ii = 0; ii<n; ii++)
		result[ii] = new double[n];

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			result[i][j] = Matrix1[i][j] - Matrix2[i][j];
		}
	}
	return result;
}

double** sumMatrix(double** Matrix1, double** Matrix2, int n)
{
	double **result;
	result = new double*[n];

	for (int ii = 0; ii<n; ii++)
		result[ii] = new double[n];

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			result[i][j] = Matrix1[i][j] + Matrix2[i][j];
		}
	}
	return result;
}

double** transposeMatrix(double** Matrix, int n)
{
	double **result;
	result = new double*[n];

	for (int ii = 0; ii<n; ii++)
		result[ii] = new double[n];

	for (int i = 0; i < n; i++)
	{
		for (int j = i; j < n; j++)
		{
			result[i][j] = Matrix[j][i];
			result[j][i] = Matrix[i][j];
		}
	}

	return result;
}

void unitMatrix(double** Matrix, int n)
{
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			if (i == j) { Matrix[i][j] = 1; }
			else  Matrix[i][j] = 0;
		}
	}
}

double discrepancy(double** extMatrix, double* solution, int n)
{
	double *b1;
	b1 = new double[n];
	double result = 0;

	b1 = multiplyMatrixVector(extMatrix, solution, n);

	for (int i = 0; i < n - 1; i++)
	{
		result += pow((b1[i] - extMatrix[i][n]), 2);
	}
	result = sqrt(result);
	return result;
}

double normVectorUnit(double* Vector, int n)
{
	double max = 0;
	for (int i = 0; i < n; i++)
	{
		max += abs(Vector[i]);
	}
	return max;
}

double normVectorInfinity(double* Vector, int n)
{
	double max = 0;
	for (int i = 0; i < n; i++)
	{
		if (abs(Vector[i]) > max) max = abs(Vector[i]);
	}
	return max;
}

double normVectorEuclid(double* Vector, int n)
{
	double result = 0;
	for (int i = 0; i < n; i++)
	{
		result += pow(Vector[i], 2);
	}
	result = pow(result, 0.5);
	return result;
}

double vectormultiplyvectorscalar(double* Vector1, double* Vector2, int n)
{
	double result = 0;
	for (int i = 0; i < n; i++)
	{
		result += Vector1[i]*Vector2[i];
	}

	return result;
}

double normMatrixUnit(double** Matrix, int n)//���� �������
{
	double max = 0;
	double *Vector;
	Vector = new double[n];
	for (int i = 0; i < n; i++)
	{
		Vector[i] = 0;
	}

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			Vector[j] += Matrix[i][j];
		}
	}
	max = normVectorInfinity(Vector, n);
	return max;
}

double normMatrixInfinity(double** Matrix, int n)//���� ������
{
	double max = 0;
	double *Vector;
	Vector = new double[n];
	for (int i = 0; i < n; i++)
	{
		Vector[i] = 0;
	}

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			Vector[i] += Matrix[i][j];
		}
	}
	max = normVectorInfinity(Vector, n);
	return max;
}

double* gaussmethod(double** extMatrix, int n)
{
	double *solution;
	solution = new double[n];
	int imax;
	double maxvalue = 0;

	for (int cnt = 0; cnt < n; cnt++)
	{
		solution[cnt] = 0;
	}
	//double det = determinant(extMatrix, n);
	//if (abs(det) < 1e-30)
	//{
	//cout << "������������ ����� 0. �� ���������� ������������� �������." << endl;
	//solution = nullptr;
	//}
	//else
	//{

	for (int i = 0; i < n - 1; i++)//���� �� �������, ������� ���������� �� �����������
	{
		//����� ���� �������� �� i-�� �������
		maxvalue = 0;
		for (int il = i; il < n; il++)
		{
			if (maxvalue < abs(extMatrix[il][i]))
			{
				maxvalue = abs(extMatrix[il][i]);
				imax = il;
			}
		}

		/*if (maxvalue < 1e-10)
		{
			cout << "�� ���������� ������������� �������." << endl;
			return nullptr;
		}*/

		if (imax != i)
		{
			double* buf = extMatrix[imax];
			extMatrix[imax] = extMatrix[i];
			extMatrix[i] = buf;
		}

		//extMatrix[i][n] = extMatrix[i][n] / extMatrix[i][i];
		double aii = extMatrix[i][i];

		/*if (abs(aii) < 1e-10)
		{
			cout << "�� ���������� ������������� �������. ��������� ������ ������������������� ������� - �������" << endl;
			return nullptr;
		}*/

		for (int j = i; j <= n; j++)//���� �� ��������� ��������, ������� ���������� �� �����������  �� i+1???
		{
			extMatrix[i][j] = extMatrix[i][j] / aii;
		}

		for (int ii = i + 1; ii < n; ii++)//��������� �� ���������� ����� i-�� ������
		{
			double a_ii_i = extMatrix[ii][i];
			for (int jj = i; jj <= n; jj++)
			{
				extMatrix[ii][jj] -= a_ii_i * extMatrix[i][jj];
			}
		}
	}
	//��������� ������ ������
	double	 aii = extMatrix[n - 1][n - 1];
	/*if (abs(aii) < 1e-10)
	{
		cout << "�� ���������� ������������� �������. ��������� ������ ������������������� ������� - �������" << endl;
		return nullptr;
	}*/
	for (int j = n - 1; j <= n; j++)//���� �� ��������� ��������, ������� ���������� �� �����������  �� i+1???
	{
		extMatrix[n - 1][j] = extMatrix[n - 1][j] / aii;
	}
	//printMatrix(extMatrix, n, n + 1, true);
	//�������� ���

	double sum = 0;
	for (int i = n - 1; i >= 0; i--)
	{
		sum = 0;
		for (int j = i + 1; j < n; j++) //��������� ��� ����� ������� ����������  ���������� �� ������������ ������� ������
		{
			sum += solution[j] * extMatrix[i][j];
		}
		solution[i] = extMatrix[i][n] - sum;//�������� �� ������ ����� 
	}

	//printMatrix(extMatrix, n);//������ ������������������� (��� ��������)
	return solution;
}

double** reverseMatrix(double** Matrix, int n)
{
	double **extMatrix;
	extMatrix = new double*[n];
	double **result;
	result = new double*[n];
	double *ordinary;
	ordinary = new double[n];

	for (int ii = 0; ii<n; ii++)
		result[ii] = new double[n];
	for (int ii = 0; ii<n; ii++)
		extMatrix[ii] = new double[n + 1];

	for (int j = 0; j < n; j++)//��������� ������ ������� extMatrix ������
	{
		extMatrix[j][n] = 0;
	}

	for (int i = 0; i < n; i++)
	{
		for (int ii = 0; ii < n; ii++)//��������� �������� �� ������� � ����������� �������
		{
			for (int jj = 0; jj < n; jj++)
			{
				extMatrix[ii][jj] = Matrix[ii][jj];
			}
		}
		for (int jj = 0; jj < n; jj++)//��������� ������ ������� extMatrix ������
		{
			extMatrix[jj][n] = 0;
		}

		extMatrix[i][n] = 1;

		ordinary = gaussmethod(extMatrix, n);
		if (ordinary == nullptr)
		{
			cout << "�� ���������� �������� �������." << endl;
			return nullptr;
		}
		for (int j = 0; j < n; j++)
		{
			result[j][i] = ordinary[j];
		}

	}

	return result;
}

double** HessenbergForm(double** Matrix, int n)
{
	double alpha, beta;
	double** tkl;
	tkl = new double*[n];
	double** tkltranspose=nullptr;
	//tkltranspose = new double*[n];
	double** MatrixNew;
	MatrixNew = new double*[n];
	double** result=nullptr;
	//result = new double*[n];
	for (int ii = 0; ii < n; ii++) {
		//result[ii] = new double[n];
		tkl[ii] = new double[n];
		MatrixNew[ii] = new double[n];
		//tkltranspose[ii] = new double[n];
	}

	//�������������� tkl � �������������� ��� ���������
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (i != j) tkl[i][j] = 0;
			else tkl[i][j] = 1;
		}
	}


	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			MatrixNew[i][j] = Matrix[i][j];
			
		}
	}

	//for (int i = 0; i < n-1; i++)
		for (int i = 2; i < n ; i++)//str
	{
		for (int j = 0; j < i-1; j++)
		{//�������� ������� i � ������� ������� tkl
			//i=l, j =k-1, ����� �������� l=i, k = j+1 
			alpha = MatrixNew[j + 1][j] / pow((MatrixNew[j + 1][j] * MatrixNew[j + 1][j] + MatrixNew[i][j] * MatrixNew[i][j]), 1.0 / 2.0);
			beta = MatrixNew[i][j] / pow((MatrixNew[j + 1][j] * MatrixNew[j + 1][j] + MatrixNew[i][j] * MatrixNew[i][j]), 1.0 / 2.0);

			tkl[j + 1][j + 1] = tkl[i][i] =alpha; //tkk = tll
			tkl[i][i] = tkl[i][i] = alpha; //tkl=tlk
			tkl[j + 1][i] = beta; //tkl
			tkl[i][j+1] = -beta; //tlk
			//printMatrix(tkl, n, n, true);

            //�������� � ����� �� tkl � ������ �� tkl �����������������
			tkltranspose = transposeMatrix(tkl, n);
			result = multiplyMatrix(tkl, MatrixNew, n);
			MatrixNew = multiplyMatrix(result, tkltranspose, n);
			//����������� ������
			if (tkltranspose != nullptr) { matrix_destroyer(tkltranspose, n); };
			if (result != nullptr) { matrix_destroyer(result, n); };
			/*cout << " A result:" << endl;
			printMatrix(MatrixNew, n, n, true);*/

			for (int i = 0; i < n; i++) {
				for (int j = 0; j < n; j++) {
					if (i != j) tkl[i][j] = 0;
					else tkl[i][j] = 1;
				}
			}
		}

	}

	return MatrixNew;
}

double** qrdecomposition(double** Matrix, int n, bool flag)
{
	double c, s;
	double *solution;
	solution = new double[n];
	double **t;
	t = new double*[n];
	double **tij;
	tij = new double*[n];
	double **r;
	r = new double*[n];
	double *b;
	b = new double[n];
	double **eiimax;
	eiimax = new double*[n];
	double **multiplyqr;
	multiplyqr = new double*[n];

	double **oldmemory;

	for (int ii = 0; ii < n; ii++)
	{
		t[ii] = new double[n];
		tij[ii] = new double[n];
		r[ii] = new double[n];
		eiimax[ii] = new double[n];
		multiplyqr[ii] = new double[n];
	}


	for (int i = 0; i < n; i++)//��������� �������� �� ������� ������� � ������� r
	{
		for (int j = 0; j < n; j++)
		{
			r[i][j] = Matrix[i][j];
		}
	}

	unitMatrix(t, n);

	double maxvalue = 0;
	int imax;

		//������ �������� ����� �� i � j
		for (int i = 0; i < n - 1; i++)
		{

			maxvalue = 0;
			for (int il = i; il < n; il++)
			{
				if (maxvalue < abs(Matrix[il][i]))
				{
					maxvalue = abs(Matrix[il][i]);
					imax = il;
				}
			}

			/*if (maxvalue < 1e-10)
			{
				cout << "�� ���������� ������������� �������." << endl;
				return nullptr;
			}*/

				if (imax != i)
				{//������������ �����, � i-�� ������ ���������� imax-������ � ��������
					double* buf = r[imax];
					r[imax] = r[i];
					r[i] = buf;
					//��������� ������� t �� ��������������� ������� eiimax

					unitMatrix(eiimax, n);//�������������� ������� t ��� ��������� ��� ������������ ������������

					eiimax[i][i] = 0;//��������� ��������������� �������
					eiimax[imax][imax] = 0;
					eiimax[i][imax] = 1;
					eiimax[imax][i] = 1;
					//���� ����������� ����� ��� ����� � qr, �� ����� ���������� ������
					oldmemory = t;
					t = multiplyMatrix(eiimax, t, n);
					matrix_destroyer(oldmemory, n);
				}


			for (int j = i + 1; j < n; j++)
			{
				c = r[i][i] / (sqrt(r[i][i] * r[i][i] + r[j][i] * r[j][i]));
				s = r[j][i] / (sqrt(r[i][i] * r[i][i] + r[j][i] * r[j][i]));

				unitMatrix(tij, n);//tij ���������������� ��������� ��������

				tij[i][i] = c;
				tij[i][j] = s;
				tij[j][i] = -s;
				tij[j][j] = c;
				/*cout << "������� t" <<i<<j<< endl;
				printMatrix(tij, n, n, false);*/
				oldmemory = t;
				t = multiplyMatrix(tij, t, n);//��������� ������� t �� ������ ����� ��������
				matrix_destroyer(oldmemory, n);
				oldmemory = r;
				r = multiplyMatrix(tij, r, n);//��������� ������� r �� ������ ����� ��������
				matrix_destroyer(oldmemory, n);
				/*cout << "������� r" << i << j << endl;
				printMatrix(r, n, n, false);*/
			}
		}
		//����� �������� ����� �� i � j
	

	//��������� ������ ������ �� ��� ������� ������������
	double	 aii = r[n - 1][n - 1];
	/*if (abs(aii) < 1e-10)
	{
		cout << "�� ���������� ������������� �������. ��������� ������ ������������������� ������� - �������" << endl;
		printMatrix(r, n, n, false);
		return nullptr;
	}*/

	double** q = transposeMatrix(t, n);
	if (flag)
	{
		cout << endl;
		cout << "������� Q:" << endl;
		cout << endl;
		printMatrix(q, n, n, false);
		cout << endl;
		cout << "������� R:" << endl;
		cout << endl;
		printMatrix(r, n, n, false);
		cout << endl;
		cout << "������������ ������ Q � R:" << endl;
		cout << endl;
		multiplyqr = multiplyMatrix(q, r, n);
		printMatrix(multiplyqr, n, n, false);
	}

	return q;
}

double** qreigenvalues(double** Matrix, int n, double eps)
{
	double** result = nullptr, **matrixtokill;
	result = new double*[n];
	double** q = nullptr;
	double** qtranspose = nullptr;
	bool convergence = false;
	int kol = 0, bigcnt=0,ibig=-1,jbig=-1;
	for (int ii = 0; ii < n; ii++)
	{
		result[ii] = new double[n];
	}

	for (int i = 0; i < n; i++)//��������� �������� �� ������� ������� � result
	{
		for (int j = 0; j < n; j++)
		{
			result[i][j] = Matrix[i][j];
		}
	}

	while (!convergence)
	{

	convergence = true;
	matrixtokill = q;
	q = qrdecomposition(result, n, false);
	if (matrixtokill != nullptr) { matrix_destroyer(matrixtokill, n); };
	matrixtokill = qtranspose;
	qtranspose = transposeMatrix(q, n);
	if (matrixtokill != nullptr) { matrix_destroyer(matrixtokill, n); };
	matrixtokill = result;
	result = multiplyMatrix(qtranspose, result, n);
	 matrix_destroyer(matrixtokill, n);
	result = multiplyMatrix(result, q, n);
	//printMatrix(result, n, n, false);
	bigcnt = 0;
		for (int i = 1; i < n && convergence; i++)//���� �� ������� ������������
			{
			for (int j = 0; j < i && convergence; j++)
			{
				if (abs(result[i][j]) > eps) {
					convergence = false; bigcnt++; ibig = i; jbig = j;
				}
			}

		}
		if (!convergence&&bigcnt == 1 && kol > 100)
		{//2 ����� ���������� �� ������ ��
		// http://pmpu.ru/vf4/algebra2/charpoly/qralgor

		 			//�������������. ��-� ��� ����. �����
			double a=1, b, c,d,lambi,lambj;
			b = -result[ibig][ibig] - result[jbig][jbig];
			c = result[ibig][ibig] * result[jbig][jbig] - result[ibig][jbig] * result[jbig][ibig];
			d = b*b - 4*a*c;
			if (d >= 0) {
				if (b > 0) { lambi = (-b + sqrt(d)) / 2; lambj = (-b - sqrt(d)) / 2; }
				else { lambj = (-b + sqrt(d)) / 2; lambi = (-b - sqrt(d)) / 2; }
				//�������� �������
				result[ibig][ibig] = lambi; result[jbig][jbig] = lambj;
				result[ibig][jbig] = 0.0; 
				return result;
			}
		}
		kol++;
    }

	cout << "���������� �������� = " << kol << endl;
	return result;
}

double** qreigenvaluesshift(double** Matrix, int n, double eps)
{
	double** result = nullptr, **matrixtokill;
	result = new double*[n];
	double** q = nullptr;
	double** qtranspose = nullptr;
	bool convergence = false;
	int kol = 0;
	double sigma = 1;
	for (int ii = 0; ii < n; ii++)
	{
		result[ii] = new double[n];
	}

	for (int i = 0; i < n; i++)//��������� �������� �� ������� ������� � result
	{
		for (int j = 0; j < n; j++)
		{
			result[i][j] = Matrix[i][j];
		}
	}
	int nullablestring = n - 1;
	while (nullablestring>=0)
	{
		sigma = result[nullablestring][nullablestring];

		for (int i = 0; i < n; i++)
		{
			result[i][i] -= sigma;
		}

		convergence = true;
		matrixtokill = q;
		q = qrdecomposition(result, n, false);
		if (matrixtokill != nullptr) { matrix_destroyer(matrixtokill, n); };
		matrixtokill = qtranspose;
		qtranspose = transposeMatrix(q, n);
		if (matrixtokill != nullptr) { matrix_destroyer(matrixtokill, n); };
		matrixtokill = result;
		result = multiplyMatrix(qtranspose, result, n);
		matrix_destroyer(matrixtokill, n);
		result = multiplyMatrix(result, q, n);
		//printMatrix(result, n, n, false);
		
			for (int j = 0; j < n&&j<=nullablestring && convergence; j++)
			{
				if (abs(result[nullablestring][j]) > eps) {
					convergence = false;
				}
			}
			if (convergence) 
			{ nullablestring--; 

			cout << "������=" << nullablestring + 1 << endl;
			printMatrix(result, n, n, false);
			}
		
		//for (int i = 1; i < n && convergence; i++)//���� �� ������� ������������
		//{
		//	for (int j = 0; j < i && convergence; j++)
		//	{
		//		if (abs(result[i][j]) > eps) {
		//			convergence = false;
		//		}
		//	}

		//}
		kol++;

		for (int i = 0; i < n; i++)
		{
			result[i][i] += sigma;
		}


	}
	/*
	convergence = false;
	while (!convergence)
	{

		convergence = true;
		matrixtokill = q;
		q = qrdecomposition(result, n, false);
		if (matrixtokill != nullptr) { matrix_destroyer(matrixtokill, n); };
		matrixtokill = qtranspose;
		qtranspose = transposeMatrix(q, n);
		if (matrixtokill != nullptr) { matrix_destroyer(matrixtokill, n); };
		matrixtokill = result;
		result = multiplyMatrix(qtranspose, result, n);
		matrix_destroyer(matrixtokill, n);
		result = multiplyMatrix(result, q, n);
		//printMatrix(result, n, n, false);

		for (int i = 1; i < n && convergence; i++)//���� �� ������� ������������
		{
			for (int j = 0; j < i && convergence; j++)
			{
				if (abs(result[i][j]) > eps) {
					convergence = false;
				}
			}

		}
		//kol++;
	}


	*/




	cout << "���������� �������� = " << kol << endl;
	return result;
}

double** qreigenvaluesshiftminor(double** Matrix, int n, double eps)
{
	double** result = nullptr, **matrixtokill;
	result = new double*[n];
	double** q = nullptr;
	double** qtranspose = nullptr;
	bool convergence = false;
	int kol = 0;
	double sigma = 1;
	for (int ii = 0; ii < n; ii++)
	{
		result[ii] = new double[n];
	}

	for (int i = 0; i < n; i++)//��������� �������� �� ������� ������� � result
	{
		for (int j = 0; j < n; j++)
		{
			result[i][j] = Matrix[i][j];
		}
	}
	int nullablestring = n - 1;
	while (nullablestring >= 0)
	{
		sigma = result[nullablestring][nullablestring];

		for (int i = 0; i < n; i++)
		{
			result[i][i] -= sigma;
		}

		convergence = true;
		matrixtokill = q;
		q = qrdecomposition(result, n, false);
		if (matrixtokill != nullptr) { matrix_destroyer(matrixtokill, n); };
		matrixtokill = qtranspose;
		qtranspose = transposeMatrix(q, n);
		if (matrixtokill != nullptr) { matrix_destroyer(matrixtokill, n); };
		matrixtokill = result;
		result = multiplyMatrix(qtranspose, result, n);
		matrix_destroyer(matrixtokill, n);
		result = multiplyMatrix(result, q, n);
		//printMatrix(result, n, n, false);

		for (int j = 0; j < n && j <= nullablestring && convergence; j++)
		{
			if (abs(result[nullablestring][j]) > eps) {
				convergence = false;
			}
		}
		if (convergence)
		{
			nullablestring--;
			n--;
			//printMatrix(result, n, n, false);
			cout << "����������� ����� L" << nullablestring + 2 << " = " << sigma << endl;
		}

		
		kol++;

		for (int i = 0; i < n; i++)
		{
			result[i][i] += sigma;
		}


	}

	cout << endl;
	cout << "���������� �������� = " << kol << endl;
	return result;
}

void printEigenValues(double** Matrix, int n)
{
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			if (i == j) cout << "����������� ����� L" << i + 1 << " = " << Matrix[i][j] << endl;
		}
	}
}

double* methodReverseIterations(double** Matrix, double eigenvalue, int n, double eps)
{
	double norm, normsubstract;
	double **extMatrix;
	extMatrix = new double*[n];
	double **boofMatrix;
	boofMatrix = new double*[n];
	double **unMatrix;
	unMatrix = new double*[n];
	double **SquareMatrix;
	SquareMatrix = new double*[n];
	double *x;
	x = new double[n];
	double *substract;
	substract = new double[n];

	double *y;//������� �������
	y = new double[n];
	bool flag = false;

	for (int ii = 0; ii < n; ii++)
	{
		extMatrix[ii] = new double[n + 1];//����������� ������� ��� ������� ���������
		SquareMatrix[ii] = new double[n];
		unMatrix[ii] = new double[n];
		boofMatrix[ii] = new double[n];
	}

	for (int i = 0; i < n; i++)//SquareMatrix
	{
		for (int j = 0; j< n; j++)
		{
			SquareMatrix[i][j] = 0;
		}
	}

	for (int i = 0; i < n; i++)//����������� x ��������� ����������� �0
	{
    	x[i] = 0;
	}
	x[0] = -1;
	//�������� ���������� ������� A-Li*E
	unitMatrix(unMatrix, n);

	/*cout << endl;
	cout << "��������� �������: " << endl;
	printMatrix(unMatrix, n, n, false);*/

	unMatrix = multiplyMatrixNumber(unMatrix, eigenvalue, n);

	/*cout << endl;
	cout << "��������� ������� ���������� �� Li: " << endl;
	printMatrix(unMatrix, n, n, false);*/

	SquareMatrix = substractMatrix(Matrix, unMatrix, n);

	/*cout << endl;
	cout << "� - ��������� ������� ���������� �� Li: " << endl;
	printMatrix(SquareMatrix, n, n, false);*/

	for (int ii = 0; ii < n; ii++)//��������� �������� �� ������� � ����������� �������
	{
		for (int jj = 0; jj < n; jj++)
		{
			extMatrix[ii][jj] = SquareMatrix[ii][jj];
		}
	}

	for (int jj = 0; jj < n; jj++)//��������� � ������ ������� ����������� ������� ��������� ����������� x0
	{
		extMatrix[jj][n] = x[jj];
	}

	/*cout << endl;
	cout << "���������� �������: " << endl;
	printMatrix(extMatrix, n, n+1, true);*/

	while (!flag)
	{
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n + 1; j++)
			{
				boofMatrix[i][j] = extMatrix[i][j];
			}

		}

		y = gaussmethod(boofMatrix, n);//������ �������
		//��������� y
		norm = normVectorEuclid(y, n);
		y = multiplyVectorNumber(y, 1.0 / norm, n);

		/*cout << "������������� �: " << endl;
		printVector(y, n);*/

		for (int i = 0; i < n; i++)
		{
			substract[i] = y[i] - x[i];
		}

		//multiplyMatrixVector(SquareMatrix, y, n);

		//normsubstract = normVectorEuclid(substract, n);
		normsubstract = normVectorEuclid(multiplyMatrixVector(SquareMatrix, y, n), n);
		flag = normsubstract < eps;

		for (int i = 0; i < n; i++)
		{
			x[i] = y[i];
		}

		for (int jj = 0; jj < n; jj++)//��������� � ������ ������� ����������� ������� ������� ������ �
		{
			extMatrix[jj][n] = x[jj];
		}
	}

	return x;
}

double* methodRayleigh(double** Matrix, double eigenvalue, int n, double eps)
{
	int kol = 0;
	double relay = eigenvalue;
	double norm, normsubstract;
	double **extMatrix;
	extMatrix = new double*[n];
	double **boofMatrix;
	boofMatrix = new double*[n];
	double **unMatrix;
	unMatrix = new double*[n];
	double **SquareMatrix;
	SquareMatrix = new double*[n];
	double *x;
	x = new double[n];
	double *substract;
	substract = new double[n];

	double *y;//������� �������
	y = new double[n];
	bool flag = false;

	for (int ii = 0; ii < n; ii++)
	{
		extMatrix[ii] = new double[n + 1];//����������� ������� ��� ������� ���������
		SquareMatrix[ii] = new double[n];
		unMatrix[ii] = new double[n];
		boofMatrix[ii] = new double[n];
	}

	for (int i = 0; i < n; i++)//SquareMatrix
	{
		for (int j = 0; j< n; j++)
		{
			SquareMatrix[i][j] = 0;
		}
	}

	for (int i = 0; i < n; i++)//����������� x ��������� ����������� �0
	{
		x[i] = 0;
	}
	x[0] = -1;
	//�������� ���������� ������� A-Li*E
	unitMatrix(unMatrix, n);

	/*cout << endl;
	cout << "��������� �������: " << endl;
	printMatrix(unMatrix, n, n, false);*/

	unMatrix = multiplyMatrixNumber(unMatrix, relay, n);

	/*cout << endl;
	cout << "��������� ������� ���������� �� Li: " << endl;
	printMatrix(unMatrix, n, n, false);*/

	SquareMatrix = substractMatrix(Matrix, unMatrix, n);

	/*cout << endl;
	cout << "� - ��������� ������� ���������� �� Li: " << endl;
	printMatrix(SquareMatrix, n, n, false);*/

	for (int ii = 0; ii < n; ii++)//��������� �������� �� ������� � ����������� �������
	{
		for (int jj = 0; jj < n; jj++)
		{
			extMatrix[ii][jj] = SquareMatrix[ii][jj];
		}
	}

	for (int jj = 0; jj < n; jj++)//��������� � ������ ������� ����������� ������� ��������� ����������� x0
	{
		extMatrix[jj][n] = x[jj];
	}

	/*cout << endl;
	cout << "���������� �������: " << endl;
	printMatrix(extMatrix, n, n+1, true);*/

	while (!flag)
	{
		
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n + 1; j++)
			{
				boofMatrix[i][j] = extMatrix[i][j];
			}

		}

		y = gaussmethod(boofMatrix, n);//������ �������
		
		norm = normVectorEuclid(y, n);//��������� y
		y = multiplyVectorNumber(y, 1.0 / norm, n);

		/*cout << "������������� �: " << endl;
		printVector(y, n);*/

		for (int i = 0; i < n; i++)
		{
			substract[i] = y[i] - x[i];
		}

		//multiplyMatrixVector(SquareMatrix, y, n);
		if (kol>3 || flag)
		{
			double axx = vectormultiplyvectorscalar(multiplyMatrixVector(Matrix, y, n), y, n);//(Ax,x)
			double xx = vectormultiplyvectorscalar(y, y, n);//(x,x)
			relay = axx / xx;
		}
		//normsubstract = normVectorEuclid(substract, n);
		normsubstract = normVectorEuclid(multiplyMatrixVector(SquareMatrix, y, n), n);
		flag = normsubstract < eps;

		for (int i = 0; i < n; i++)
		{
			x[i] = y[i];
		}

			unitMatrix(unMatrix, n);
			unMatrix = multiplyMatrixNumber(unMatrix, relay, n);
			SquareMatrix = substractMatrix(Matrix, unMatrix, n);
		
		for (int ii = 0; ii < n; ii++)//��������� �������� �� ������� � ����������� �������
		{
			for (int jj = 0; jj < n; jj++)
			{
				extMatrix[ii][jj] = SquareMatrix[ii][jj];
			}
		}

		for (int jj = 0; jj < n; jj++)//��������� � ������ ������� ����������� ������� ������� ������ �
		{
			extMatrix[jj][n] = x[jj];
		}
		kol++;
	}

	cout << "���������� �������� = " << kol << endl;
	cout << "C���������� ����� = " <<  relay << endl;

	return x;
}

int main()
{
	setlocale(LC_ALL, "rus");
	int k = 4;
	double eps = 1e-4;
	double *solution;
	solution = new double[k];
	double *eigenvaluemas;//������ ����������� �����
	eigenvaluemas = new double[k];

	double** test = readfromfilesquare("test.txt", k);
	cout << "������� ��������� �: " << endl;
	printMatrix(test, k, k, true);
	cout << endl;
	cout << "�������� ������� ��������� � � ����� �����������: " << endl;
	test = HessenbergForm(test, k);
	printMatrix(test, k, k, true);
	cout << endl;
	cout << "���� ����������� �������� � ������� QR-������: " << endl;
	cout << endl;
	cout << "�������� eps = " << eps << endl;
	cout << endl;
	//qrdecomposition(test, k, false);
	test = qreigenvalues(test, k, eps);
	printMatrix(test, k, k, false);
	cout << endl; 
	printEigenValues(test, k);
	cout << endl;

	//��������� ������ ����������� ����� ��� ������������ ���������� ����������� �������� ������� �������� ��������
	for (int i = 0; i < k; i++)
	{
		for (int j = 0; j < k; j++)
		{
			if (i == j) eigenvaluemas[i] = test[i][j];
		}
	}

	cout << "���� ����������� �������� � ������� QR-������ c� ��������: " << endl;
	double** test2 = readfromfilesquare("test.txt", k);
	cout << endl;
	cout << "�������� eps = " << eps << endl;
	cout << endl;
	test2 = qreigenvaluesshiftminor(test2, k, eps);
	cout << endl;

	double** test3 = readfromfilesquare("test.txt", k);
	cout << "���� ����������� ������� � ������� ������ �������� ��������: " << endl;
	cout << endl;
	cout << "������� ��������� �: " << endl;
	printMatrix(test3, k, k, false);
	cout << endl;
	
	//���� ����������� ������� ������� �������� ��������
	for (int i = 0; i < k; i++)
	{
		solution = methodReverseIterations(test3, eigenvaluemas[i], k, eps);
		cout << "����������� ������, ��������������� ������������ ����� L" << i << " = " << eigenvaluemas[i] << ": " << endl;
		printVector(solution, k);
		cout << endl;
	}
	//��������, ��� ������ ������������� �����������
	cout << "��������, ��� ������� ������������� ����������� (Ax=Lx): " << endl;
	cout << endl;
	double* res1 = multiplyMatrixVector(test3, solution, k);
	double* res2 = multiplyVectorNumber(solution, eigenvaluemas[3], k);
	cout << "Ax= " << endl;
	printVector(res1, k);
	cout << endl;

	cout << "Lx= " << endl;
	printVector(res2, k);
	cout << endl;
	double normdiscrep = normVectorEuclid(substractVector(res1, res2, k), k);
	cout << "����� ������� ������� ����� " << normdiscrep << endl;
	cout << endl;

	cout << "���� ����������� ����� � ����������� ������� ������� ��������� � ������� ��������� �����: " << endl;
	cout << endl;
	cout << "������� ��������� �: " << endl;
	printMatrix(test3, k, k, false);
	cout << endl;

	solution = methodRayleigh(test3, eigenvaluemas[0]+0.05, k, eps);
	cout << endl;
	cout << "����������� ������: " << endl;
	printVector(solution, k);
	cout << endl;
	//����� �������

	double** testnew = readfromfilesquare("testnew.txt", k);
	cout << "������� ��������� �: " << endl;
	printMatrix(testnew, k, k, true);
	cout << endl;
	cout << "�������� ������� ��������� � � ����� �����������: " << endl;
	testnew = HessenbergForm(testnew, k);
	//�����������, ������ �������� �������������, ��� ��������������� ���������� � ����� �����������
	printMatrix(testnew, k, k, true);
	cout << endl;
	cout << "���� ����������� �������� � ������� QR-������: " << endl;
	cout << endl;
	cout << "�������� eps = " << eps << endl;
	cout << endl;
	//qrdecomposition(test, k, false);
	testnew = qreigenvalues(testnew, k, eps);
	printMatrix(testnew, k, k, false);
	cout << endl;
	printEigenValues(testnew, k);
	cout << endl;

	system("pause");
	return 0;
}



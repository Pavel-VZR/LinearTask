#include <iostream>
#include <fstream>
#include <ctime>
#include <Eigen/Dense>

void ShowMatrix(const Eigen::MatrixXd& Matrix, std::string = "");
void ShowVector(const Eigen::VectorXd& Vector, std::string = "");
void ReadMatrix(Eigen::MatrixXd&, std::fstream&, std::string);
void ReadVector(Eigen::VectorXd&, std::fstream&, std::string);
void WriteVector(Eigen::VectorXd&, std::fstream&, std::string);
void WriteMatrix(Eigen::MatrixXd&, std::fstream&, std::string);
void LinearSolve(Eigen::VectorXd&, const Eigen::MatrixXd&, const Eigen::VectorXd&);
void GMRES(const Eigen::MatrixXd&, const Eigen::VectorXd&, Eigen::VectorXd&);

double RightF(double, double);
double K(double, double, double);
double ExactSolution(double, double);
double Error(const Eigen::VectorXd&, const Eigen::VectorXd&);
void SystemBuilder_mod_one(Eigen::MatrixXd&, Eigen::VectorXd&, double, int);

Eigen::VectorXd borderS, borderW, borderN, borderE;

int main()
{
	unsigned int StartTime = clock();

	//Размер области, число узлов, шаг
	double L = 1.;
	int Num = 20;
	double h = L / (Num - 1);
	std::cout << "h: " << h << "\n";
	std::cout << "Num: " << Num << "\n";

	borderS.resize(Num);
	borderN.resize(Num);
	borderW.resize(Num - 2);
	borderE.resize(Num - 2);

	//Граничные условия
	double x, y;
	x = h;
	for (int index = 0; index < Num - 2; index++)
	{
		borderW(index) = 1.;
		borderE(index) = sin(x) + 1.;

		x += h;
	}

	y = 0.;
	for (int index = 0; index < Num; index++)
	{
		borderS(index) = sin(y) + 1.;
		borderN(index) = 1.;

		y += h;
	}

	//Матрицы коэффициетов системы уравнений
	Eigen::MatrixXd LHS = Eigen::MatrixXd::Zero((Num - 2) * (Num - 2), (Num - 2) * (Num - 2));
	Eigen::VectorXd RHS = Eigen::VectorXd::Zero((Num - 2) * (Num - 2));

	//Значения на предыдущем шаге
	Eigen::VectorXd XPrev;
	//Точное решение
	Eigen::VectorXd XExact((Num - 2) * (Num - 2));
	//Текущее решение 
	Eigen::VectorXd X = Eigen::VectorXd::Constant((Num - 2) * (Num - 2), 1.);

	x = 0., y = 0.;
	for (int i = 0, index = 0; i < Num - 2; i++)
	{
		x += h;
		y = h;
		for (int j = 0; j < Num - 2; j++, index++)
		{
			XExact(index) = ExactSolution(x, y);
			y += h;
		}
	}

	//Построение системы
	SystemBuilder_mod_one(LHS, RHS, h, Num);
	
	//Параметры сходимости
	double eps = .1, res = 0, resPrev = 0;

	//Влияние решения с предыдущего шага
	double alpha = .8;
	do
	{
		XPrev = X;

		GMRES(LHS, RHS, X);

		res = Error(X, XExact);

		//Условия на сходимость
		if (!abs(res - resPrev))
			break;
		if (res < eps)
			break;

		resPrev = res;
		
		std::cout << "res: " << res << "\n";

		X = (1. - alpha) * XPrev + alpha * X;

	} while (true);


	//Построение матрицы для вывода
	Eigen::MatrixXd Answer(Num, Num);
	Answer.row(0) = borderN;
	Answer.row(Answer.rows() - 1) = borderS;
	Answer.block(1, 0, Answer.cols() - 2, 1) = borderW;
	Answer.block(1, Answer.rows() - 1, Answer.cols() - 2, 1) = borderE;
	for (int i = 1, index = 0; i < Answer.rows() - 1; i++)
		for (int j = 1; j < Answer.cols() - 1; j++, index++)
			Answer(i, j) = X(index);

	y = 0.;
	std::fstream Stream;
	Stream.open("Data\\Answer.txt", std::ios::out);
	for (int i = 0, index = 0; i < Num; i++)
	{
		x = 0.;
		for (int j = 0; j < Num; j++, index++)
		{
			Stream << x << " " << y << " " << Answer(i, j) << "\n";
			x += h;
		}
		y += h;
	}
	Stream.close();

	x = 0.; y = 0.;
	Stream.open("Data\\Error30.txt", std::ios::out);
	Eigen::MatrixXd Er((Num - 2), (Num - 2));
	for (int i = 0, index = 0; i < (Num - 2); i++)
	{
		y += h;
		x = h;
		for (int j = 0; j < Num - 2; j++, index++)
		{
			Stream << x << " " << y << " " << abs(X(index) - XExact(index)) / XExact(index) << "\n";
			x += h;
		}
	}
	Stream.close();

	unsigned int EndTime = clock();
	std::cout << "Program running time: " << (EndTime - StartTime) / 1000.;
}

void ShowMatrix(const Eigen::MatrixXd& Matrix, std::string text)
{
	std::cout << text;
	for (int i = 0; i < Matrix.rows(); i++)
	{
		for (int j = 0; j < Matrix.cols(); j++)
		{
			std::cout << Matrix(i, j) << " ";
		}
		std::cout << "\n";
	}
	std::cout << "\n";
}

void ShowVector(const Eigen::VectorXd& Vector, std::string text)
{
	std::cout << text;
	for (int i = 0; i < Vector.size(); i++)
	{
		std::cout << Vector(i) << " ";
	}
	std::cout << "\n";
}

void ReadMatrix(Eigen::MatrixXd& Matrix, std::fstream& RS, std::string FileName)
{
	int NumRow, NumCol;
	RS.open(FileName);

	if (!RS.is_open())
	{
		std::cout << "File open error!\n";
		exit(1);
	}

	RS >> NumRow >> NumCol;

	Matrix.resize(NumRow, NumCol);
	for (int i = 0; i < NumRow; i++)
	{
		for (int j = 0; j < NumCol; j++)
		{
			RS >> Matrix(i, j);
		}
	}
	RS.close();
}

void ReadVector(Eigen::VectorXd& Vector, std::fstream& RS, std::string FileName)
{
	int length;
	RS.open(FileName);

	if (!RS.is_open())
	{
		std::cout << "File open error!\n";
		exit(1);
	}

	RS >> length;

	Vector.resize(length);
	for (int i = 0; i < length; i++)
	{
		RS >> Vector(i);
	}
	RS.close();
}

void WriteVector(Eigen::VectorXd& v, std::fstream& stream, std::string FileName)
{
	stream.open(FileName, std::ios::out | std::ios::trunc);
	if (!stream.is_open())
	{
		std::cout << "File open error!\n";
		exit(1);
	}

	for (int i = 0; i < v.size(); i++)
	{
		stream << v(i) << "\n";
	}
	stream.close();
}

void WriteMatrix(Eigen::MatrixXd& Matrix, std::fstream& stream, std::string FileName)
{
	stream.open(FileName, std::ios::out | std::ios::trunc);
	if (!stream.is_open())
	{
		std::cout << "File open error!\n";
		exit(1);
	}

	for (int i = 0; i < Matrix.rows(); i++)
	{
		for (int j = 0; j < Matrix.cols(); j++)
			stream << Matrix(i, j) << " ";
		stream << "\n";
	}
	stream.close();
}

void LinearSolve(Eigen::VectorXd& x, const Eigen::MatrixXd& Left, const Eigen::VectorXd& Right)
{
	x.resize(Right.size());
	double temp;
	for (int i = Right.size() - 1; i >= 0; i--)
	{
		temp = 0.;
		for (int j = Right.size() - 1; j > i; j--)
			temp += Left(i, j) * x(j);

		x(i) = (Right(i) - temp) / Left(i, i);
	}
}

void GMRES(const Eigen::MatrixXd& A, const Eigen::VectorXd& b, Eigen::VectorXd& x)
{
	//Инициализация вектора начальных значений x_0
	Eigen::VectorXd x0;
	x0 = Eigen::VectorXd::Zero(b.size());

	//r0 = b - A * x0;
	Eigen::VectorXd r0(b - A * x0);

	//v=r0 / ||r0||
	Eigen::VectorXd v(r0 / r0.norm());

	double residual = 1.;
	double tol = 1.e-2;
	size_t k = 1;

	//Матрица ортогональных базисных векторов Крыловского подпространства 
	Eigen::MatrixXd V(b.size(), 1);

	//Матрица Хессенберга (размерность k+1 на k)
	Eigen::MatrixXd H;

	Eigen::MatrixXd QT = Eigen::MatrixXd::Identity(2, 2);
	Eigen::MatrixXd J, Htemp, R;

	Eigen::VectorXd e1, cc, w;
	V.col(0) = v;

	while (residual > tol)
	{
		H.conservativeResize(k + 1, k);
		H.row(k) = Eigen::VectorXd::Zero(k);

		w = A * v;

		for (int j = 0; j < k; j++)
		{
			H(j, k - 1) = V.col(j).transpose() * w;
			w = w - H(j, k - 1) * V.col(j);
		}

		H(k, k - 1) = w.norm();
		v = w / H(k, k - 1);
		V.conservativeResize(V.rows(), k + 1);
		V.col(k) = v;

		if (k == 1)
		{
			Htemp = H;
		}
		else
		{
			QT.conservativeResize(k + 1, k + 1);

			QT.col(k) = Eigen::VectorXd::Zero(k + 1);
			QT.row(k) = Eigen::VectorXd::Zero(k + 1);
			QT(k, k) = 1.;
			Htemp = QT * H;
		}
		J = Eigen::MatrixXd::Identity(k + 1, k + 1);

		J(k - 1, k - 1) = Htemp(k - 1, k - 1) / sqrt(pow(Htemp(k - 1, k - 1), 2) + pow(Htemp(k, k - 1), 2));
		J(k - 1, k) = Htemp(k, k - 1) / sqrt(pow(Htemp(k - 1, k - 1), 2) + pow(Htemp(k, k - 1), 2));
		J(k, k - 1) = -J(k - 1, k);
		J(k, k) = J(k - 1, k - 1);

		QT = J * QT;

		residual = fabs(QT(QT.rows() - 1, 0) * r0.norm());

		k++;
	}

	std::cout << "GMRES iteration converged in " << k - 1 << " steps\n\n";

	R = QT * H;
	R.conservativeResize(R.rows() - 1, R.cols());

	e1 = r0.norm() * Eigen::VectorXd::Unit(k, 0);
	cc = QT * e1;
	cc.conservativeResize(k - 1);

	Eigen::VectorXd y;
	LinearSolve(y, R, cc);

	V.conservativeResize(V.rows(), V.cols() - 1);

	//x.resize(x0.size());
	x = x0 + V * y;
}

void SystemBuilder_mod_one(Eigen::MatrixXd& left, Eigen::VectorXd& right, double h, int Num)
{
	double x = 0., y = 0;
	//std::cout << right << "\n";
	for (int i = 0, index = 0; i < Num - 2; i++)
	{
		x += h;
		y = h;
		for (int j = 0; j < Num - 2; j++, index++)
		{
			right(index) = h * RightF(x, y);
			y += h;
		}
	}
	//std::cout << right << "\n";

	double A, B, D, E;
	double revStep = 1. / h;
	double halfStep = h / 2.;
	x = h, y = h;
	//Левый верхний угол
	A = revStep * K(x - halfStep, y, 1.);
	B = revStep * K(x + halfStep, y, 1.);
	D = revStep * K(x, y - halfStep, 1.);
	E = revStep * K(x, y + halfStep, 1.);
	left(0, 0) = -(A + B + D + E);	//C
	left(0, 1) = E;				//E
	left(0, Num - 2) = B;		//B
	right(0) += -A * borderS(1) - D * borderW(0);

	for (int i = 1; i < Num - 3; i++)
	{
		y += h;
		//Первая строка
		A = revStep * K(x - halfStep, y, 1.);
		B = revStep * K(x + halfStep, y, 1.);
		D = revStep * K(x, y - halfStep, 1.);
		E = revStep * K(x, y + halfStep, 1.);
		left(i, i) = -(A + B + D + E);	//C
		left(i, i - 1) = D;				//D
		left(i, i + 1) = E;				//E
		left(i, i + Num - 2) = B;		//B
		right(i) += -A * borderS(i + 1);
	}

	y += h;
	//Правый верхний угол
	A = revStep * K(x - halfStep, y, 1.);
	B = revStep * K(x + halfStep, y, 1.);
	D = revStep * K(x, y - halfStep, 1.);
	E = revStep * K(x, y + halfStep, 1.);
	left(Num - 3, Num - 3) = -(A + B + D + E);				//C
	left(Num - 3, Num - 4) = D;								//D
	left(Num - 3, Num - 3 + Num - 2) = B;					//B
	right(Num - 3) += -A * borderS(Num - 2) - E * borderE(0);

	//Цикл по строкам (ячейки в центре)
	for (int i = 1; i < Num - 3; i++)
	{
		x += h; y = h;
		//Первый элемент в каждой ячейке
		A = revStep * K(x - halfStep, y, 1.);
		B = revStep * K(x + halfStep, y, 1.);
		D = revStep * K(x, y - halfStep, 1.);
		E = revStep * K(x, y + halfStep, 1.);
		left((Num - 2) * i, (Num - 2) * i) = -(A + B + D + E);				//C
		left((Num - 2) * i, (Num - 2) * i - (Num - 2)) = A;					//A
		left((Num - 2) * i, (Num - 2) * i + (Num - 2)) = B;					//B
		left((Num - 2) * i, (Num - 2) * i + 1) = E;							//E
		right((Num - 2) * i) += -D * borderW(i);

		//Центральные элементы ячейки
		for (int j = (Num - 2) * i + 1; j < (Num - 2) * i + Num - 3; j++)
		{
			y += h;
			A = revStep * K(x - halfStep, y, 1.);
			B = revStep * K(x + halfStep, y, 1.);
			D = revStep * K(x, y - halfStep, 1.);
			E = revStep * K(x, y + halfStep, 1.);
			left(j, j) = -(A + B + D + E);				//C
			left(j, j - (Num - 2)) = A;					//A
			left(j, j + (Num - 2)) = B;					//B
			left(j, j - 1) = D;							//D
			left(j, j + 1) = E;							//E
		}

		y += h;
		//Последний элемент в каждой ячейке
		A = revStep * K(x - halfStep, y, 1.);
		B = revStep * K(x + halfStep, y, 1.);
		D = revStep * K(x, y - halfStep, 1.);
		E = revStep * K(x, y + halfStep, 1.);
		left((Num - 2) * i + (Num - 3), (Num - 2) * i + (Num - 3)) = -(A + B + D + E);				//C
		left((Num - 2) * i + (Num - 3), (Num - 2) * i + (Num - 3) - (Num - 2)) = A;					//A
		left((Num - 2) * i + (Num - 3), (Num - 2) * i + (Num - 3) + (Num - 2)) = B;					//B
		left((Num - 2) * i + (Num - 3), (Num - 2) * i + (Num - 3) - 1) = D;							//D
		right((Num - 2) * i + (Num - 3)) += -E * borderE(i);
	}

	x += h; y = h;
	//Левый нижний угол
	A = revStep * K(x - halfStep, y, 1.);
	B = revStep * K(x + halfStep, y, 1.);
	D = revStep * K(x, y - halfStep, 1.);
	E = revStep * K(x, y + halfStep, 1.);
	left((Num - 2) * (Num - 3), (Num - 2) * (Num - 3)) = -(A + B + D + E);	//C
	left((Num - 2) * (Num - 3), (Num - 2) * (Num - 3) + 1) = E;				//E
	left((Num - 2) * (Num - 3), (Num - 2) * (Num - 3) - (Num - 2)) = A;		//A
	right((Num - 2) * (Num - 3)) += -D * borderW(Num - 3) - B * borderN(1);

	for (int i = (Num - 2) * (Num - 3) + 1; i < (Num - 2) * (Num - 2) - 1; i++)
	{
		y += h;
		//Последняя строка
		A = revStep * K(x - halfStep, y, 1.);
		B = revStep * K(x + halfStep, y, 1.);
		D = revStep * K(x, y - halfStep, 1.);
		E = revStep * K(x, y + halfStep, 1.);
		left(i, i) = -(A + B + D + E);	//C
		left(i, i - 1) = D;				//D
		left(i, i + 1) = E;				//E
		left(i, i - (Num - 2)) = A;		//A
		right(i) += -B * borderN(i - (Num - 2) * (Num - 3) + 1);
	}

	y += h;
	//Правый нижний угол
	A = revStep * K(x - halfStep, y, 1.);
	B = revStep * K(x + halfStep, y, 1.);
	D = revStep * K(x, y - halfStep, 1.);
	E = revStep * K(x, y + halfStep, 1.);
	left((Num - 2) * (Num - 2) - 1, (Num - 2) * (Num - 2) - 1) = -(A + B + D + E);	//C
	left((Num - 2) * (Num - 2) - 1, (Num - 2) * (Num - 2) - 1 - 1) = D;				//D
	left((Num - 2) * (Num - 2) - 1, (Num - 2) * (Num - 2) - 1 - (Num - 2)) = A;		//A
	right((Num - 2) * (Num - 2) - 1) += -B * borderN(Num - 2) - E * borderE(Num - 3);

	//std::cout << "left:\n" << left << "\n";
	//std::cout << "right:\n" << right << "\n";
}

double RightF(double xx, double yy)
{
	return 4. * xx * yy * cos(xx * yy) - pow(xx * xx + yy * yy, 2) * sin(xx * yy);			//u=sin(xy)+1 //k(x,y)=x^2+y^2
}

double K(double xx, double yy, double uu)
{
	return xx * xx + yy * yy;
}

double ExactSolution(double xx, double yy)
{
	return sin(xx * yy) + 1.;
}

double Error(const Eigen::VectorXd& NumSol, const Eigen::VectorXd& ExSol)
{
	double max = 0., temp;
	for (int i = 0; i < NumSol.size(); i++)
	{
		temp = abs(NumSol(i) - ExSol(i));
		if (temp > max)
			max = temp;
	}
	return max;
}
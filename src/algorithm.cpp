#include "LinkedList.h"

void Hill_Climbing(std::vector<int> vec_int, int& value_max, int& posi_max)
{
	//可能找到局部最優解
	int Bigger = 0;
	int posi = -1;
	value_max = -1;
	//若算出比最大值小5代表正在走下坡，下坡走五次後停止
	while (((++posi) < vec_int.size()) && (Bigger < 5)) {
		if (vec_int.at(posi) > value_max) {
			value_max = vec_int.at(posi);
			posi_max = posi;
		}
		else if (vec_int.at(posi) < value_max) {
			//走下坡
			Bigger++;
		}
	}
}

void Simulated_Annealing(std::vector<int> vec_int, int& value_max, int& posi_max) {
	int posi = -1, P_L = 0, P_L_max = 10;
	while ((++posi) < vec_int.size()) {
		if (vec_int.at(posi) > value_max) {
			value_max = vec_int.at(posi);
			posi_max = posi;
		}
		else if (vec_int.at(posi) < value_max) {	//若沒有更優，按一定概率接受
			double rd = rand() / (RAND_MAX + 1.0);//隨機產生概率：0~1

			double d_exp = exp(((double)vec_int.at(posi) - value_max) / ((double)vec_int.at(vec_int.size() - 1) - vec_int.at(0)));//移動機率要越來越小
			if (d_exp > rd)
			{
				P_L++; //一定概率接受次數 
				if (P_L > P_L_max) {
					posi = vec_int.size();
				}
			}
		}
	}
}

void Fibonacci(std::vector<int> vec_int, int& value_max, int& posi_max) {
	int array_value_max = vec_int.size() - 1;

	//新增Fibonacci array
	std::vector<int> Fib_array;
	int Fib_a = 1, Fib_b = 2;
	Fib_array.push_back(Fib_a);
	Fib_array.push_back(Fib_b);
	while ((Fib_a + Fib_b) < array_value_max) {
		Fib_a += Fib_b;
		Fib_array.push_back(Fib_a);
		Fib_b += Fib_a;
		Fib_array.push_back(Fib_b);
	}

	//從n-2的位置開始
	Fib_array.pop_back();
	Fib_array.pop_back();

	//顛倒array，以利後續auto取值
	std::reverse(Fib_array.begin(), Fib_array.end());
	bool dircetion = false;		//方向是否向負
	bool delta_symbol = true;	//前後MTF差之符號(+或-)
	int before_sharpness = 0;

	int posi = 0;
	for (auto array_value : Fib_array) {
		if (dircetion) {
			array_value = -array_value;
		}

		posi += array_value;

		if (vec_int.at(posi) > value_max) {
			value_max = vec_int.at(posi);
			posi_max = posi;
		}

		//前後兩點差值為delta S，若delta S產生變號，則開始反向搜尋
		if (((vec_int.at(posi) - before_sharpness) > 0) != delta_symbol) {
			dircetion = !dircetion;
		}
		delta_symbol = (vec_int.at(posi) - before_sharpness) > 0 ? true : false;
		before_sharpness = vec_int.at(posi);
	}
}

//矩陣乘法 https://docs.microsoft.com/zh-tw/cpp/parallel/amp/walkthrough-matrix-multiplication?view=vs-2019
bool MutipleWithoutAMP(std::vector<std::vector<float>> aMatrix, std::vector<std::vector<float>> bMatrix, std::vector<std::vector<float>>& product)
{
	int inner_max = aMatrix.at(0).size();
	if (inner_max != bMatrix.size())
		return false;

	int row_max = aMatrix.size();
	int col_max = bMatrix.at(0).size();

	product = std::vector<std::vector<float>>(row_max, std::vector<float>(col_max, 0));

	for (int row = 0; row < row_max; row++) {
		for (int col = 0; col < col_max; col++) {
			// Multiply the row of A by the column of B to get the row, column of product.
			for (int inner = 0; inner < inner_max; inner++) {
				product[row][col] += aMatrix[row][inner] * bMatrix[inner][col];
			}
		}
	}

	return true;
}

#include <amp.h>
bool MutipleWithAMP(std::vector<std::vector<float>> aMatrix, std::vector<std::vector<float>> bMatrix, std::vector<std::vector<float>>& product)
{
	int inner_max = aMatrix.at(0).size();
	if (inner_max != bMatrix.size())
		return false;

	int row_max = aMatrix.size();
	int col_max = bMatrix.at(0).size();

	//2d vector -> 1d vector
	std::vector<float> a_1d;
	for (auto tmp : aMatrix)
		for (auto tmp_int : tmp)
			a_1d.push_back(tmp_int);

	std::vector<float> b_1d;
	for (auto tmp : bMatrix)
		for (auto tmp_int : tmp)
			b_1d.push_back(tmp_int);

	std::vector<float> product_1d = std::vector<float>(row_max * col_max, 0);

	Concurrency::array_view<float, 2> a(row_max, inner_max, a_1d);
	Concurrency::array_view<float, 2> b(inner_max, col_max, b_1d);
	Concurrency::array_view<float, 2> productMatrix(row_max, col_max, product_1d);

	Concurrency::parallel_for_each(productMatrix.extent,
		[=](concurrency::index<2> idx) restrict(amp) {
			int row = idx[0];
			int col = idx[1];
			for (int inner = 0; inner < inner_max; inner++) {
				productMatrix[idx] += a(row, inner) * b(inner, col);
			}
		});

	productMatrix.synchronize();

	//1d vector -> 2d vector
	for (int i = 0; i < row_max; i++) {
		std::vector<float> tmp;
		for (int j = 0; j < col_max; j++) {
			tmp.push_back(product_1d[i * col_max + j]);
		}
		product.push_back(tmp);
	}

	return true;
}

//行列式 determinant
// Function to get cofactor of mat[p][q] in temp[][]. n is current 
// dimension of mat[][] 
void getCofactor(std::vector<std::vector<float>> mat, std::vector<std::vector<float>>& temp, int p, int q, int n)
{
	int i = 0, j = 0;

	// Looping for each element of the matrix 
	for (int row = 0; row < n; row++)
	{
		for (int col = 0; col < n; col++)
		{
			//  Copying into temporary matrix only those element 
			//  which are not in given row and column 
			if (row != p && col != q)
			{
				temp[i][j++] = mat[row][col];

				// Row is filled, so increase row index and 
				// reset col index 
				if (j == n - 1)
				{
					j = 0;
					i++;
				}
			}
		}
	}
}

//matrix is NxN matrix
float determinantOfMatrix(std::vector<std::vector<float>> matrix, int n)
{
	float D = 0; // Initialize result 

	//  Base case : if matrix contains single element 
	if (n == 1)
		return matrix[0][0];

	std::vector<std::vector<float>> temp(n, std::vector<float>(n, 0)); // To store cofactors 

	int sign = 1;  // To store sign multiplier 

	 // Iterate for each element of first row 
	for (int f = 0; f < n; f++)
	{
		// Getting Cofactor of mat[0][f] 
		getCofactor(matrix, temp, 0, f, n);
		D += sign * matrix[0][f] * determinantOfMatrix(temp, n - 1);

		// terms are to be added with alternate sign 
		sign = -sign;
	}

	return D;
}

//伴隨矩陣
//because irrational number can't save by float
//so it will exists deviation
void ADJ_matrix(std::vector<std::vector<float>> matrix, int n, std::vector<std::vector<float>>& adj_matrix)//to find adjoint matrix
{
	adj_matrix = std::vector<std::vector<float>>(n, std::vector<float>(n, 0));
	if (n == 1) {
		adj_matrix[0][0] = 1; return;
	}

	int s = 1;
	std::vector<std::vector<float>> temp(n, std::vector<float>(n, 0)); // To store cofactors 
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			//To get cofactor of M[i][j]
			getCofactor(matrix, temp, i, j, n);
			s = ((i + j) % 2 == 0) ? 1 : -1; //sign of adj[j][i] positive if sum of row and column indexes is even.
			adj_matrix[j][i] = (s) * (determinantOfMatrix(temp, n - 1)); //Interchange rows and columns to get the transpose of the cofactor matrix
		}
	}
}

//反矩陣，Inverse Matrix
//because irrational number can't save by float
//so it will exists deviation
bool Inverse_matrix(std::vector<std::vector<float>>matrix, std::vector<std::vector<float>>& inv_matrix)
{
	int n = matrix.size();
	float det = determinantOfMatrix(matrix, n);
	if (det == 0) {
		std::cout << "determinant is 0 , return." << std::endl;
		return false;
	}

	std::vector<std::vector<float>>adj_matrix;
	ADJ_matrix(matrix, matrix.size(), adj_matrix);

	inv_matrix = std::vector<std::vector<float>>(n, std::vector<float>(n, 0));
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			inv_matrix[i][j] = adj_matrix[i][j] / det;

	return true;
}

//Transpose matrix，轉置矩陣
std::vector<std::vector<float>> transpose_matrix(std::vector<std::vector<float>>matrix) 
{
	int row = matrix.size();
	int col = matrix.at(0).size();
	std::vector<std::vector<float>> transpose_mat = std::vector<std::vector<float>>(col, std::vector<float>(row, 0));

	for (int i = 0; i < row; i++)
		for (int j = 0; j < col; j++)
			transpose_mat[j][i] = matrix[i][j];

	return transpose_mat;
}

//(Moore-Penrose) pseudo inverse matrix
bool Pseudo_Inverse_Matrix(std::vector<std::vector<float>> matrix, std::vector<std::vector<float>>& pim)
{
	// M+ = (M.T.dot(M))^(-1).dot(M.T)
	std::vector<std::vector<float>> MT = transpose_matrix(matrix);

	std::vector<std::vector<float>> MT_dotM;
	if (!MutipleWithoutAMP(MT, matrix, MT_dotM))
		return false;

	std::vector<std::vector<float>> inv_MT_dotM;
	if (!Inverse_matrix(MT_dotM, inv_MT_dotM))
		return false;

	if (!MutipleWithoutAMP(inv_MT_dotM, MT, pim))
		return false;

	return true;
}

//離散差分方程 https://blog.csdn.net/qq_24548569/article/details/85838201
class DDEPM {
public:
	DDEPM();
	~DDEPM();

	bool InitDDEPM(std::vector<float> vec_sharpness);
};

bool DDEPM::InitDDEPM(std::vector<float> vec_sharpness)
{
	if (vec_sharpness.size() < 3) {
		std::cout << "the length of sequence is less than 3" << std::endl;
		return false;
	}

	//累加生成運算 Accumulated Generating Operation
	std::vector<int> vec_AGO;
	int xi = 0;
	for (auto i : vec_sharpness) {
		xi += i;
		vec_AGO.push_back(xi);
	}

	//DDE(2,1) 二接離散差分方程式
	std::vector<std::vector<float>> X;
	std::vector<std::vector<float>> Y;

	for (int i = 0; i < vec_AGO.size() - 2; i++) {
		std::vector<float> X_tmp;
		X_tmp.push_back(vec_AGO.at(i));
		X_tmp.push_back(vec_AGO.at(i + 1));
		X.push_back(X_tmp);

		std::vector<float> Y_tmp;
		Y_tmp.push_back(vec_AGO.at(i + 2));
		Y.push_back(Y_tmp);
	}

	//python code : np.linalg.pinv(X.T.dot(X)).dot(X.T).dot(Y)
	//X.T
	std::vector<std::vector<float>> XT = transpose_matrix(X);

	//X.T.dot(X)
	std::vector<std::vector<float>> XT_dotX;
	if (!MutipleWithAMP(XT, X, XT_dotX))
		return false;

	//np.linalg.pinv(X.T.dot(X)) = pim
	std::vector<std::vector<float>> pim;
	if (!Pseudo_Inverse_Matrix(XT_dotX, pim))
		return false;

	//pim.dot(X.T)
	std::vector<std::vector<float>> pim_dotXT;
	if (!MutipleWithAMP(pim, XT, pim_dotXT)) 
		return false;

	//pim.dot(X.T).dot(Y) = phi
	std::vector<std::vector<float>> phi;
	if (!MutipleWithAMP(pim_dotXT, Y, phi))
		return false;

	float a = phi[0][0];
	float b = phi[1][0];
}

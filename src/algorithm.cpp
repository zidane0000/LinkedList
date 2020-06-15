#include "LinkedList.h"

void Hill_Climbing(std::vector<int> vec_int, int& value_max, int& posi_max)
{
	//�i���짽�����u��
	int Bigger = 0;
	int posi = -1;
	value_max = -1;
	//�Y��X��̤j�Ȥp5�N�����b���U�Y�A�U�Y�������ᰱ��
	while (((++posi) < vec_int.size()) && (Bigger < 5)) {
		if (vec_int.at(posi) > value_max) {
			value_max = vec_int.at(posi);
			posi_max = posi;
		}
		else if (vec_int.at(posi) < value_max) {
			//���U�Y
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
		else if (vec_int.at(posi) < value_max) {	//�Y�S�����u�A���@�w���v����
			double rd = rand() / (RAND_MAX + 1.0);//�H�����ͷ��v�G0~1

			double d_exp = exp(((double)vec_int.at(posi) - value_max) / ((double)vec_int.at(vec_int.size() - 1) - vec_int.at(0)));//���ʾ��v�n�V�ӶV�p
			if (d_exp > rd)
			{
				P_L++; //�@�w���v�������� 
				if (P_L > P_L_max) {
					posi = vec_int.size();
				}
			}
		}
	}
}

void Fibonacci(std::vector<int> vec_int, int& value_max, int& posi_max) {
	int array_value_max = vec_int.size() - 1;

	//�s�WFibonacci array
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

	for (auto posi : Fib_array) {
		if (vec_int.at(posi) > value_max) {
			value_max = vec_int.at(posi);
			posi_max = posi;
		}
	}
}
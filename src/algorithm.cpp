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
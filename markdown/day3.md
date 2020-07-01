# 以C++實現DDEPM之路[3]

程式碼參考DDEPM实例 : [https://blog.csdn.net/qq_24548569/article/details/85838201](https://blog.csdn.net/qq_24548569/article/details/85838201)

程式碼實作： [https://github.com/zidane0000/LinkedList/blob/master/src/algorithm.cpp](https://github.com/zidane0000/LinkedList/blob/master/src/algorithm.cpp)

## 實作Moore-Penrose Pseudo Inverse Matrix
1. 實作Singular Value Decomposition

- 實驗後發現，在特定情況下，參考文獻與python上的svd之結果會有不同
- 特定情況為當InitDDEPM時vec_sharpness的size大於3，則結果svd計算出來的U以及V會不同(數字一樣，正負號有差異)
- 經多方查詢svd計算方法發現，共六個(連結已放在參考)，第一個為純數學，二至四為c library，前四個結果大致相同，最後兩個為python library，結果與前三個相異
- 實驗可以矩陣A[[9857,3131],[3131,1470]]為例，矩陣A以vec_sharpness[7,8,19,40,53]生成
- python code 在 https://github.com/zidane0000/LinkedList/blob/master/src/DDEPM.txt
- python online IDE 使用 https://repl.it/languages/python3
- 結論為相信前四個參考

```cpp
bool dsvd(std::vector<std::vector<float>>& a, int m, int n, std::vector<std::vector<float>>& w, std::vector<std::vector<float>>& v)
```

2. 實作class DDEPM

```cpp
class DDEPM {
public:
	DDEPM();
	~DDEPM();

	bool InitDDEPM(std::vector<float> vec_sharpness);
	bool predict(double posi, double& value);
private:
	float a, b;
	int type;
	double r1, r2, C1, C2;
};
```

### 參考

> https://atozmath.com/MatrixEv.aspx?q=svd

> http://svn.lirec.eu/libs/magicsquares/src/SVD.cpp

> http://www.mymathlib.com/c_source/matrices/linearsystems/singular_value_decomposition.c

> https://keisan.casio.com/exec/system/15076953160460

> https://github.com/numpy/numpy/blob/v1.19.0/numpy/linalg/linalg.py#L1483-L1675

> https://github.com/scipy/scipy/blob/v1.5.0/scipy/linalg/decomp_svd.py#L13-L136

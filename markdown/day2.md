# 以C++實現DDEPM之路[2]

程式碼參考DDEPM实例 : [https://blog.csdn.net/qq_24548569/article/details/85838201](https://blog.csdn.net/qq_24548569/article/details/85838201)

程式碼實作： [https://github.com/zidane0000/LinkedList/blob/master/src/algorithm.cpp](https://github.com/zidane0000/LinkedList/blob/master/src/algorithm.cpp)

## 實作Moore-Penrose Pseudo Inverse Matrix
1. 實作轉置矩陣	transpose_matrix

```cpp
std::vector<std::vector<float>> transpose_matrix(std::vector<std::vector<float>>matrix)
```

2. 實作矩陣相乘	mutiply_matrix
- 是否使用c++ amp 加速
- 實驗後發現沒有比較快速，猜測可能是矩陣大小問題

```cpp
bool MutipleWithoutAMP(std::vector<std::vector<float>> aMatrix, std::vector<std::vector<float>> bMatrix, std::vector<std::vector<float>>& product)
bool MutipleWithAMP(std::vector<std::vector<float>> aMatrix, std::vector<std::vector<float>> bMatrix, std::vector<std::vector<float>>& product)
```

3. 實作餘因子矩陣 cofactor matrix

```cpp
void getCofactor(std::vector<std::vector<float>> mat, std::vector<std::vector<float>>& temp, int p, int q, int n)
```

4. 實作行列式 determinant

```cpp
float determinantOfMatrix(std::vector<std::vector<float>> matrix , int n)
```

5. 實作伴隨矩陣 ajoint matrix

```cpp
void ADJ_matrix(std::vector<std::vector<float>> matrix, int n, std::vector<std::vector<float>>& adj_matrix)//to find adjoint matrix
```

6. 實作反矩陣	inverse matrix

```cpp
typedef std::pair<int, int> Double_int;
```

7. 求Pseudo_Inverse_Matrix(有誤，需變換解法)

```cpp
bool Pseudo_Inverse_Matrix(std::vector<std::vector<float>> matrix, std::vector<std::vector<float>>& pim)
```

8. 根據np.linalg.pinv，是使用SVD(Singular value decomposition)來解，放在第三天

下次見

### 參考

> https://www.geeksforgeeks.org/adjoint-inverse-matrix/
> https://docs.microsoft.com/zh-tw/cpp/parallel/amp/walkthrough-matrix-multiplication?view=vs-2019
> https://fractalytics.io/moore-penrose-matrix-optimization-cuda-c 程式有誤，僅供參考!!!
> https://numpy.org/doc/stable/reference/generated/numpy.linalg.pinv.html
> https://github.com/numpy/numpy/blob/v1.18.4/numpy/linalg/linalg.py#L1881-L1970

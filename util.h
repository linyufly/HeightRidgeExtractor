// Author: Mingcheng Chen (linyufly@gmail.com)

#ifndef UTIL_H_
#define UTIL_H_

template<class T>
T **create_matrix(int num_of_rows, int num_of_cols) {
  T *data_array = new T[num_of_rows * num_of_cols];
  T **matrix = new T *[num_of_rows];

  for (int i = 0; i < num_of_rows; i++) {
    matrix[i] = data_array + i * num_of_cols;
  }

  return matrix;
}

template<class T>
void delete_matrix(T **matrix) {
  delete [] matrix[0];
  delete [] matrix;
}

#endif  // UTIL_H_

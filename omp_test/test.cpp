#include <omp.h>
#include <cmath>
#include <iostream>

int main(int argc, char** argv) {
  int N = 1000;
  int nthread = 1;

  if (argc > 1) {
    N = atoi(argv[1]);
  }

  if (argc > 2) {
    nthread = atoi(argv[2]);
  }

  omp_set_num_threads(nthread);

  double res = 0;
  double dx = 1 / double(N);

  #pragma omp parallel
  {
    #pragma omp for
    for (int i = 0; i < N; ++i) {
      #pragma omp atomic
      res += 1 / (1. + std::pow(double(i) * dx, 2));
    }
  }

  res *= 4 * dx;

  double pi = std::acos(-1);
  std::cout << std::endl << (res - pi) / pi << std::endl;

  return 0;
}

// // #include <omp.h>
// static long num_steps = 100000;
// double step;
// #define NUM_THREADS 10

// int main() {
//   int i, nthreads;
//   double pi, sum[NUM_THREADS];
//   step = 1.0 / (double)num_steps;
//   omp_set_num_threads(NUM_THREADS);
// #pragma omp parallel
//   {
//     int i, id, nthrds;
//     double x;
//     id = omp_get_thread_num();
//     nthrds = omp_get_num_threads();
//     if (id == 0) nthreads = nthrds;
//     for (i = id, sum[id] = 0.0; i < num_steps; i = i + nthrds) {
//       x = (i + 0.5) * step;
//       sum[id] += 4.0 / (1.0 + x * x);
//     }
//   }
//   for (i = 0, pi = 0.0; i < nthreads; i++) pi += sum[i] * step;

//   double real_pi = std::acos(-1);
//     std::cout << std::endl << (pi-real_pi)/real_pi << std::endl;

//   return 0;
// }

// #include <filesystem>
// #include <iostream>

// int main() {

//   std::string("hello world");
//   std::cout << "Hello world" << std::endl;

//   return 0;
// }

// #include <iostream>
// #include <armadillo>
// #include <vector>
// #include <array>
// #include <valarray>
// #include <chrono>

// #define n 100000000

// int main(){
//   std::time_t start_time, end_time;

//   {
//     std::vector<double> v1(3,0);
//     std::vector<double> v2(3, 0);

//     start_time = std::time(nullptr);
//     for (int i=0; i<n; ++i){
//       std::vector<double> v3(3,0);
//       for (int j=0; j<v1.size(); ++j){
//         v3[j] = v1[j]+v2[j];
//       }
//     }
//     end_time = std::time(nullptr);
//     std::cout << "\nvector: " << std::difftime(end_time, start_time) << " seconds" << std::endl;
//   }

//   {
//     std::valarray<double> v1(3, 0);
//     std::valarray<double> v2(3, 0);

//     start_time = std::time(nullptr);
//     for (int i = 0; i < n; ++i) {
//       std::valarray<double> v3 = v1 + v2;
//     }
//     end_time = std::time(nullptr);
//     std::cout << "\nvalarray: " << std::difftime(end_time, start_time) << " seconds" << std::endl;
//   }

//   {
//     std::array<double,3> v1{0,0,0};
//     std::array<double,3> v2{0,0,0};

//     start_time = std::time(nullptr);
//     for (int i = 0; i < n; ++i) {
//       std::array<double, 3> v3{0,0,0};
//       for (int j = 0; j < v1.size(); ++j) {
//         v3[j] = v1[j] + v2[j];
//       }
//     }
//     end_time = std::time(nullptr);
//     std::cout << "\narray: " << std::difftime(end_time, start_time) << " seconds" << std::endl;
//   }

//   {
//     arma::vec v1{0,0,0};
//     arma::vec v2{0, 0, 0};

//     start_time = std::time(nullptr);
//     for (int i = 0; i < n; ++i) {
//       arma::vec v3 = v1 + v2;
//     }
//     end_time = std::time(nullptr);
//     std::cout << "\narma: " << std::difftime(end_time, start_time) << " seconds" << std::endl;
//   }

//   return 0;
// }
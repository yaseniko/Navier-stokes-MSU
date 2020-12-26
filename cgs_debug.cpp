#include <iostream>
#include <math.h>
#include <memory>
#include <stdio.h>
#include <stdlib.h>

using std::unique_ptr;

constexpr double precision = 1e-9;

int solve_gauss (uint size, double *non_zeroes, double *rhs)
{
  for (int i = 0; i < size; i++)
    {
      double tmp = non_zeroes[i * size + i];
      for (int j = i; j < size; j++)
        non_zeroes[i * size + j] /= tmp;
      rhs[i] /= tmp;

      for (int j = i + 1; j < size; j++)
        {
          for (int k = i + 1; k < size; k++)
            non_zeroes[j * size + k] -= non_zeroes[i * size + k] * non_zeroes[j * size + i];

          rhs[j] -= rhs[i] * non_zeroes[j * size + i];
        }
    }

  for (int i = size - 1; i >= 0; i--)
    {
      for (int j = i + 1; j < size; j++)
        rhs[i] -= rhs[j] * non_zeroes[i * size + j];
    }

  return 0;
}

void copy_vector (uint size, double *src, double *dst)
{
  for (uint i = 0; i < size; i++)
    dst[i] = src[i];
}

void set_ort_vector (int size, double *v)
{
  v[0] = 1;
  for (int i = 1; i < size; i++)
    v[i] = 0;
}

void set_unit_vector (int size, double *v)
{
  uint norm = size * (size + 1) * (2 * size + 1) / 6;
  for (int i = 0; i < size; i++)
    v[i] = (i + 1.) / size;
}

void print_vector (uint size, double *v)
{
  for (uint i = 0; i < size; i++)
    printf ("%.2f ", v[i]);
  printf("\n");
}

double sc_prod (uint size, double *v_1, double *v_2)
{
  double s_p = 0.;
  for (uint i = 0; i < size; i++)
    s_p += v_1[i] * v_2[i];

  return s_p;
}

// res = b + alpha * a
void lin_operation_1 (uint size, double alpha, double *b, double *a, double *res)
{
  for (uint i = 0; i < size; i++)
    res[i] = b[i] + alpha * a[i];
}

// res = b + alpha * a + beta * c
void lin_operation_2 (uint size, double alpha, double beta, double *b, double *a, double *c, double *res)
{
  for (uint i = 0; i < size; i++)
    res[i] = b[i] + alpha * a[i] + beta * c[i];
}

void matrix_mult_vector_csr (uint size, double *non_zeroes, int *cols, int *ind, double *vect, double *res)
{
  for (uint i = 0; i < size; i++)
    {
      res[i] = 0.;
      for (int j = ind[i]; j < ind[i + 1]; j++)
        res[i] += non_zeroes[j] * vect[cols[j]];
    }
}

//res = A * (v_1 + v_2)
void matrix_mult_sum_vectors_csr (uint size, double *non_zeroes, int *cols, int *ind, double *v_1, double *v_2, double *res)
{
  unique_ptr<double []> sum (new double[size]);
  lin_operation_1 (size, 1., v_1, v_2, sum.get ());
  matrix_mult_vector_csr (size, non_zeroes, cols, ind, sum.get (), res);
}

int cgs_solve_csr_system (uint size, double *non_zeroes, int *cols, int *ind, double *rhs, double *solution)
{
  unique_ptr<double []> z_a (new double[size]);
  unique_ptr<double []> d (new double[size]);
  unique_ptr<double []> Ap (new double[size]);
  unique_ptr<double []> p (new double[size]);
  unique_ptr<double []> r (new double[size]);
  unique_ptr<double []> r_next (new double[size]);
  unique_ptr<double []> u (new double[size]);
  unique_ptr<double []> q (new double[size]);
  unique_ptr<double []> tmp (new double[size]);

  set_unit_vector (size, z_a.get ());
  copy_vector (size, solution, d.get ());

  // Ad (tmp) = A * d
  matrix_mult_vector_csr (size, non_zeroes, cols, ind, d.get (), tmp.get ());
  // r = rhs - Ad (tmp)
  lin_operation_1 (size, -1., rhs, tmp.get (), r.get ());
  // p = r
  copy_vector (size, r.get (), p.get ());
  // u = r
  copy_vector (size, r.get (), u.get ());

  double begin_norm = sqrt (sc_prod (size, r.get (), r.get ()));
  std::cout << "begin residual: " << begin_norm << std::endl;
  if (begin_norm < precision)
    {
      std::cout << "WARNING! Initial approximation is solution!" << std::endl;
      return 0;
    }

  unsigned int max_iter = 1000;
  double alpha, beta, norm;
  uint iter;
  for (iter = 0; iter < max_iter; iter++)
    {
      // Ap = A * p
      matrix_mult_vector_csr (size, non_zeroes, cols, ind, p.get (), Ap.get ());
      // alpha = (r, z_a) / (Ap, z_a)
      alpha = sc_prod (size, r.get (), z_a.get ()) / sc_prod (size, Ap.get (), z_a.get ());
      // q = u - alpha * Ap
      lin_operation_1 (size, -alpha, u.get (), Ap.get (), q.get ());
      // d = d + alpha (u + q)
      lin_operation_2 (size, alpha, alpha, d.get (), u.get (), q.get (), d.get ());
      // tmp = A * (u + q)
      matrix_mult_sum_vectors_csr (size, non_zeroes, cols, ind, u.get (), q.get (), tmp.get ());
      // r_next = r - alpha * A * (u + q)
      lin_operation_1 (size, -alpha, r.get (), tmp.get (), r_next.get ());
      // beta = (r_next, z_a) / (r, z_a)
      beta = sc_prod (size, r_next.get (), z_a.get ()) / sc_prod (size, r.get (), z_a.get ());
      // u = r_next + beta * q
      lin_operation_1 (size, beta, r_next.get (), q.get (), u.get ());
      // p = u + beta * (q + beta * p)
      lin_operation_2 (size, beta, beta * beta, u.get (), q.get (), p.get (), p.get ());
      // r = r_next
      copy_vector (size, r_next.get (), r.get ());
      norm = sqrt (sc_prod (size, r.get (), r.get ()));
      std::cout << "iter = " << iter << ", residual: " << norm << std::endl;
      if (fabs (norm / begin_norm) < precision)
        break;
    }

  if (iter == max_iter)
    {
      std::cout << "ERROR! NOT CONVERGED" << std::endl;
      return 1;
    }

  copy_vector (size, d.get (), solution);
  return 0;
}

/*int main ()
{
  int size = 5;
  double min_for_divison = 1e-64;

  unique_ptr<double []> rnz (new double[size * size]);
  unique_ptr<int []> cols (new int[size * size]);
  unique_ptr<int []> ind (new int[size + 1]);
  unique_ptr<double []> rhs (new double[size]);
  unique_ptr<double []> sol (new double[size]);
  unique_ptr<double []> diag (new double[size]);

  // fill matrix
  ind[0] = 0;
  for (int i = 0; i < size; i++)
    {
      ind[i + 1] = ind[i] + size;
      for (int j = 0; j < size; j++)
        {
          rnz[i * size + j] = fabs (i - j + 2);
          cols[i * size + j] = j;
        }
    }

  // fill rhs
  for (int i = 0; i < size; i++)
    {
      rhs[i] = 0.;
      double coeff = 1.33;
      for (int j = 0; j < size; j++)
        {
          rhs[i] += fabs (i - j + 2) * coeff;
          coeff += 0.5;
        }
    }

  set_unit_vector (size, sol.get ());
  cgs_solve_csr_system (size, rnz.get (), cols.get (), ind.get (), rhs.get (), sol.get ());
  print_vector (size, sol.get ());
  return 0;
}*/

int main ()
{
  int size = 20;

  unique_ptr<double []> rnz (new double[size * size]);
  unique_ptr<double []> rhs (new double[size]);
  unique_ptr<double []> sol (new double[size]);

  for (int i = 0; i < size; i++)
    {
      rhs[i] = 0.;
      for (int j = 0; j < size; j++)
        {
          double tmp = fabs (i - j + 1);
          rnz[i * size + j] = tmp;
          rhs[i] += (j % 2 == 0) ? 3 * tmp : -2 * tmp;
        }
    }

  solve_gauss (size, rnz.get (), rhs.get ());
  print_vector (size, rhs.get ());
  return 0;
}

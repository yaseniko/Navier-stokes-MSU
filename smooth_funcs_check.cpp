#include <iostream>
#include <math.h>
#include <memory>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

#define FIX_UNUSED(X) (void) (X);

using std::unique_ptr;

constexpr double two_pi = 2. * M_PI;
constexpr double four_sq_pi = 4. * M_PI * M_PI;

constexpr int max_print = 20;
constexpr int max_iter = 1000;
constexpr double precision = 1e-9;
constexpr double min_cmp = 1e-16;
constexpr double rhs_coeff = 0.01;

constexpr double len_t = 1.;
constexpr double len_x = 3.;
constexpr double len_y = 3.;
constexpr double p_ro = 10.;
constexpr double mu = 0.01;
constexpr double mu_div_three = mu / 3.;
constexpr double four_mu_div_three = 4. * mu_div_three;


double ro (double t, double x_1, double x_2)
{
  return exp (t) * (cos (two_pi * x_1) + 1.5) * (sin (two_pi * x_2) + 1.5) * rhs_coeff;
}

double ro_t (double t, double x_1, double x_2)
{
  return ro (t, x_1, x_2);
}

double ro_x (double t, double x_1, double x_2)
{
  return -two_pi * exp (t) * sin (two_pi * x_1) * (sin (two_pi * x_2) + 1.5) * rhs_coeff;
}

double ro_y (double t, double x_1, double x_2)
{
  return +two_pi * exp (t) * (cos (two_pi * x_1) + 1.5) * cos (two_pi * x_2) * rhs_coeff;
}

double u_1 (double t, double x_1, double x_2)
{
  return exp (t) * sin (two_pi * x_1) * sin (two_pi * x_2) * rhs_coeff;
}

double u_1_t (double t, double x_1, double x_2)
{
  return u_1 (t, x_1, x_2);
}

double u_1_x (double t, double x_1, double x_2)
{
  return two_pi * exp (t) * cos (two_pi * x_1) * sin (two_pi * x_2) * rhs_coeff;
}

double u_1_y (double t, double x_1, double x_2)
{
  return two_pi * exp (t) * sin (two_pi * x_1) * cos (two_pi * x_2) * rhs_coeff;
}

double u_1_xx (double t, double x_1, double x_2)
{
  return -four_sq_pi * u_1 (t, x_1, x_2);
}

double u_1_xy (double t, double x_1, double x_2)
{
  return four_sq_pi * exp (t) * cos (two_pi * x_1) * cos (two_pi * x_2) * rhs_coeff;
}

double u_1_yy (double t, double x_1, double x_2)
{
  return -four_sq_pi * u_1 (t, x_1, x_2);
}

double u_2 (double t, double x_1, double x_2)
{
  return exp (-t) * sin (two_pi * x_1) * sin (two_pi * x_2) * rhs_coeff;
}

double u_2_t (double t, double x_1, double x_2)
{
  return -u_2 (t, x_1, x_2);
}

double u_2_x (double t, double x_1, double x_2)
{
  return two_pi * exp (-t) * cos (two_pi * x_1) * sin (two_pi * x_2) * rhs_coeff;
}

double u_2_y (double t, double x_1, double x_2)
{
  return two_pi * exp (-t) * sin (two_pi * x_1) * cos (two_pi * x_2) * rhs_coeff;
}

double u_2_xx (double t, double x_1, double x_2)
{
  return -four_sq_pi * u_2 (t, x_1, x_2);
}

double u_2_xy (double t, double x_1, double x_2)
{
  return four_sq_pi * exp (-t) * cos (two_pi * x_1) * cos (two_pi * x_2) * rhs_coeff;
}

double u_2_yy (double t, double x_1, double x_2)
{
  return -four_sq_pi * u_2 (t, x_1, x_2);
}

double f_0 (double t, double x_1, double x_2)
{
  double rho = ro (t, x_1, x_2);
  double rho_t = ro_t (t, x_1, x_2);
  double rho_x = ro_x (t, x_1, x_2);
  double rho_y = ro_y (t, x_1, x_2);
  double u = u_1 (t, x_1, x_2);
  double u_x = u_1_x (t, x_1, x_2);
  double v = u_2 (t, x_1, x_2);
  double v_y = u_2_y (t, x_1, x_2);

  return rho_t + u * rho_x + rho * u_x + v * rho_y + rho * v_y;
}

double f_1 (double t, double x_1, double x_2)
{
  double rho = ro (t, x_1, x_2);
  double rho_t = ro_t (t, x_1, x_2);
  double rho_x = ro_x (t, x_1, x_2);
  double rho_y = ro_y (t, x_1, x_2);
  double u = u_1 (t, x_1, x_2);
  double u_t = u_1_t (t, x_1, x_2);
  double u_x = u_1_x (t, x_1, x_2);
  double u_y = u_1_y (t, x_1, x_2);
  double v = u_2 (t, x_1, x_2);
  double v_x = u_2_x (t, x_1, x_2);
  double v_y = u_2_y (t, x_1, x_2);

  return (u * rho_t + rho * u_t + u * u * rho_x + 2. * rho * u * u_x +
          u * v * rho_y + u * rho * v_y + rho * v * u_y + p_ro * rho_x -
          mu * u_1_yy (t, x_1, x_2) - four_mu_div_three * u_1_xx (t, x_1, x_2) - mu_div_three * u_2_xy (t, x_1, x_2)) / rho;
}

double f_2 (double t, double x_1, double x_2)
{
  double rho = ro (t, x_1, x_2);
  double rho_t = ro_t (t, x_1, x_2);
  double rho_x = ro_x (t, x_1, x_2);
  double rho_y = ro_y (t, x_1, x_2);
  double u = u_1 (t, x_1, x_2);
  double u_x = u_1_x (t, x_1, x_2);
  double u_y = u_1_y (t, x_1, x_2);
  double v = u_2 (t, x_1, x_2);
  double v_t = u_2_t (t, x_1, x_2);
  double v_x = u_2_x (t, x_1, x_2);
  double v_y = u_2_y (t, x_1, x_2);

  return (v * rho_t + rho * v_t + v * v * rho_y + 2. * rho * v * v_y +
          u * v * rho_x + u * rho * v_x + rho * v * u_x + p_ro * rho_y -
          mu * u_2_xx (t, x_1, x_2) - four_mu_div_three * u_2_yy (t, x_1, x_2) - mu_div_three * u_1_xy (t, x_1, x_2)) / rho;
}

void copy_vector (int size, double *src, double *dst)
{
  for (int i = 0; i < size; i++)
    dst[i] = src[i];
}

void set_zero_vector (int size, double *v)
{
  for (int i = 0; i < size; i++)
    v[i] = 0.;
}

void flush_small_values (int size, double *v)
{
  for (int i = 0; i < size; i++)
    if (fabs (v[i]) < min_cmp)
      v[i] = 0.;
}

void print_vector (int size, double *v, bool if_debug = false)
{
  int print_size = if_debug ? size : std::min (size, max_print);
  for (int i = 0; i < print_size; i++)
    {
      if (i % static_cast<int> (sqrt (size)) == 0)
        printf ("\n");
      printf ("%6.2f ", v[i]);
    }

  printf("\n");
}


void set_ort_vector (int size, double *v)
{
  v[0] = 1.;
  for (int i = 1; i < size; i++)
    v[i] = 0.;
}

void set_unit_vector (int size, double *v)
{
  int norm = size * (size + 1) * (2 * size + 1) / 6;
  int sign = 1;
  for (int i = 0; i < size; i++)
    {
      v[i] = sign * (i + 1.) / size;
      sign = -sign;
    }
}

void print_csr_matr (int size, double *rnz, int *cols, int *ind)
{
  int nz = 0;
  for (int i = 0; i < size; i++)
    {
      for (int j = ind[i]; j < ind[i + 1]; j++)
        {
          printf ("m[%d][%d] = %7.2e ", i, cols[nz], rnz[nz]);
          nz++;
        }
      printf ("\n");
    }
}

double sc_prod (int size, double *v_1, double *v_2)
{
  double s_p = 0.;
  for (int i = 0; i < size; i++)
    s_p += v_1[i] * v_2[i];

  return s_p;
}

// res = b + alpha * a
void lin_operation_1 (int size, double alpha, double *b, double *a, double *res)
{
  for (int i = 0; i < size; i++)
    res[i] = b[i] + alpha * a[i];
}

// res = b + alpha * a + beta * c
void lin_operation_2 (int size, double alpha, double beta, double *b, double *a, double *c, double *res)
{
  for (int i = 0; i < size; i++)
    res[i] = b[i] + alpha * a[i] + beta * c[i];
}

void matrix_mult_vector_csr (int size, double *non_zeroes, int *cols, int *ind, double *vect, double *res)
{
  for (int i = 0; i < size; i++)
    {
      res[i] = 0.;
      for (int j = ind[i]; j < ind[i + 1]; j++)
        res[i] += non_zeroes[j] * vect[cols[j]];
    }
}

//res = A * (v_1 + v_2)
void matrix_mult_sum_vectors_csr (int size, double *non_zeroes, int *cols, int *ind, double *v_1, double *v_2, double *res)
{
  unique_ptr<double []> sum (new double[size]);
  lin_operation_1 (size, 1., v_1, v_2, sum.get ());
  matrix_mult_vector_csr (size, non_zeroes, cols, ind, sum.get (), res);
}

void system_residual (int size, double *non_zeroes, int *cols, int *ind, double *rhs, double *solution)
{
  unique_ptr<double []> Ax (new double[size]);
  matrix_mult_vector_csr (size, non_zeroes, cols, ind, solution, Ax.get ());

  double res = 0;
  double value = 0;
  for (int i = 0; i < size; i++)
    {
      res = std::max (res, fabs (rhs[i] - Ax[i]));
    }

  if (fabs (res) > 1e-5)
    printf ("CGS residual: %e\n", res);
}


int cgs_solve_csr_system (int size, double *non_zeroes, int *cols, int *ind, double *rhs, double *solution)
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
  //set_unit_vector (size, d.get ());
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

#if 0
  std::cout << "begin residual: " << begin_norm << std::endl;
#endif

  if (begin_norm < precision)
    {
      //std::cout << "WARNING! Initial approximation is solution!" << std::endl;
      return 0;
    }

  double alpha, beta, norm;
  int iter;
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
#if 0
      std::cout << "iter = " << iter << ", residual: " << norm << std::endl;
#endif
      if (fabs (norm / begin_norm) < precision)
        break;
    }

  if (iter == max_iter)
    {
      //std::cout << "ERROR! NOT CONVERGED" << std::endl;
      return 1;
    }
#if 0
  else
    {
      std::cout << "CONVERGED" << std::endl;
    }
#endif

  copy_vector (size, d.get (), solution);
#if 0
  system_residual (size, non_zeroes, cols, ind, rhs, solution);
#endif

  return 0;
}

void fill_grid (int m_x, int m_y, double h_x, double h_y, bool *is_rho_valid,
                int *st, int *lt, int *rt, double *x_coord, double *y_coord)
{
  int nodes_before_hole = (m_y + 1) * (m_x / 3 + 1);
  int nodes_near_hole = (m_y / 3 + 1) * (m_x / 3 - 1);

  for (int curr_node = 0; curr_node < nodes_before_hole; curr_node++)
    {
      int row = curr_node % (m_y + 1);
      int col = curr_node / (m_y + 1);

      x_coord[curr_node] = h_x * col;
      y_coord[curr_node] = h_y * row;

      if (col == 0)
        {
          lt[curr_node] = -1;
          rt[curr_node] = curr_node + (m_y + 1);

          if (row == 0)
            {
              is_rho_valid[curr_node] = true;
              st[curr_node] = 5;
            }
          else if (row == m_y)
            {
              is_rho_valid[curr_node] = false;
              st[curr_node] = 7;
            }
          else
            {
              is_rho_valid[curr_node] = true;
              st[curr_node] = 17;
            }
        }
      else if (col == m_x / 3)
        {
          lt[curr_node] = curr_node - (m_y + 1);

          if (row == 0)
            {
              is_rho_valid[curr_node] = false;
              st[curr_node] = 13;
              rt[curr_node] = -1;
            }
          else if (row == m_y)
            {
              is_rho_valid[curr_node] = false;
              st[curr_node] = 15;
              rt[curr_node] = -1;
            }
          else
            {
              if (row == m_y / 3)
                {
                  is_rho_valid[curr_node] = true;
                  st[curr_node] = 9;
                  rt[curr_node] = nodes_before_hole;
                }
              else if (row > m_y / 3 && row < 2 * m_y / 3)
                {
                  is_rho_valid[curr_node] = true;
                  st[curr_node] = 0;
                  rt[curr_node] = rt[curr_node - 1] + 1;
                }
              else if (row == 2 * m_y / 3)
                {
                  is_rho_valid[curr_node] = false;
                  st[curr_node] = 11;
                  rt[curr_node] = rt[curr_node - 1] + 1;
                }
              else
                {
                  is_rho_valid[curr_node] = false;
                  st[curr_node] = 4;
                  rt[curr_node] = -1;
                }
            }
        }
      else
        {
          lt[curr_node] = curr_node - (m_y + 1);
          rt[curr_node] = curr_node + (m_y + 1);

          if (row == 0)
            {
              is_rho_valid[curr_node] = true;
              st[curr_node] = 1;
            }
          else if (row == m_y)
            {
              is_rho_valid[curr_node] = false;
              st[curr_node] = 2;
            }
          else
            {
              is_rho_valid[curr_node] = true;
              st[curr_node] = 0;
            }
        }
    }

  for (int curr_node = 0; curr_node < nodes_near_hole; curr_node++)
    {
      int row = curr_node % (m_y / 3 + 1);
      int col = curr_node / (m_y / 3 + 1);
      int node_id = nodes_before_hole + curr_node;

      x_coord[node_id] = len_x / 3 + h_x * (col + 1);
      y_coord[node_id] = len_y / 3 + h_y * row;

      lt[node_id] = col == 0 ? node_id - (2 * m_y / 3 + 1) : node_id - (m_y / 3 + 1);
      rt[node_id] = col == (m_x / 3 - 2) ? node_id + (2 * m_y / 3 + 1) : node_id + (m_y / 3 + 1);

      if (row == 0)
        {
          is_rho_valid[node_id] = true;
          st[node_id] = 1;
        }
      else if (row == m_y / 3)
        {
          is_rho_valid[node_id] = false;
          st[node_id] = 2;
        }
      else
        {
          is_rho_valid[node_id] = true;
          st[node_id] = 0;
        }
    }

  for (int curr_node = 0; curr_node < nodes_before_hole; curr_node++)
    {
      int row = curr_node % (m_y + 1);
      int col = curr_node / (m_y + 1);
      int node_id = nodes_before_hole + nodes_near_hole + curr_node;

      x_coord[node_id] = 2 * len_x / 3 + h_x * col;
      y_coord[node_id] = h_y * row;

      if (col == m_x / 3)
        {
          lt[node_id] = node_id - (m_y + 1);
          rt[node_id] = -1;
          is_rho_valid[node_id] = false;

          if (row == 0)
            st[node_id] = 6;
          else if (row == m_y)
            st[node_id] = 8;
          else
            st[node_id] = 4;
        }
      else if (col == 0) // done
        {
          rt[node_id] = node_id + (m_y + 1);

          if (row == 0)
            {
              is_rho_valid[node_id] = true;
              st[node_id] = 14;
              lt[node_id] = -1;
            }
          else if (row == m_y)
            {
              is_rho_valid[node_id] = false;
              st[node_id] = 16;
              lt[node_id] = -1;
            }
          else
            {
              is_rho_valid[node_id] = true;

              if (row == m_y / 3)
                {
                  st[node_id] = 10;
                  lt[node_id] = node_id - (2 * m_y / 3 + 1);
                }
              else if (row > m_y / 3 && row < 2 * m_y / 3)
                {
                  st[node_id] = 0;
                  lt[node_id] = lt[node_id - 1] + 1;
                }
              else if (row == 2 * m_y / 3)
                {
                  st[node_id] = 12;
                  lt[node_id] = lt[node_id - 1] + 1;
                }
              else
                {
                  st[node_id] = 3;
                  lt[node_id] = -1;
                }
            }
        }
      else
        {
          lt[node_id] = node_id - (m_y + 1);
          rt[node_id] = node_id + (m_y + 1);

          if (row == 0)
            {
              is_rho_valid[node_id] = true;
              st[node_id] = 18;
            }
          else if (row == m_y)
            {
              is_rho_valid[node_id] = false;
              st[node_id] = 2;
            }
          else
            {
              is_rho_valid[node_id] = true;
              st[node_id] = 0;
            }
        }
    }
}


void check_grid (int dim, bool *is_rho_valid, int *st,
                 int *lt, int *rt, double *x_coord, double *y_coord)
{
  printf("\ndim = %d\n", dim);
  for (int id = 0; id < dim; id++)
    printf("%d %d %d %d %d %f %f\n", id, is_rho_valid[id], st[id],
           lt[id], rt[id], x_coord[id], y_coord[id]);
  printf ("\n");
}

void fill_zero_layer (int dim, double h_x, double h_y, double *v_1, double *v_2,
                      bool *is_rho_valid, double *h, double *x_coord, double *y_coord)
{
  for (int id = 0; id < dim; id++)
    {
      v_1[id] = u_1 (0., x_coord[id], y_coord[id]);
      v_2[id] = u_2 (0., x_coord[id], y_coord[id]);
      h[id] = is_rho_valid[id] ? ro (0., x_coord[id] + h_x / 2., y_coord[id] + h_y / 2.) : 0.;
    }
}

void build_sle_for_v_1 (int dim, double curr_t, double h_x, double h_y, double tau,
                        int *st, double *x_coord, double *y_coord,
                        double *v_1, double *v_2, double *h, int *lt, int *rt,
                        double *non_zeroes, int *cols, int *ind, double *b)
{
  double tau_div_hx = tau / h_x;
  double tau_div_hy = tau / h_y;
  double tau_div_two_hx = tau / (2. * h_x);
  double tau_div_two_hy = tau / (2. * h_y);
  double four_tau_mu_div_three_sq_hx = (four_mu_div_three * tau) / (h_x * h_x);
  double tau_mu_div_sq_hy = tau * mu / (h_y * h_y);
  double tau_mu_div_twelve_hx_hy = tau * mu / (12. * h_x * h_y);

  auto inner = [&] (int id) -> bool { return st[id] == 0; };

  auto H1 = [&] (int id) -> double { return (h[id] + h[id - 1]) / 2.; };
  auto HT = [&] (int id) -> double { return (h[id] + h[id - 1] + h[lt[id]] + h[lt[id - 1]]) / 4.; };
  auto f1 = [&] (int id) -> double { return f_1 (curr_t, x_coord[id], y_coord[id]); };

  auto v_plus_fabs_v  = [&] (double v) -> double {return v + fabs (v); };
  auto v_minus_fabs_v = [&] (double v) -> double {return v - fabs (v); };

  int non_zeroes_count = 0;
  ind[0] = 0;

  for (int id = 0; id < dim; id++)
    {
      if (!inner (id) || fabs (HT (id)) < min_cmp)
        {
          non_zeroes[non_zeroes_count] = 1.;
          cols[non_zeroes_count++] = id;
          ind[id + 1] = ind[id] + 1;
          b[id] = 0.;
        }
      else
        {
          ind[id + 1] = ind[id];

          double H_T = HT (id);
          non_zeroes[non_zeroes_count] = -tau_div_two_hx * v_plus_fabs_v (v_1[id]) * H_T - four_tau_mu_div_three_sq_hx;
          cols[non_zeroes_count++] = lt[id];
          ind[id + 1]++;

          non_zeroes[non_zeroes_count] = -tau_div_two_hy * v_plus_fabs_v (v_2[id]) * H_T - tau_mu_div_sq_hy;
          cols[non_zeroes_count++] = id - 1;
          ind[id + 1]++;

          non_zeroes[non_zeroes_count] = +H_T * (1. + tau_div_hx * fabs (v_1[id]) + tau_div_hy * fabs (v_2[id])) + 2 * (four_tau_mu_div_three_sq_hx + tau_mu_div_sq_hy);
          cols[non_zeroes_count++] = id;
          ind[id + 1]++;

          non_zeroes[non_zeroes_count] = +tau_div_two_hy * v_minus_fabs_v (v_2[id]) * H_T - tau_mu_div_sq_hy;
          cols[non_zeroes_count++] = id + 1;
          ind[id + 1]++;

          non_zeroes[non_zeroes_count] = +tau_div_two_hx * v_minus_fabs_v (v_1[id]) * H_T - four_tau_mu_div_three_sq_hx;
          cols[non_zeroes_count++] = rt[id];
          ind[id + 1]++;

          b[id] = H_T * (v_1[id] + tau * f1 (id)) - tau_div_hx * p_ro * (H1 (id) - H1 (lt[id]))
                + tau_mu_div_twelve_hx_hy * (v_2[rt[id + 1]] - v_2[rt[id - 1]] - v_2[lt[id + 1]] + v_2[lt[id - 1]]);
        }
    }
}

void build_sle_for_v_2 (int dim, double curr_t, double h_x, double h_y, double tau,
                        int *st, double *x_coord, double *y_coord,
                        double *v_1, double *v_2, double *h, int *lt, int *rt,
                        double *non_zeroes, int *cols, int *ind, double *b)
{
  double tau_div_hx = tau / h_x;
  double tau_div_hy = tau / h_y;
  double tau_div_two_hx = tau / (2. * h_x);
  double tau_div_two_hy = tau / (2. * h_y);
  double tau_mu_div_sq_hx = tau * mu / (h_x * h_x);
  double four_tau_mu_div_three_sq_hy = (four_mu_div_three * tau) / (h_y * h_y);
  double tau_mu_div_twelve_hx_hy = tau * mu / (12. * h_x * h_y);

  auto inner = [&] (int id) -> bool { return st[id] == 0; };

  auto H2 = [&] (int id) -> double { return (h[id] + h[lt[id]]) / 2.; };
  auto HT = [&] (int id) -> double { return (h[id] + h[id - 1] + h[lt[id]] + h[lt[id - 1]]) / 4.; };
  auto f2 = [&] (int id) -> double { return f_2 (curr_t, x_coord[id], y_coord[id]); };

  auto v_plus_fabs_v  = [&] (double v) -> double {return v + fabs (v); };
  auto v_minus_fabs_v = [&] (double v) -> double {return v - fabs (v); };

  int non_zeroes_count = 0;
  ind[0] = 0;

  for (int id = 0; id < dim; id++)
    {
      if (!inner (id) || fabs (HT (id)) < min_cmp)
        {
          non_zeroes[non_zeroes_count] = 1.;
          cols[non_zeroes_count++] = id;
          ind[id + 1] = ind[id] + 1;
          b[id] = 0.;
        }
      else
        {
          ind[id + 1] = ind[id];

          double H_T = HT (id);
          non_zeroes[non_zeroes_count] = -tau_div_two_hx * v_plus_fabs_v (v_1[id]) * H_T - tau_mu_div_sq_hx;
          cols[non_zeroes_count++] = lt[id];
          ind[id + 1]++;

          non_zeroes[non_zeroes_count] = -tau_div_two_hy * v_plus_fabs_v (v_2[id]) * H_T - four_tau_mu_div_three_sq_hy;
          cols[non_zeroes_count++] = id - 1;
          ind[id + 1]++;

          non_zeroes[non_zeroes_count] = +H_T * (1. + tau_div_hx * fabs (v_1[id]) + tau_div_hy * fabs (v_2[id])) + 2 * (tau_mu_div_sq_hx + four_tau_mu_div_three_sq_hy);
          cols[non_zeroes_count++] = id;
          ind[id + 1]++;

          non_zeroes[non_zeroes_count] = +tau_div_two_hy * v_minus_fabs_v (v_2[id]) * H_T - four_tau_mu_div_three_sq_hy;
          cols[non_zeroes_count++] = id + 1;
          ind[id + 1]++;

          non_zeroes[non_zeroes_count] = +tau_div_two_hx * v_minus_fabs_v (v_1[id]) * H_T - tau_mu_div_sq_hx;
          cols[non_zeroes_count++] = rt[id];
          ind[id + 1]++;

          b[id] = H_T * (v_2[id] + tau * f2 (id)) - tau_div_hy * p_ro * (H2 (id) - H2 (id - 1))
                + tau_mu_div_twelve_hx_hy * (v_1[rt[id + 1]] - v_1[rt[id - 1]] - v_1[lt[id + 1]] + v_1[lt[id - 1]]);
        }
    }
}


void build_sle_for_h (int dim, double curr_t, double h_x, double h_y, double tau,
                      bool *is_rho_valid, int *st, double *x_coord, double *y_coord,
                      double *v_1, double *v_2, double *h, int *lt, int *rt,
                      double *non_zeroes, int *cols, int *ind, double *b)
{
  double half_hx = h_x / 2.;
  double half_hy = h_y / 2.;
  double tau_div_two_hx = tau / (2. * h_x);
  double tau_div_two_hy = tau / (2. * h_y);

  auto not_left  = [&] (int id) -> bool { return st[id] != 3 && st[id] != 5 && st[id] != 7 && st[id] != 14 && st[id] != 16 && st[id] != 17; };
  auto not_down  = [&] (int id) -> bool { return st[id] != 1 && st[id] != 5 && st[id] != 6 && st[id] != 13 && st[id] != 14 && st[id] != 18; };
  auto not_up    = [&] (int id) -> bool { return is_rho_valid[id + 1]; };
  auto not_right = [&] (int id) -> bool { return is_rho_valid[rt[id]]; };

  auto VT1 = [&] (int id) -> double { return (v_1[id] + v_1[id + 1]) / 2.; };
  auto VT2 = [&] (int id) -> double { return (v_2[id] + v_2[rt[id]]) / 2.; };

  auto f0  = [&] (int id) -> double { return f_0 (curr_t, x_coord[id] + half_hx, y_coord[id] + half_hy); };

  auto v_plus_fabs_v  = [&] (double v) -> double {return v + fabs (v); };
  auto v_minus_fabs_v = [&] (double v) -> double {return v - fabs (v); };

  int non_zeroes_count = 0;
  ind[0] = 0;

  for (int id = 0; id < dim; id++)
    {
      if (!is_rho_valid[id])
        {
          non_zeroes[non_zeroes_count] = 1.;
          cols[non_zeroes_count++] = id;
          ind[id + 1] = ind[id] + 1;
          b[id] = 0.;
        }
      else
        {
          ind[id + 1] = ind[id];

          if (not_left (id))
            {
              non_zeroes[non_zeroes_count] = -tau_div_two_hx * v_plus_fabs_v (VT1 (id));
              cols[non_zeroes_count++] = lt[id];
              ind[id + 1]++;
            }

          if (not_down (id))
            {
              non_zeroes[non_zeroes_count] = -tau_div_two_hy * v_plus_fabs_v (VT2 (id));
              cols[non_zeroes_count++] = id - 1;
              ind[id + 1]++;
            }

          non_zeroes[non_zeroes_count] = 1.
                                       + tau_div_two_hx * (v_plus_fabs_v (VT1 (rt[id])) - v_minus_fabs_v (VT1 (id)))
                                       + tau_div_two_hy * (v_plus_fabs_v (VT2 (id + 1)) - v_minus_fabs_v (VT2 (id)));

          cols[non_zeroes_count++] = id;
          ind[id + 1]++;

          if (not_up (id))
            {
              non_zeroes[non_zeroes_count] = +tau_div_two_hy * v_minus_fabs_v (VT2 (id + 1));
              cols[non_zeroes_count++] = id + 1;
              ind[id + 1]++;
            }

          if (not_right (id))
            {
              non_zeroes[non_zeroes_count] = +tau_div_two_hx * v_minus_fabs_v (VT1 (rt[id]));
              cols[non_zeroes_count++] = rt[id];
              ind[id + 1]++;
            }

          b[id] = h[id] + tau * f0 (id);
        }
    }
}

void find_residual (int dim, double *h, bool *is_rho_valid, int *st, double *v_1, double *v_2, double *x_coord, double *y_coord,
                    int *lt, int *rt, double h_x, double h_y, double curr_t, double &norm_c, double &norm_l, double &norm_w)
{
  //std::cout << "final time = " << curr_t << std::endl;

  unique_ptr<double []> residual_h (new double[dim]);
  unique_ptr<double []> residual_v_1 (new double[dim]);
  unique_ptr<double []> residual_v_2 (new double[dim]);

  unique_ptr<double []> real_h (new double[dim]);
  unique_ptr<double []> real_v_1 (new double[dim]);
  unique_ptr<double []> real_v_2 (new double[dim]);

  for (int i = 0; i < dim; i++)
    {
      real_v_1[i] = u_1 (curr_t, x_coord[i], y_coord[i]);
      real_v_2[i] = u_2 (curr_t, x_coord[i], y_coord[i]);
      real_h[i] = is_rho_valid[i] ? ro (curr_t, x_coord[i] + h_x / 2., y_coord[i] + h_y / 2.) : 0.;

      residual_h[i] = is_rho_valid[i] ? h[i] - real_h[i] : 0.;
      residual_v_1[i] = v_1[i] - real_v_1[i];
      residual_v_2[i] = v_2[i] - real_v_2[i];
    }

#if 0
  print_vector (dim, h, true);
  print_vector (dim, real_h.get (), true);
#endif

#if 0
  print_vector (dim, v_1, true);
  print_vector (dim, real_v_1.get (), true);
#endif

#if 0
  print_vector (dim, v_2, true);
  print_vector (dim, real_v_2.get (), true);
#endif

  double tmp;
  /*norm_c = norm_l = norm_w = 0.;
  for (int id = 0; id < dim; id++)
    {
      norm_c = std::max (norm_c, fabs (residual_h[id]));

      tmp = h_x * h_y * residual_h[id] * residual_h[id];
      norm_l += st[id] == 0 ? tmp : tmp / 2.;

      if (st[id] == 0)
        norm_w += (residual_h[id + 1] - residual_h[id]) * (residual_h[id + 1] - residual_h[id]) / h_x
                + (residual_h[rt[id]] - residual_h[id]) * (residual_h[rt[id]] - residual_h[id]) / h_y;
    }*/

  norm_c = norm_l = norm_w = 0.;
  for (int id = 0; id < dim; id++)
    {
      norm_c = std::max (norm_c, fabs (residual_v_1[id]));

      tmp = h_x * h_y * residual_v_1[id] * residual_v_1[id];
      norm_l += st[id] == 0 ? tmp : tmp / 2.;

      if (st[id] == 0)
        norm_w += (residual_v_1[id + 1] - residual_v_1[id]) * (residual_v_1[id + 1] - residual_v_1[id]) / h_x
                + (residual_v_1[rt[id]] - residual_v_1[id]) * (residual_v_1[rt[id]] - residual_v_1[id]) / h_y;

    }

  /*norm_c = norm_l = norm_w = 0.;
  for (int id = 0; id < dim; id++)
    {
      norm_c = std::max (norm_c, fabs (residual_v_2[id]));

      tmp = h_x * h_y * residual_v_2[id] * residual_v_2[id];
      norm_l += st[id] == 0 ? tmp : tmp / 2.;

      if (st[id] == 0)
        norm_w += (residual_v_2[id + 1] - residual_v_2[id]) * (residual_v_2[id + 1] - residual_v_2[id]) / h_x
                + (residual_v_2[rt[id]] - residual_v_2[id]) * (residual_v_2[rt[id]] - residual_v_2[id]) / h_y;
    }*/

  norm_w = sqrt (norm_l + sqrt (norm_w));
  norm_l = sqrt (norm_l);

  norm_c /= rhs_coeff;
  norm_l /= rhs_coeff;
  norm_w /= rhs_coeff;
  return;
}


bool run_programm (int m_x, int m_y, int n, double &norm_c, double &norm_l, double &norm_w)
{
  double h_x = len_x / m_x;
  double h_y = len_y / m_y;
  double tau = len_t / n;

  //std::cout << "len_x, len_y = " << len_x << ", " << len_y << std::endl;
  //std::cout << "h1, h2 = " << h_x << ", " << h_y << std::endl;
  //std::cout << "tau = " << tau << std::endl;

  int nodes_before_hole = (m_y + 1) * (m_x / 3 + 1);
  int nodes_near_hole   = (m_y / 3 + 1) * (m_x / 3 - 1);
  int nodes_after_hole  = (m_y + 1) * (m_x / 3 + 1);

  int dim = nodes_before_hole + nodes_near_hole + nodes_after_hole;

  unique_ptr<bool []> is_rho_valid (new bool[dim]);
  unique_ptr<int []> st (new int[dim]);
  unique_ptr<int []> lt (new int[dim]);
  unique_ptr<int []> rt (new int[dim]);
  unique_ptr<double []> x_coord (new double[dim]);
  unique_ptr<double []> y_coord (new double[dim]);

  unique_ptr<double []> h (new double[dim]);
  unique_ptr<double []> copy_h (new double[dim]);
  unique_ptr<double []> v_1 (new double[dim]);
  unique_ptr<double []> copy_v_1 (new double[dim]);
  unique_ptr<double []> v_2 (new double[dim]);
  unique_ptr<double []> copy_v_2 (new double[dim]);

  std::vector<double> non_zeroes (5 * dim);
  std::vector<int> cols (5 * dim);
  std::vector<int> ind (dim + 1);
  unique_ptr<double []> b (new double[dim]);

  fill_grid (m_x, m_y, h_x, h_y, is_rho_valid.get (), st.get (),
             lt.get (), rt.get (), x_coord.get (), y_coord.get ());

#if 0
  check_grid (dim, is_rho_valid.get (), st.get (), lt.get (), rt.get (), x_coord.get (), y_coord.get ());
#endif

  fill_zero_layer (dim, h_x, h_y, v_1.get (), v_2.get (), is_rho_valid.get (), h.get (), x_coord.get (), y_coord.get ());

#if 0
  print_vector (dim, h.get (), true);
#endif
#if 0
  print_vector (dim, v_1.get (), true);
#endif
#if 0
  print_vector (dim, v_2.get (), true);
#endif

  double curr_t = tau;
  for (int iter = 0; curr_t <= len_t + precision; iter++, curr_t += tau)
    {
      copy_vector (dim, v_1.get (), copy_v_1.get ());
      copy_vector (dim, v_2.get (), copy_v_2.get ());
      copy_vector (dim, h.get (), copy_h.get ());

      //--------------begin_v_1------------------------------------
      non_zeroes.clear (); cols.clear (); ind.clear ();
      build_sle_for_v_1 (dim, curr_t, h_x, h_y, tau, st.get (), x_coord.get (), y_coord.get (),
                         copy_v_1.get (), copy_v_2.get (), h.get (), lt.get (), rt.get (),
                         non_zeroes.data (), cols.data (), ind.data (), b.get ());
#if 0
      print_csr_matr (dim, non_zeroes.data (), cols.data (), ind.data ());
      print_vector (dim, b.get (), true);
#endif

      if (cgs_solve_csr_system (dim, non_zeroes.data (), cols.data (), ind.data (), b.get (), v_1.get ()))
        return true;
      //--------------end_v_1------------------------------------

      //--------------begin_v_2------------------------------------
      non_zeroes.clear (); cols.clear (); ind.clear ();
      build_sle_for_v_2 (dim, curr_t, h_x, h_y, tau, st.get (), x_coord.get (), y_coord.get (),
                         copy_v_1.get (), copy_v_2.get (), h.get (), lt.get (), rt.get (),
                         non_zeroes.data (), cols.data (), ind.data (), b.get ());
#if 0
      print_csr_matr (dim, non_zeroes.data (), cols.data (), ind.data ());
      print_vector (dim, b.get (), true);
#endif

      if (cgs_solve_csr_system (dim, non_zeroes.data (), cols.data (), ind.data (), b.get (), v_2.get ()))
        return true;
      //--------------end_v_2------------------------------------

      //--------------begin_h------------------------------------
      non_zeroes.clear (); cols.clear (); ind.clear ();
      build_sle_for_h (dim, curr_t, h_x, h_y, tau, is_rho_valid.get (), st.get (), x_coord.get (), y_coord.get (),
                       v_1.get (), v_2.get (), h.get (), lt.get (), rt.get (),
                       non_zeroes.data (), cols.data (), ind.data (), b.get ());

#if 0
      print_csr_matr (dim, non_zeroes.data (), cols.data (), ind.data ());
      print_vector (dim, b.get (), true);
#endif
      if (cgs_solve_csr_system (dim, non_zeroes.data (), cols.data (), ind.data (), b.get (), h.get ()))
        return true;
      //--------------end_h-----------------------------------

     //printf ("time: %f\n", curr_t);
    }

  find_residual (dim, h.get (), is_rho_valid.get (), st.get (), v_1.get (), v_2.get (),
                 x_coord.get (), y_coord.get (), lt.get (), rt.get (), h_x, h_y, len_t,
                 norm_c, norm_l, norm_w);

  return false;
}


int main ()
{
  int range_N[4] {20, 40, 80, 160};
  int range_M[4] {30, 60, 120, 240};
  double norm_c, norm_l, norm_w;

  for (auto M : range_M)
    {
      if (M % 3 != 0)
        {
          M = (M % 3 == 1) ? M + 2 : M + 1;
          std::cout << "M must be divided by 3\n" << "Rounding to the nearest larger divisible by 3: " << M << std::endl;
        }
    }

  printf ("$\\tau$ / h & 0.05 & 0.025 & 0.0125 & 0.00625 \\\\\n\\hline\n");
  for (auto N : range_N)
    {
      printf ("%f", 1. / N);
      for (auto M : range_M)
        {
          norm_c = norm_l = norm_w = 0.;
          if (run_programm (M, M, N, norm_c, norm_l, norm_w))
            printf (" & nan ");
          else
            printf (" & %e ", norm_c);
        }
      printf ("\\\\\n");
      for (auto M : range_M)
        {
          norm_c = norm_l = norm_w = 0.;
          if (run_programm (M, M, N, norm_c, norm_l, norm_w))
            printf (" & nan ");
          else
            printf (" & %e ", norm_l);
        }
      printf ("\\\\\n");
      for (auto M : range_M)
        {
          norm_c = norm_l = norm_w = 0.;
          if (run_programm (M, M, N, norm_c, norm_l, norm_w))
            printf (" & nan ");
          else
            printf (" & %e ", norm_w);
        }

      printf ("\\\\\n\\hline\n");
    }

  /*run_programm (60, 60, 40, norm_c, norm_l, norm_w);
  printf ("C: %e \n", norm_c);
  printf ("L: %e \n", norm_l);
  printf ("W: %e \n", norm_w);*/

  return 0;
}

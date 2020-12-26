#include <iostream>
#include <math.h>
#include <memory>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

#define FIX_UNUSED(X) ((void) X)

using std::unique_ptr;

constexpr int max_print = 20;
constexpr double len_t = 1.;
constexpr double len_x = 3.;
constexpr double len_y = 3.;
constexpr double p_ro = 1.;
constexpr double mu = 0.1;
constexpr double omega = 0.5;
constexpr double precision = 1e-9;
constexpr double min_cmp = 1e-12;
constexpr int max_iter = 10000;

double ro (double t, double x_1, double x_2)
{
  return exp (t) * (cos (2 * M_PI * x_1) + 1.5) * (sin (2 * M_PI * x_2) + 1.5);
}

double ro_t (double t, double x_1, double x_2)
{
  return ro (t, x_1, x_2);
}

double ro_x (double t, double x_1, double x_2)
{
  return -2 * M_PI * exp (t) * sin (2 * M_PI * x_1) * (sin (2 * M_PI * x_2) + 1.5);
}

double ro_y (double t, double x_1, double x_2)
{
  return +2 * M_PI * exp (t) * (cos (2 * M_PI * x_1) + 1.5) * cos (2 * M_PI * x_2);
}

double u_1 (double t, double x_1, double x_2)
{
  return exp (t) * sin (2 * M_PI * x_1) * sin (2 * M_PI * x_2);
}

double u_1_t (double t, double x_1, double x_2)
{
  return u_1 (t, x_1, x_2);
}

double u_1_x (double t, double x_1, double x_2)
{
  return 2 * M_PI * exp (t) * cos (2 * M_PI * x_1) * sin (2 * M_PI * x_2);
}

double u_1_y (double t, double x_1, double x_2)
{
  return 2 * M_PI * exp (t) * sin (2 * M_PI * x_1) * cos (2 * M_PI * x_2);
}

double u_1_xx (double t, double x_1, double x_2)
{
  return -4 * M_PI * M_PI * u_1 (t, x_1, x_2);
}

double u_1_xy (double t, double x_1, double x_2)
{
  return 4 * M_PI * M_PI * exp (t) * cos (2 * M_PI * x_1) * cos (2 * M_PI * x_2);
}

double u_1_yy (double t, double x_1, double x_2)
{
  return -4 * M_PI * M_PI * u_1 (t, x_1, x_2);
}

double u_2 (double t, double x_1, double x_2)
{
  return exp (-t) * sin (2 * M_PI * x_1) * sin (2 * M_PI * x_2);
}

double u_2_t (double t, double x_1, double x_2)
{
  return -u_2 (t, x_1, x_2);
}

double u_2_x (double t, double x_1, double x_2)
{
  return 2 * M_PI * exp (-t) * cos (2 * M_PI * x_1) * sin (2 * M_PI * x_2);
}

double u_2_y (double t, double x_1, double x_2)
{
  return 2 * M_PI * exp (-t) * sin (2 * M_PI * x_1) * cos (2 * M_PI * x_2);
}

double u_2_xx (double t, double x_1, double x_2)
{
  return -4 * M_PI * M_PI * u_2 (t, x_1, x_2);
}

double u_2_xy (double t, double x_1, double x_2)
{
  return 4 * M_PI * M_PI * exp (-t) * cos (2 * M_PI * x_1) * cos (2 * M_PI * x_2);
}

double u_2_yy (double t, double x_1, double x_2)
{
  return -4 * M_PI * M_PI * u_2 (t, x_1, x_2);
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
          mu * u_1_yy (t, x_1, x_2) - 4. * mu * u_1_xx (t, x_1, x_2) / 3. - mu * u_2_xy (t, x_1, x_2) / 3.) / rho;
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
          mu * u_2_xx (t, x_1, x_2) - 4. * mu * u_2_yy (t, x_1, x_2) / 3. - mu * u_1_xy (t, x_1, x_2) / 3.) / rho;
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

void print_vector (int size, double *v, bool if_debug = false)
{
  int print_size = if_debug ? size : std::min (size, max_print);
  for (int i = 0; i < print_size; i++)
    {
      if (i % static_cast<int> (sqrt (size)) == 0)
        printf ("\n");
      printf ("%7.3f ", v[i]);
    }

  printf("\n");
}


void print_matrix (int size, double *matrix, double *b)
{
  for (int i = 0; i < size; i++)
    {
      for (int j = 0; j < size; j++)
        {
          if (fabs (matrix[i * size + j]) > min_cmp)
            printf ("mm[%d][%d] = %.3f ", i, j, matrix[i * size + j]);
        }
      printf ("b[%d] = %.3f\n", i, b[i]);
    }

  printf("\n");
}

int residual_gauss (int size, double *matrix, double *rhs, double *sol)
{
  for (int i = 0; i < size; i++)
    {
      double tmp = 0.;
      for (int j = 0; j < size; j++)
        tmp += matrix[i * size + j] * sol[j];
      rhs[i] -= tmp;
    }

#if 1
  print_vector (size, rhs, true);
#endif
}


void solve_gauss (int size, double *matrix, double *rhs, double *sol)
{
#if 0
  unique_ptr<double []> matrix_copy (new double[size * size]);
  unique_ptr<double []> rhs_copy (new double[size]);
  copy_vector (size * size, matrix, matrix_copy.get ());
  copy_vector (size, rhs, rhs_copy.get ());
#endif

#if 0
  print_matrix (size, matrix, rhs);
#endif

  for (int i = 0; i < size; i++)
    {
      double tmp = matrix[i * size + i];
      for (int j = i; j < size; j++)
        matrix[i * size + j] /= tmp;
      rhs[i] /= tmp;

      for (int j = i + 1; j < size; j++)
        {
          for (int k = i + 1; k < size; k++)
            matrix[j * size + k] -= matrix[i * size + k] * matrix[j * size + i];

          rhs[j] -= rhs[i] * matrix[j * size + i];
        }
    }

  for (int i = size - 1; i >= 0; i--)
    {
      for (int j = i + 1; j < size; j++)
        rhs[i] -= rhs[j] * matrix[i * size + j];
    }

  copy_vector (size, rhs, sol);

#if 0
  residual_gauss (size, matrix_copy.get (), rhs_copy.get (), sol);
#endif
}

void fill_grid_rect (int m_x, int m_y, double h_x, double h_y, bool *is_rho_valid,
                     int *st, int *lt, int *rt, double *x_coord, double *y_coord)
{
  int dim = (m_x + 1) * (m_y + 1);

  for (int curr_node = 0; curr_node < dim; curr_node++)
    {
      int row = curr_node % (m_y + 1);
      int col = curr_node / (m_y + 1);

      x_coord[curr_node] = h_x * col;
      y_coord[curr_node] = h_y * row;

      if (col == 0) // done
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
              st[curr_node] = 3;
            }
        }
      else if (col == m_x) // done
        {
          lt[curr_node] = curr_node - (m_y + 1);
          rt[curr_node] = -1;
          is_rho_valid[curr_node] = false;

          if (row == 0)
            {
              st[curr_node] = 6;
            }
          else if (row == m_y)
            {
              st[curr_node] = 8;
            }
          else
            {
              st[curr_node] = 4;
            }
        }
      else // done
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
                        const int *st, const double *x_coord, const double *y_coord,
                        const double *v_1, const double *v_2, const double *h,
                        const int *lt, const int *rt,
                        double *matrix, double *b)
{
  double two_hx = 2. * h_x;
  double two_hy = 2. * h_y;
  double four_mu_div_three_sq_hx = (4. * mu) / (3. * h_x * h_x);
  double mu_div_sq_hy = mu / (h_y * h_y);
  double mu_div_twelve_hx_hy = mu / (12. * h_x * h_y);

  auto inner = [&] (int id) -> bool { return st[id] == 0; };
  auto not_inner = [&] (int id) -> bool { return st[id] != 0; };

  auto H1 = [&] (int id) -> double { return (h[id] + h[id - 1]) / 2.; };
  auto HT = [&] (int id) -> double { return (h[id] + h[id - 1] + h[lt[id]] + h[lt[id - 1]]) / 4.; };
  auto f1 = [&] (int id) -> double { return f_1 (curr_t, x_coord[id], y_coord[id]); };

  for (int id = 0; id < dim; id++)
    {
      if (not_inner (id) || fabs (HT (id)) < min_cmp)
        {
          matrix[id * dim + id] = 1.;
          b[id] = 0.;
        }
      else
        {
          // diag
          matrix[id * dim + id] = +HT (id) * (1. / tau + fabs (v_1[id]) / h_x + fabs (v_2[id]) / h_y) + 2 * (four_mu_div_three_sq_hx + mu_div_sq_hy);

          // offdiag (left, down, up, right)
          if (inner (lt[id]))
            matrix[id * dim + lt[id]] = -HT (id) * (v_1[id] + fabs (v_1[id])) / two_hx - four_mu_div_three_sq_hx;
          if (inner (id - 1))
            matrix[id * dim + id - 1] = -HT (id) * (v_2[id] + fabs (v_2[id])) / two_hy - mu_div_sq_hy;
          if (inner (id + 1))
            matrix[id * dim + id + 1] = +HT (id) * (v_2[id] - fabs (v_2[id])) / two_hy - mu_div_sq_hy;
          if (inner (rt[id]))
            matrix[id * dim + rt[id]] = +HT (id) * (v_1[id] - fabs (v_1[id])) / two_hx - four_mu_div_three_sq_hx;

          // rhs
          b[id] = HT (id) * (v_1[id] / tau + f1 (id)) - p_ro * (H1 (id) - H1 (lt[id])) / h_x
                + mu_div_twelve_hx_hy * (v_2[rt[id + 1]] - v_2[rt[id - 1]] - v_2[lt[id + 1]] + v_2[lt[id - 1]]);
        }
    }
}

void build_sle_for_v_2 (int dim, double curr_t, double h_x, double h_y, double tau,
                        const int *st, const double *x_coord, const double *y_coord,
                        const double *v_1, const double *v_2, const double *h,
                        const int *lt, const int *rt,
                        double *matrix, double *b)
{
  double two_hx = 2. * h_x;
  double two_hy = 2. * h_y;
  double mu_div_sq_hx = mu / (h_x * h_x);
  double four_mu_div_three_sq_hy = (4. * mu) / (3. * h_y * h_y);
  double mu_div_twelve_hx_hy = mu / (12. * h_x * h_y);

  auto inner = [&] (int id) -> bool { return st[id] == 0; };
  auto not_inner = [&] (int id) -> bool { return st[id] != 0; };

  auto H2 = [&] (int id) -> double { return (h[id] + h[lt[id]]) / 2.; };
  auto HT = [&] (int id) -> double { return (h[id] + h[id - 1] + h[lt[id]] + h[lt[id - 1]]) / 4.; };
  auto f2 = [&] (int id) -> double { return f_2 (curr_t, x_coord[id], y_coord[id]); };

  for (int id = 0; id < dim; id++)
    {
      if (not_inner (id) || fabs (HT (id)) < min_cmp)
        {
          matrix[id * dim + id] = 1.;
          b[id] = 0.;
        }
      else
        {
          // diag
          matrix[id * dim + id] = HT (id) * (1. / tau + fabs (v_1[id]) / h_x + fabs (v_2[id]) / h_y) + 2 * (mu_div_sq_hx + four_mu_div_three_sq_hy);

          // offdiag (left, down, up, right)
          if (inner (lt[id]))
            matrix[id * dim + lt[id]] = -HT (id) * (v_1[id] + fabs (v_1[id])) / two_hx - mu_div_sq_hx;
          if (inner (id - 1))
            matrix[id * dim + id - 1] = -HT (id) * (v_2[id] + fabs (v_2[id])) / two_hy - four_mu_div_three_sq_hy;
          if (inner (id + 1))
            matrix[id * dim + id + 1] = +HT (id) * (v_2[id] - fabs (v_2[id])) / two_hy - four_mu_div_three_sq_hy;
          if (inner (rt[id]))
            matrix[id * dim + rt[id]] = +HT (id) * (v_1[id] - fabs (v_1[id])) / two_hx - mu_div_sq_hx;

          // rhs
          b[id] = HT (id) * (v_2[id] / tau + f2 (id)) - p_ro * (H2 (id) - H2 (id - 1)) / h_y
                + mu_div_twelve_hx_hy * (v_1[rt[id + 1]] - v_1[rt[id - 1]] - v_1[lt[id + 1]] + v_1[lt[id - 1]]);
        }
    }
}

void build_sle_for_h (int dim, double curr_t, double h_x, double h_y, double tau,
                      const bool *is_rho_valid, const int *st, const double *x_coord, const double *y_coord,
                      const double *v_1, const double *v_2, const double *h,
                      const int *lt, const int *rt,
                      double *matrix, double *b)
{
  double half_hx = h_x / 2.;
  double half_hy = h_y / 2.;
  double two_hx = 2. * h_x;
  double two_hy = 2. * h_y;

  auto not_left  = [&] (int id) -> bool { return fabs (x_coord[id]) > min_cmp; };
  auto not_down  = [&] (int id) -> bool { return fabs (y_coord[id]) > min_cmp; };
  auto not_up    = [&] (int id) -> bool { return is_rho_valid[id + 1]; };
  auto not_right = [&] (int id) -> bool { return is_rho_valid[rt[id]]; };

  auto VT1 = [&] (int id) -> double { return (v_1[id] + v_1[id + 1]) / 2.; };
  auto VT2 = [&] (int id) -> double { return (v_2[id] + v_2[rt[id]]) / 2.; };

  auto f0  = [&] (int id) -> double { return f_0 (curr_t, x_coord[id] + half_hx, y_coord[id] + half_hy); };

  for (int id = 0; id < dim; id++)
    {
      if (!is_rho_valid[id])
        {
          matrix[id * dim + id] = 1.;
          b[id] = 0.;
        }
      else
        {
          // diag
          matrix[id * dim + id] = 1. / tau
                                   + (VT1 (rt[id]) + fabs (VT1 (rt[id])) - VT1 (id) + fabs (VT1 (id))) / two_hx
                                   + (VT2 (id + 1) + fabs (VT2 (id + 1)) - VT2 (id) + fabs (VT2 (id))) / two_hy;

          // offdiag (left, down, up, right)
          if (not_left (id))
            matrix[id * dim + lt[id]] = -(VT1 (id) + fabs (VT1 (id))) / two_hx;
          if (not_down (id))
            matrix[id * dim + id - 1] = -(VT2 (id) + fabs (VT2 (id))) / two_hy;
          if (not_up (id))
            matrix[id * dim + id + 1] = +(VT2 (id + 1) - fabs (VT2 (id + 1))) / two_hy;
          if (not_right (id))
            matrix[id * dim + rt[id]] = +(VT1 (rt[id]) - fabs (VT1 (rt[id]))) / two_hx;

          // rhs
          b[id] = h[id] / tau + f0 (id);
        }
    }
}

void find_residual (int dim, double *h, bool *is_rho_valid, int *st, double *v_1, double *v_2, double *x_coord, double *y_coord,
                    int *lt, int *rt, double h_x, double h_y)
{
  unique_ptr<double []> residual_h (new double[dim]);
  unique_ptr<double []> residual_v_1 (new double[dim]);
  unique_ptr<double []> residual_v_2 (new double[dim]);

  unique_ptr<double []> real_h (new double[dim]);
  unique_ptr<double []> real_v_1 (new double[dim]);
  unique_ptr<double []> real_v_2 (new double[dim]);

  for (int i = 0; i < dim; i++)
    {
      real_h[i] = is_rho_valid[i] ? ro (len_t, x_coord[i] + h_x / 2., y_coord[i] + h_y / 2.) : 0.;
      real_v_1[i] = u_1 (len_t, x_coord[i], y_coord[i]);
      real_v_2[i] = u_2 (len_t, x_coord[i], y_coord[i]);

      residual_h[i] = is_rho_valid[i] ? h[i] - real_h[i] : 0.;
      residual_v_1[i] = v_1[i] - real_v_1[i];
      residual_v_2[i] = v_2[i] - real_v_2[i];
    }

#if 1
  print_vector (dim, h, true);
  print_vector (dim, real_h.get (), true);
#endif

#if 1
  print_vector (dim, v_1, true);
  print_vector (dim, real_v_1.get (), true);
#endif

#if 1
  print_vector (dim, v_2, true);
  print_vector (dim, real_v_2.get (), true);
#endif

  double tmp;
  double norm_c = 0., norm_l = 0., norm_w = 0.;

  for (int id = 0; id < dim; id++)
    {
      norm_c = std::max (norm_c, fabs (residual_h[id]));

      tmp = h_x * h_y * residual_h[id] * residual_h[id];
      norm_l += st[id] == 0 ? tmp : tmp / 2.;

      if (st[id] == 0)
        norm_w += (residual_h[id + 1] - residual_h[id]) * (residual_h[id + 1] - residual_h[id]) / h_x
                + (residual_h[rt[id]] - residual_h[id]) * (residual_h[rt[id]] - residual_h[id]) / h_y;
    }

  printf ("H: C_Norma: %e\n", norm_c);
  printf ("H: L_Norma: %e\n", sqrt (norm_l));
  printf ("H: W_Norma: %e\n\n", sqrt (norm_l + sqrt (norm_w)));

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

  printf ("V_1: C_Norma: %e\n", norm_c);
  printf ("V_1: L_Norma: %e\n", sqrt (norm_l));
  printf ("V_1: W_Norma: %e\n\n", sqrt (norm_l + sqrt (norm_w)));

  norm_c = norm_l = norm_w = 0.;
  for (int id = 0; id < dim; id++)
    {
      norm_c = std::max (norm_c, fabs (residual_v_2[id]));

      tmp = h_x * h_y * residual_v_2[id] * residual_v_2[id];
      norm_l += st[id] == 0 ? tmp : tmp / 2.;

      if (st[id] == 0)
        norm_w += (residual_v_2[id + 1] - residual_v_2[id]) * (residual_v_2[id + 1] - residual_v_2[id]) / h_x
                + (residual_v_2[rt[id]] - residual_v_2[id]) * (residual_v_2[rt[id]] - residual_v_2[id]) / h_y;
    }

  printf ("V_2: C_Norma: %e\n", norm_c);
  printf ("V_2: L_Norma: %e\n", sqrt (norm_l));
  printf ("V_2: W_Norma: %e\n\n", sqrt (norm_l + sqrt (norm_w)));
  return;
}


void run_programm (int m_x, int m_y, int n)
{
  double h_x = len_x / m_x;
  double h_y = len_y / m_y;
  double tau = len_t / n;

  std::cout << "len_x, len_y = " << len_x << ", " << len_y << std::endl;
  std::cout << "h1, h2 = " << h_x << ", " << h_y << std::endl;
  std::cout << "tau = " << tau << std::endl;

  int dim = (m_y + 1) * (m_x + 1);

  unique_ptr<bool []> is_rho_valid (new bool[dim]);
  unique_ptr<int []> st (new int[dim]);
  unique_ptr<int []> lt (new int[dim]);
  unique_ptr<int []> rt (new int[dim]);
  unique_ptr<double []> x_coord (new double[dim]);
  unique_ptr<double []> y_coord (new double[dim]);

  unique_ptr<double []> h (new double[dim]);
  unique_ptr<double []> h_copy (new double[dim]);
  unique_ptr<double []> v_1 (new double[dim]);
  unique_ptr<double []> v_1_copy (new double[dim]);
  unique_ptr<double []> v_2 (new double[dim]);
  unique_ptr<double []> v_2_copy (new double[dim]);

  unique_ptr<double []> matrix (new double[dim * dim]);
  unique_ptr<double []> b (new double[dim]);

  fill_grid_rect (m_x, m_y, h_x, h_y, is_rho_valid.get (), st.get (),
                  lt.get (), rt.get (), x_coord.get (), y_coord.get ());

#if 0
  check_grid (dim, is_rho_valid.get (), st.get (), lt.get (), rt.get (), x_coord.get (), y_coord.get ());
#endif

  fill_zero_layer (dim, h_x, h_y, v_1.get (), v_2.get (), is_rho_valid.get (), h.get (), x_coord.get (), y_coord.get ());

#if 1
  print_vector (dim, h.get (), true);
#endif
#if 1
  print_vector (dim, v_1.get (), true);
#endif
#if 1
  print_vector (dim, v_2.get (), true);
#endif

  double curr_t = tau;
  for (int iter = 0; iter < n; iter++, curr_t += tau)
    {
#if 0
      print_vector (dim, h.get (), true);
      print_vector (dim, v_1.get (), true);
      print_vector (dim, v_2.get (), true);
#endif

      copy_vector (dim, h.get (), h_copy.get ());
      copy_vector (dim, v_1.get (), v_1_copy.get ());
      copy_vector (dim, v_2.get (), v_2_copy.get ());

      //--------------begin_v_1------------------------------------
      set_zero_vector (dim * dim, matrix.get ());
      set_zero_vector (dim, b.get ());
      build_sle_for_v_1 (dim, curr_t, h_x, h_y, tau,
                         st.get (), x_coord.get (), y_coord.get (),
                         v_1_copy.get (), v_2_copy.get (), h_copy.get (),
                         lt.get (), rt.get (),
                         matrix.get (), b.get ());

      solve_gauss (dim, matrix.get (), b.get (), v_1.get ());
      //--------------end_v_1------------------------------------

      //--------------begin_v_2------------------------------------
      set_zero_vector (dim * dim, matrix.get ());
      set_zero_vector (dim, b.get ());
      build_sle_for_v_2 (dim, curr_t, h_x, h_y, tau,
                         st.get (), x_coord.get (), y_coord.get (),
                         v_1_copy.get (), v_2_copy.get (), h_copy.get (),
                         lt.get (), rt.get (),
                         matrix.get (), b.get ());

      solve_gauss (dim, matrix.get (), b.get (), v_2.get ());
      //--------------end_v_2------------------------------------

      //--------------begin_h------------------------------------
      set_zero_vector (dim * dim, matrix.get ());
      set_zero_vector (dim, b.get ());
      build_sle_for_h (dim, curr_t, h_x, h_y, tau,
                       is_rho_valid.get (), st.get (), x_coord.get (), y_coord.get (),
                       v_1.get (), v_2.get (), h_copy.get (),
                       lt.get (), rt.get (),
                       matrix.get (), b.get ());

      solve_gauss (dim, matrix.get (), b.get (), h.get ());
      //--------------end_h------------------------------------

      printf ("iter = %d from %d, curr_t = %f\n", iter + 1, n, curr_t);
    }

  find_residual (dim, h.get (), is_rho_valid.get (), st.get (), v_1.get (), v_2.get (), x_coord.get (), y_coord.get (), lt.get (), rt.get (), h_x, h_y);
}


int main ()
{
/*  int range_N[4] {100, 200, 400, 800};
  int range_M[4] {100, 200, 400, 800};
  double norm_c, norm_l, norm_w;

  printf ("N / M & 100 & 200 & 400 & 800 \\\\\n\\hline\n");
  for (auto N : range_N)
    {
      printf ("%u", N);
      for (auto M : range_M)
        {
          norm_c = norm_l = norm_w = 0.;
          run_programm (M, N, norm_c, norm_l, norm_w);
          printf (" & %e ", norm_c);
        }
      printf ("\\\\\n");
      for (auto M : range_M)
        {
          norm_c = norm_l = norm_w = 0.;
          run_programm (M, N, norm_c, norm_l, norm_w);
          printf (" & %e ", norm_l);
        }
      printf ("\\\\\n");
      for (auto M : range_M)
        {
          norm_c = norm_l = norm_w = 0.;
          run_programm (M, N, norm_c, norm_l, norm_w);
          printf (" & %e ", norm_w);
        }

      printf ("\\\\\n\\hline\n");
    } */


  int M_x, M_y, N;
  std::cout << "Input m_x, m_y and N" << std::endl;
  std::cin >> M_x >> M_y >> N;

  if (M_x % 3 != 0)
    {
      M_x = (M_x % 3 == 1) ? M_x + 2 : M_x + 1;
      std::cout << "M_x must be divided by 3\n" << "Rounding to the nearest larger divisible by 3: " << M_x << std::endl;
    }

  if (M_y % 3 != 0)
    {
      M_y = (M_y % 3 == 1) ? M_y + 2 : M_y + 1;
      std::cout << "M_y must be divided by 3\n" << "Rounding to the nearest larger divisible by 3: " << M_y << std::endl;
    }

  run_programm (M_x, M_y, N);

  return 0;
}


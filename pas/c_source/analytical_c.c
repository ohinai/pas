

#include<stdio.h>
#include<math.h>
#include<string.h>
#include<stdlib.h>

#include <Python.h>


#define PI 3.141592653		

static PyObject *
analytical_pointsource_all_d(PyObject *self, PyObject *args)
{

  /*
    x_in, y_in, z_in, t_in, 
    x_well, y_well, z_well, 
    n_x, n_y, n_z, 
    a, b, d, 
    n_max, m_max, w_max 
   */
    double pressure;
    double n_x, n_y, n_z;
    double x_well, y_well, z_well;
    double a, b, d;
    int x_i, y_i, z_i, t_i;
    double x_in, y_in, z_in;
    int m, n, w;
    double t_in;

    int n_max, m_max, w_max;

    if (!PyArg_ParseTuple(args, "dddddddddddddiii",
			  &x_in, &y_in, &z_in, &t_in, 
			  &x_well, &y_well, &z_well, 
			  &n_x, &n_y, &n_z, 
			  &a, &b, &d, 
			  &n_max, &m_max, &w_max))
        return NULL;

    double current_n; 
    double current_n_m;
    double current_value;

    double pow_pi_n_a, pow_pi_m_b;

    pressure = 0.;

    for (n=1; n<n_max; n++){
      current_n = sin(n*PI*x_in/a);
      current_n *= sin(n*PI*x_well/a);
      pow_pi_n_a = pow(PI*n/a,2)*n_x;
      for (m=1; m<m_max; m++){
	current_n_m = current_n*sin(m*PI*y_in/b);
	current_n_m *= sin(m*PI*y_well/b);
	pow_pi_m_b = pow(PI*m/b,2)*n_y;
	for (w=1; w<w_max; w++){
	  current_value = current_n_m*sin(w*PI*z_in/d);
	  current_value *= sin(w*PI*z_well/d);
	  current_value /= pow_pi_n_a+pow_pi_m_b+pow(PI*w/d,2)*n_z;
	  current_value *= -exp(-t_in*(pow_pi_n_a+pow_pi_m_b+pow(PI*w/d,2)*n_z))+1.;
	  pressure += current_value;
	}
      }
    }
    
    pressure *= 64.;
    
    return Py_BuildValue("d", pressure); 
}


static PyObject *
analytical_pointsource_all_n(PyObject *self, PyObject *args)
{

  /*
    x_in, y_in, z_in, t_in, 
    x_well, y_well, z_well, 
    n_x, n_y, n_z, 
    a, b, d, 
    n_max, m_max, w_max 
   */

    double pressure;
    double n_x, n_y, n_z;
    double x_well, y_well, z_well;
    double a, b, d;
    int x_i, y_i, z_i, t_i;
    double x_in, y_in, z_in;
    int m, n, w;
    double t_in;

    int n_max, m_max, w_max;

    if (!PyArg_ParseTuple(args, "dddddddddddddiii",
			  &x_in, &y_in, &z_in, &t_in, 
			  &x_well, &y_well, &z_well, 
			  &n_x, &n_y, &n_z, 
			  &a, &b, &d, 
			  &n_max, &m_max, &w_max))
        return NULL;

    double current_n; 
    double current_n_m;
    double current_value;

    double pow_pi_n_a, pow_pi_m_b;

    pressure = 0.;

    for (n=0; n<n_max; n++){
      current_n = cos(n*PI*x_in/a);
      current_n *= cos(n*PI*x_well/a);
      if (n==0) 
	current_n *= .5 ;
      pow_pi_n_a = pow(PI*n/a,2)*n_x;
      for (m=0; m<m_max; m++){
	current_n_m = current_n*cos(m*PI*y_in/b);
	current_n_m *= cos(m*PI*y_well/b);
	if (m==0)
	  current_n_m *= .5;
	pow_pi_m_b = pow(PI*m/b,2)*n_y;
	for (w=0; w<w_max; w++){
	  current_value = current_n_m*cos(w*PI*z_in/d);
	  current_value *= cos(w*PI*z_well/d);
	  if (w==0)
	    current_value *= .5;
	  if (n==0 && m == 0 && w == 0)
	    pressure += current_value*t_in;
	  else{
	    current_value /= pow_pi_n_a+pow_pi_m_b+pow(PI*w/d,2)*n_z;
	    current_value *= -exp(-t_in*(pow_pi_n_a+pow_pi_m_b+pow(PI*w/d,2)*n_z))+1.;
	    pressure += current_value;
	  }
	}
      }
    }
    
    pressure *= 64.;
    
    return Py_BuildValue("d", pressure);
}

static PyObject *
analytical_cincoley_sum_1(PyObject *self, PyObject *args)
{

  /* Computes the the first summation in equation 11 in Cinco Ley et al 
     "Transient Pressure Behavior for a Well With Finite-Conductivity 
     Vertical Fracture."
   */
  double t_DK; 
  double x_DJ; 
  double eta_fD;
  int n_max; 
  
  if (!PyArg_ParseTuple(args, "dddi",
			&t_DK, &x_DJ, &eta_fD, &n_max))
    return NULL;  

  long int n; 
  double summation = 0.;
  for (n=1; n<n_max; n++){
    summation += 1./(n*n)*(1.-exp(-PI*PI*n*n*eta_fD*t_DK))*cos(n*PI*x_DJ);
  }
  return Py_BuildValue("d", summation);
}

static PyObject *
analytical_cincoley_sum_2(PyObject *self, PyObject *args)
{

  /* Computes the the second summation in equation 11 in Cinco Ley et al 
     "Transient Pressure Behavior for a Well With Finite-Conductivity 
     Vertical Fracture."
   */
  double delta_t_K_l; 
  double delta_t_K_l_min1; 
  double x_DJ; 
  double x_DI; 
  double eta_fD;
  int N; 
  int n_max; 
  
  if (!PyArg_ParseTuple(args, "dddddii",
			&delta_t_K_l, 
			&delta_t_K_l_min1,
			&x_DJ, 
			&x_DI, 
			&eta_fD, 
			&N, 
			&n_max))
    return NULL;  

  long int n; 
  double summation = 0.;
  double current_entry; 
  for (n=1; n<n_max; n++){
    current_entry = 1./(n*n*n);
    current_entry *= (exp(-PI*PI*n*n*eta_fD*delta_t_K_l)-\
		      exp(-PI*PI*n*n*eta_fD*delta_t_K_l_min1)); 
    current_entry *= cos(n*PI*x_DJ);
    current_entry *= cos(n*PI*x_DI);
    current_entry *= sin(n*PI/(2.*N));
    summation += current_entry; 
  }
  return Py_BuildValue("d", summation);
}



static PyMethodDef AnalyticalMethods[] = {
  {"pointsource_all_d",  analytical_pointsource_all_d, METH_VARARGS,
   "Find solution to Dirichlet problem."},
  {"pointsource_all_n",  analytical_pointsource_all_n, METH_VARARGS,
   "Find solution to Neumann problem."},
  {"cincoley_sum_1",  analytical_cincoley_sum_1, METH_VARARGS,
   "Find first summation in Cinco Ley solution."},
  {"cincoley_sum_2",  analytical_cincoley_sum_2, METH_VARARGS,
   "Find second summation in Cinco Ley solution."},
    
};


PyMODINIT_FUNC
initanalytical(void)
{
  (void) Py_InitModule("analytical", AnalyticalMethods);
}



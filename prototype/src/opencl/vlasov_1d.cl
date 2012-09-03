#define REORDER

// system_size is defined during program building

// solves a bunch of tridiagonal linear systems
// much better performance when doing data reordering before
// so that all memory accesses are coalesced (who-ho!)
__kernel void vlasov_1d( int num_systems,
                         int system_size,
                        __global float *a_d, 
                        __global float *b_d, 
                        __global float *c_d, 
                        __global float *d_d, 
                        __global float *x_d )
{
	
	int i = get_global_id(0);
	
	// need to check for in-bounds because of the thread block size
    if (i >= num_systems) return;

#ifndef REORDER
	int stride = 1;
	int base_idx = i * system_size;
#else
	int stride = num_systems;
	int base_idx = i;
#endif

	// local memory
	float a[128];

	float c1, c2, c3;
	float f_i, x_prev, x_next;
	
	// solving next system:	
	// c1 * u_i+1 + c2 * u_i + c3 * u_i-1 = f_i
	
	c1 = c_d[base_idx];
	c2 = b_d[base_idx];
	f_i = d_d[base_idx];

	a[1] = - c1 / c2;
	x_prev = f_i / c2;

	// forward trace
	int idx = base_idx;
	x_d[base_idx] = x_prev;
	for (int k = 1; k < system_size-1; k++)
	{
		idx += stride;
	
		c1 = c_d[idx];
		c2 = b_d[idx];
		c3 = a_d[idx];
		f_i = d_d[idx];
		
		float q = (c3 * a[k] + c2);
		float t = 1 / q; 
		x_next = (f_i - c3 * x_prev) * t;
		x_d[idx] = x_prev = x_next;
		
		a[k+1] = - c1 * t;
	}
	
	idx += stride;

	c2 = b_d[idx];
	c3 = a_d[idx];
	f_i = d_d[idx];

	float q = (c3 * a[system_size-1] + c2);
	float t = 1 / q; 
	x_next = (f_i - c3 * x_prev) * t;
	x_d[idx] = x_prev = x_next;

	// backward trace
	for (int k = system_size-2; k >= 0; k--)
	{
		idx -= stride;
		x_next = x_d[idx];
		x_next += x_prev * a[k+1];
		x_d[idx] = x_prev = x_next;
	}
}

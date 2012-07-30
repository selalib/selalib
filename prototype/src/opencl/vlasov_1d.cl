__kernel void vlasov_1d(const int size_system, 
                        const int num_systems, 
                        __global float * a, 
                        __global float * b, 
                        __global float * c, 
                        __global float * d, 
                        __global float * x){

  int i = get_global_id(0);

  if(i < size_system) x[i] = a[i]+b[i]+c[i]+d[i];

}

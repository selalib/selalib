In the build directory run

make all
make test_pic_1d

Run with 2**7 mesh cells, and spline degree sdeg
Output is saved in path

Landau damping with delta-f and Runge Kutta 2
mpirun -np 2   ./bin/test_pic_1d nmark=10000,femp=7,sdeg=3,tstepw=0.1,tsteps=100,ppusher=\"rk2\",psolver=\"fem\",lalpha=0.1,path=\'/tmp/\',scenario=\'landau\',deltaf=1

Error estimates with mersons pusher
mpirun -np 2   ./bin/test_pic_1d nmark=10000,femp=7,sdeg=3,tstepw=0.1,tsteps=100,ppusher=\"merson\",psolver=\"fem\",lalpha=0.1,path=\'/tmp/\',scenario=\'landau\',deltaf=1


To view results go to output directory and type

gnuplot -persistent ./pic1dresult.gnu 



#Numerical Testcase #1
Linear Landau damping with decay rate -0.17146  (without Control Variate)
Number of particles 5e6
Spline degree: 3
tstepw= dt = 0.1
tsteps=100;
pusher= verlet
fieldsolver= fem
excitation= lalpha=0.1
deltaf=0  Delta f is turned off
Number of cells = 2^femp= 2^5=32

mpirun -np 2   ./bin/test_pic_1d nmark=500000,femp=5,sdeg=3,tstepw=0.1,tsteps=100,ppusher=\"verlet\",psolver=\"fem\",lalpha=0.1,scenario=\'landau\',deltaf=0

The output should be:
Damping Factor (fit): -0.17146 in [-0.20815,-0.13477] by 0.95%


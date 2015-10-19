In the build directory run

make all
make test_pic_1d

Run with 2**7 mesh cells, and spline degree sdeg
Output is saved in path

Landau damping with delta-f and Runge Kutta 2, and 2^5 elements
mpirun -np 2   ./bin/test_pic_1d nmark=100000,femp=5,sdeg=3,tstepw=0.1,tsteps=100,ppusher=\"rk2\",psolver=\"fem\",lalpha=0.1,path=\'/tmp/\',scenario=\'landau\',deltaf=1

With live output

mpirun -np 2   ./bin/test_pic_1d nmark=100000,femp=5,sdeg=3,tstepw=0.1,tsteps=100,ppusher=\"rk2\",psolver=\"fem\",lalpha=0.1,path=\'/tmp/\',scenario=\'landau\',deltaf=1,gpinline=1 | gnuplot


Error estimates with mersons pusher
mpirun -np 2   ./bin/test_pic_1d nmark=10000,femp=7,sdeg=3,tstepw=0.1,tsteps=100,ppusher=\"merson\",psolver=\"fem\",lalpha=0.1,path=\'/tmp/\',scenario=\'landau\',deltaf=1


To view results go to output directory and type

gnuplot -persistent ./pic1dresult.gnu




#Numerical Testcase #1 (Landau damping)
Linear Landau damping with decay rate -0.17161  (without Control Variate)
Number of particles 5e5
Spline degree: 3
tstepw= dt = 0.1
tsteps=150;
pusher= verlet
fieldsolver= fem
excitation= lalpha=0.1
deltaf=0  Delta f is turned off
Number of cells = 2^femp= 2^5=32

mpirun -np 2   ./bin/test_pic_1d nmark=500000,femp=5,sdeg=3,tstepw=0.1,tsteps=150,ppusher=\"verlet\",psolver=\"fem\",lalpha=0.1,scenario=\'landau\',deltaf=0

The output should be:
Damping Factor (fit): -0.17161 in [-0.19563,-0.14758] by 0.95%

The demo result file is:
demoresult.gnu


#Numerical Testcase #2 (Bump on tail)
%FEM
mpirun -np 2   ./bin/test_pic_1d nmark=500000,femp=6,sdeg=3,tstepw=0.1,tsteps=500,ppusher=\"verlet\",psolver=\"fem\",lalpha=0.04,scenario=\'bump\',deltaf=0

%Fourier
mpirun -np 2   ./bin/test_pic_1d nmark=500000,femp=6,sdeg=3,tstepw=0.1,tsteps=500,ppusher=\"verlet\",psolver=\"fourier\",lalpha=0.04,scenario=\'bump\',deltaf=0

%Finite differences
mpirun -np 2   ./bin/test_pic_1d nmark=500000,femp=6,sdeg=3,tstepw=0.1,tsteps=500,ppusher=\"verlet\",psolver=\"fd\",lalpha=0.04,scenario=\'bump\',deltaf=0

#Numerical Testcase #3 (Twostream)




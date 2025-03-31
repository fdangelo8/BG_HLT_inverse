# BG_HLT_inverse

This is a modified version of the HLT Backus--gilbert code provided by M. Naviglio.

## Table of Contents
- [Installation](#installation)
- [Usage](#usage)
- [Output](#output)
- [Contributors](#contributors)

## Installation <a name="installation"></a>
```bash
# Clone the repository
git clone https://github.com/your-username/project-name.git

# Navigate into the project directory
cd BG_HLT_inverse
```
The following dipendencies are required. You can find them in the `package` directory. In the following, the procedure to correctly install them.
### gmp 
```bash
wget https://ftp.gnu.org/gnu/gmp/gmp-6.3.0.tar.gz
tar -xzvf gmp-6.3.0.tar.gz
cd gmp-6.3.0
./configure --prefix=${HOME}
make
make check
make install
```
### mpfr
```bash 
wget https://www.mpfr.org/mpfr-current/mpfr-4.2.1.tar.gz
tar -xzvf mpfr-4.2.1.tar.gz
cd mpfr-4.2.1
./configure --prefix=${HOME}
make
make check
make install
```
### gmpfrxx
```bash
wget https://math.berkeley.edu/~wilken/code/gmpfrxx/gmpfrxx.zip
unzip file.zip
cd gmpfrxx
vi README #to see how to install and test gmpfrxx 
```
With the version of `gmp` and `mpfr` I reported, comment lines 74-75 of `example.cpp` before doing `make`, otherwise you get the error:
```bash
example.cpp:74:25: error: ‘GMP_RND_MAX’ was not declared in this scope; did you mean ‘GMP_RNDD’?
   74 |     mpfr_class::set_rnd(GMP_RND_MAX);
      |                         ^~~~~~~~~~~
```
You can also try to install the version of `gmp` and `mpfr` reported in the `README` file and try if everything works without commenting the two lines.  
Finally, remember to specify the paths of `mpfr` and `gmp` in the `Makefile` before doing `make`.
## Usage <a name="usage"></a>
To compile the code, in the `BG_HLT_inverse` directory
```bash
cd build
bash compiler.sh
```
All the available flags have are properly commented in`compiler.sh` file.  

Remember to specify the path of `include` directory in `compiler.sh`
```bash
INVERSE_PATH='${YOUR_PATH}/BG_HLT_inverse/include'
```
and specify the path of `gmpfrxx` directory and all the libraries previously installed
```
g++ -O$optimization -std=c++14 -o main ../src/main.cpp ../lib/params.c ../lib/smear.c ../lib/statistical.c -I${YOUR_PATH}$/gmpfrxx -L${YOUR_PATH}/gmpfrxx -I${INCLUDE_PATH_LIBRARIES} -L${LIB_PATH_LIBRARIES} -lgmpfrxx -lmpfr -lgmpxx -lgmp -lm -lgsl -lgslcblas $Method $Basis $Target $Cov

```
After the compilation, you are ready to use the code:
```bash
./main input_file
```
You can also launch the code without the `input_file` path to see how it was compiled, for instance
```
> ./main
Compiled with COS_SPHAL basis
Compiled with PSEUDO_GAUSS target
Usage: ./main <input_file>
```
## Output <a name="output"></a>
You can choose the path of the output files in the last section of the `input_file`
```
# output files

smear_func_file         Delta_Smear.dat                     
rho_vs_lambda_file      rho_vs_lambda.dat         
output_file             final_res.dat                
rs_vs_lambda_file       RS.dat                              

```
   -  `smear_func_file` contains information about target function reconstruction, this file has 4 columns 
      <div align="center">

      | $\omega$ |$\Delta(\omega,\bar{\omega})$ | $\delta(\omega,\bar{\omega})$ |$\delta(\omega,\bar{\omega})-\Delta(\omega,\bar{\omega})$|
      |-----------|-----------|-----------|-----------|

      </div>

      where $\Delta(\omega,\bar{\omega})$ is the __reconstructed function__ and $\delta(\omega,\bar{\omega})$ the __target function__.

   -  `rho_vs_lambda_file` gives the final result as a function of $\lambda$ and $\delta A_0[g_t(\lambda)]$, where 
      $$
      \delta A_\alpha[g_t(\lambda)]\equiv\frac{\int_{E_0}^\infty d\omega [\Delta(\omega,\bar{\omega})-\delta(\omega,\bar{\omega})]^2 e^{\alpha\omega}}{\int_{E_0}^\infty d\omega [\delta(\omega,\bar{\omega})]^2e^{\alpha\omega}} ,
      $$
      and it has 4 columns  
      <div align="center">

      | $\lambda$ |$\delta A_0[g_t(\lambda)]$ | $\bar{\rho}_\lambda(\bar{\omega})$ |$d\bar{\rho}_\lambda(\bar{\omega})$|
      |-----------|-----------|-----------|-----------| 

      </div>

      where $\bar{\rho}_\lambda(\bar{\omega})$ and $d\bar{\rho}_\lambda(\bar{\omega})$ are the mean value and the error of the smeared reconstructed quantity we are interested in, as a function of $\lambda$.

   - `final_res.dat` is a file with just 1 row and 3 column file of the type
      <div align="center">

      | $\sigma$ |$\bar{\rho}(\bar{\omega})$ | $d\bar{\rho}(\bar{\omega})$ |
      |-----------|-----------|-----------|

      </div>

      where $\sigma$ is the width of $\delta(\omega,\bar{\omega})$ while $\bar{\rho}(\bar{\omega})$ and $d\bar{\rho}(\bar{\omega})$ are the mean value and the error of the final estimation.
   - `rs_vs_lambda_file` is a 4 columns file
      <div align="center">

      | $\lambda$ |$\delta A_\alpha[g_t(\lambda)]$ | $B[g_t(\lambda)]$ |$(B / \delta A_\alpha) [g_t(\lambda)]$|
      |-----------|-----------|-----------|-----------| 

      </div>

      where
      $$
      B[g_t(\lambda)]\equiv \sum_{t,t'=0}^{1/T}\mathrm{Cov}_{tt'} g_t g_{t'}
      $$
      and it is important to fix the parameters `Lside` and `Rside` in the `input_file`.

## Contributors <a name="contributors"></a>
 - M. Naviglio
 - F. D'Angelo

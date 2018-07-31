%1D explicit CFL
C=[-8.8818e-16 - 6.3755e+01i
 0.0000e+00 - 1.2792e+02i
 0.0000e+00 - 2.5588e+02i
 2.8866e-15 - 5.1198e+02i
 3.5527e-15 - 1.0239e+03i
 0.0000e+00 - 2.0478e+03i
 2.1316e-14 - 4.0956e+03i
-4.2633e-14 - 8.1911e+03i];

n=[4:11];
n=n';
C=abs(C);
s=1.73205./C./(2.^(-n/3));

CFL=[6.8457e-02
   4.2987e-02
   2.7076e-02
   1.7049e-02
   1.0741e-02
   6.7665e-03
   4.2626e-03
   2.6853e-03];

dt=[2.7167e-02
   1.3540e-02
   6.7690e-03
   3.3830e-03
   1.6916e-03
   8.4581e-04
   4.2291e-04
   2.1146e-04];

%Lev=9 time=0.5
%----------------------------explicit-------------------------%
2.0000e+00   9.0000e+00   8.0000e-04   4.0471e-04   2.4317e-03  %4.154s

%--------------------------semi-implicit----------------------%
2.0000e+00   9.0000e+00   8.0000e-04   3.1345e-04   2.3766e-03  %4.186s
2.0000e+00   9.0000e+00   9.0000e-04   3.5228e-04   2.3834e-03  %3.661s
2.0000e+00   9.0000e+00   1.0000e-03  3.4279e+268  4.2845e+269  %3.623s

%-------------------------backward euler---------------------%
2.0000e+00   9.0000e+00   8.0000e-04   3.5206e-04   2.3792e-03  %26.240s
2.0000e+00   9.0000e+00   1.0000e-03   4.4042e-04   2.4061e-03  %19.453s
2.0000e+00   9.0000e+00   5.0000e-03   2.2031e-03   3.7889e-03  %11.251s
2.0000e+00   9.0000e+00   1.0000e-02   4.3896e-03   6.3567e-03  %11.814s
2.0000e+00   9.0000e+00   5.0000e-02   2.1537e-02   2.9673e-02  %10.264s

%-------------------------trapezoidal rule-------------------%
2.0000e+00   9.0000e+00   8.0000e-04   1.7828e-04   2.3350e-03  %18.765s
2.0000e+00   9.0000e+00   1.0000e-03   1.7829e-04   2.3352e-03  %16.615s
2.0000e+00   9.0000e+00   5.0000e-03   1.7915e-04   2.3464e-03  %7.355s
2.0000e+00   9.0000e+00   1.0000e-02   1.8182e-04   2.3812e-03  %6.089s
2.0000e+00   9.0000e+00   5.0000e-02   4.8385e-04   2.9890e-03  %5.189s

%-------------------------gauss-legendre---------------------%
2.0000e+00   9.0000e+00   8.0000e-04   1.7825e-04   2.3347e-03  %34.568s
2.0000e+00   9.0000e+00   1.0000e-03   1.7825e-04   2.3347e-03  %29.955s
2.0000e+00   9.0000e+00   5.0000e-03   1.7825e-04   2.3347e-03  %12.809s
2.0000e+00   9.0000e+00   1.0000e-02   1.7826e-04   2.3347e-03  %10.556s
2.0000e+00   9.0000e+00   5.0000e-02   1.7951e-04   2.3511e-03  %8.991s
2.0000e+00   9.0000e+00   8.0000e-02   1.6410e-04   2.1487e-03  %8.840s
2.0000e+00   9.0000e+00   1.0000e-01   1.9478e-04   2.5476e-03  %8.925s
2.0000e+00   9.0000e+00   5.0000e-01   4.0980e-03   5.9492e-03  %8.618s




  




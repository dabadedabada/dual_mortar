Made sure the program works on Windows using MinGW64 compiler:

1.make sure your compiler is MinGW64 Compiler (C++). Go to https://www.mathworks.com/matlabcentral/fileexchange/52848-matlab-support-for-mingw-w64-c-c-compiler and download MinGW-w64 C/C++ Compiler for matlab. Install it.

2. Launch matlab and set the workspace to 'clipper2\private'

3. Run 
mex -setup C++ 

in command window and click MinGW64 Compiler (C++) to use it as your c++ complier

4. Run 
mex('-D__int64=__int64_t','clipper.cpp','mexclipper.cpp')

in MATLAB command window.


Thanks to Yuhan Liu's comment on mathworks.
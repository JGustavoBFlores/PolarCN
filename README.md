    This is a side project for the study of the evolution of the wave function of a particle inside 
    a circular potential, by using the Crank Nicholson method on  Schrodinger's equation in polar
    coordinates.

    This code uses LAPACK to invert matrixes.


    To install LAPACK you will require to install macports
[MacPorts](https://www.macports.org/install.php)
    To install LAPACK go to:
[LAPACK](https://ports.macports.org/port/lapack/)

    To use a LAPACK subroutine, call it from your program, (defining it as an external first),
    and add: '-llapack -lblas' to the end of your compilation command.

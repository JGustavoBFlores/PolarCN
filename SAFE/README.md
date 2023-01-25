    This is a side project for the study of the evolution of the wave function of a particle inside 
    a circular potential, by using the Crank Nicholson method on  Schrodinger's equation in polar
    coordinates.

    Instead of using the regular Schrodinger equation, the atomic units variant will be used.

    This code uses LAPACK to invert matrixes.

    To install LAPACK you will require to install macports:
[MacPorts](https://www.macports.org/install.php)

    To install LAPACK go to:
[LAPACK](https://ports.macports.org/port/lapack/)

    To use a LAPACK subroutine, call it from your program, (defining it as an external first),
    and add: '-llapack -lblas' to the end of your compilation command.

    To run this code you can either run the simulator executable or runit, runit requires a directory 
    of name EVOLUTION and will delete the current version of simulator and compile a newer version 
    from the program.f 

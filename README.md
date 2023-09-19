# 3-body-stab
## To compile use the command:

### Linux or Mac:
gfortran -o g3r_stab g3r_stab.f95 -Ofast

### Windows:
gfortran -o g3r_stab.exe g3r_stab.f95 -Ofast

## To run program use file g3r_input.txt with input parameters

### Linux or Mac:
./g3r_stab < g3r_input.txt

### Windows:
g3r_stab.exe < g3r_input.txt

### For output to file use:
g3r_stab.exe < g3r_input.txt > outputfilename

## Content of the g3r_input.txt file
### FIRST LINE:
     Iwr =writing index, +1 => coordnates are written , -1 no coords
     NS = number of steps / period of inner binary
     IST = output after every IST steps
     max_outer = maximum time (in periods of outer binary)
     a1 =1      these "a" are method parameters. 
     a2 =15     see Mikkola 1997 CemDA 67:145-165

     clight= speed of light. If M_sun=1 , G=1 ,a_Earth=1, then approx c=10^4
     NOTE: if you set: clight=0 => no relativistic terms
    
### SECOND LINE:
     m1, m2, m3 -- masses of component
     
### THIRD LINE:    
     a , e , i, M_0, \Omega, \omega for the inner orbit

### FOURTH LINE:
     Q,  e , i, M_0, \Omega , \omega for the outer orbit

(To make possible to check the relativistic precession
for Mercury the code gave 43"/100 years, i.e. correct)

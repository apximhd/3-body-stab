# 3-body-stab
To compile use the command:

Linux or Mac:
gfortran -o g3r_stab g3r_stab.f95 -Ofast

Windows:
gfortran -o g3r_stab.exe g3r_stab.f95 -Ofast

To run program use file g3r_input.txt with input parameters

Linux or Mac:
./g3r_stab < g3r_input.txt

Windows:
g3r_stab.exe < g3r_input.txt

For output to file use:
g3r_stab.exe < g3r_input.txt > outpitfilename

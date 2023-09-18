set terminal wxt
set grid
set xlabel 'e_{12}'
set ylabel 'cos ι_{mut}'
plot 'TAEI_RK45_J' using 4:9 with lines lw 1 title "cos ι_{mut} — e_{12}"
pause mouse close

set title "Numerical integrator accuracy"

set key autotitle columnhead

set grid

set style line 1 lc rgb '#0060ad' lw 2 pt 5 ps 1.5   # --- blue
set style line 2 lc rgb '#dd181f' lw 2 pt 7 ps 1.5   # --- red
set style line 3 lc rgb '#339900' lw 2 pt 9 ps 1.5   # --- green
set style line 4 lc rgb '#990066' lw 2 pt 11 ps 1.5   # --- magenta
set style line 5 lc rgb '#FF6633' lw 2 pt 13 ps 1.5   # --- orange

plot file index 0 with linespoints ls 1, \
     file index 1 with linespoints ls 2, \
     file index 2 with linespoints ls 3, \
     file index 3 with linespoints ls 4, \
     file index 4 with linespoints ls 5

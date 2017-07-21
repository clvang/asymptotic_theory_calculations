set xlabel "\epsilon"
set ylabel "f(\epsilon)"
m="./asymTheoryData"
set nokey
set grid
set title 'Plot of Asymptotic Theory'
plot m using 1:2 with linespoints
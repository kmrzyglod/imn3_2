set terminal png size 800,600 enhanced font 'Helvetica,12'
set o 'IntegralZad1.png'
set xl 'iteracja'
set yl 'wrtosc calki'
set size ratio -1
set title 'Wartosc calki w kolejnych iteracjach'
set log x
plot 'IntegralZad1.txt' with line 

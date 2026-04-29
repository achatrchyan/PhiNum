run = int(system("cat lastrun.txt"))
base = sprintf("./output/%d", run)

set terminal gif animate delay 10 size 1700, 500
set output sprintf("%s/2dz-0-animation.gif", base)

unset key
set xlabel 'x'
set ylabel 'y'

stats sprintf("%s/2dz-0.txt", base) nooutput

j = 0
while (j <= 392) {
  set multiplot layout 1,3

  set cbrange [-3:8]
  set title 'Dilaton field'
  plot sprintf("%s/2dz-%d.txt", base, j) u 1:2:4 w image t ''

  set cbrange [-6:6]
  set title 'Axion field'
  plot sprintf("%s/2dz-%d.txt", base, j) u 1:2:5 w image t ''

  set cbrange [0:20]
  set title 'Axion energy density'
  plot sprintf("%s/2dz-%d.txt", base, j) u 1:2:6 w image t ''

  unset multiplot
  j = j + 1
}

set output

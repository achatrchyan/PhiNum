run = int(system("cat lastrun.txt"))
base = sprintf("./output/%d", run)

tmpfile = sprintf("%s/2dz-0-animation.tmp.gif", base)
outfile = sprintf("%s/2dz-0-animation.gif", base)
set terminal gif animate delay 10 size 1700, 500
set output tmpfile

unset key
set xlabel 'x'
set ylabel 'y'

# Count available files
nframes = int(system(sprintf("ls %s/2dz-*.txt 2>/dev/null | wc -l", base))) - 1

stats sprintf("%s/2dz-0.txt", base) nooutput

j = 0
dt = 0.2*0.1*20
while (j <= nframes) {
  set multiplot layout 1,3

  t = real(system(sprintf("awk 'NR==1{print $3}' %s/2dz-%d.txt", base, j)))
  set label 1 sprintf("t = %.2f", t) at screen 0.01, 0.95 front

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

  unset label 1

  j = j + 1
}

set output
system(sprintf("mv %s %s", tmpfile, outfile))
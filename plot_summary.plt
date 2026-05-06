run = int(system("cat lastrun.txt"))
base = sprintf("./output/%d", run)

# --- find latest 2dz snapshot ---
nframes = int(system(sprintf("ls %s/2dz-*.txt 2>/dev/null | sed 's/.*2dz-//' | sed 's/\\.txt//' | sort -n | tail -1", base)))
t = real(system(sprintf("awk 'NR==1{print $3}' %s/2dz-%d.txt", base, nframes)))

# --- find latest Distrib-1-x.txt (largest integer x) ---
latest_distrib = int(system(sprintf("ls %s/Distrib-1-*.txt 2>/dev/null | sed 's/.*Distrib-1-//' | sed 's/\\.txt//' | sort -n | tail -1", base)))

# --- output ---
set terminal pngcairo size 1500, 900
set output sprintf("%s/summary.png", base)

set multiplot layout 2,3 title sprintf("t = %.4f", t)
unset key
set xlabel 'x'
set ylabel 'y'

# --- Row 1: animation panels ---
set cbrange [-3:8]
set title 'Dilaton field'
plot sprintf("%s/2dz-%d.txt", base, nframes) u 1:2:4 w image t ''

set cbrange [-6:6]
set title 'Axion field'
plot sprintf("%s/2dz-%d.txt", base, nframes) u 1:2:5 w image t ''

set cbrange [0:20]
set title 'Axion energy density'
plot sprintf("%s/2dz-%d.txt", base, nframes) u 1:2:6 w image t ''

set title sprintf("momentum distribution")
set xlabel 'p'
set ylabel 'f(p)xp^4'
set key bottom left
set log
unset cbrange
plot \
  sprintf("%s/Distrib-1-0.txt", base)                   u 1:6 w l lw 2 t 't=0', \
  sprintf("%s/Distrib-1-%d.txt", base, latest_distrib)  u 1:6 w l lw 2 t sprintf("t=%d", latest_distrib)

# --- Panel 3: SO file, col 1 vs col ?? ---
unset log
unset key
set title 'Pocket volume'
set xlabel 't'
set ylabel 'col ?'         # <-- change ylabel
plot sprintf("%s/SO-1.txt", base) u 1:14 w l lw 2 t ''   # <-- change :2 to your column

# --- Panel 4: SO file, col 1 vs col ?? ---
set title 'Energy'
set xlabel 't'
set ylabel 'col ?'         # <-- change ylabel
plot sprintf("%s/SO-1.txt", base) u 1:11 w l lw 2 t ''   # <-- change :3 to your column


unset multiplot
set output
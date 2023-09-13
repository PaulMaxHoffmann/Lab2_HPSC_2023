set term png
set output "plot.png"
set xlabel "X-axis label" offset 1,-2
set ylabel "Y-axis label" offset 0,-2
set zlabel "Z-axis label"
set title "Diffusion Equation" offset 0,3.5
set view 60, 120

set lmargin at screen 0.30
set rmargin at screen 0.75
set bmargin at screen 0.30
set tmargin at screen 0.75
set size ratio -1  # Preserve aspect ratio


# Set the position of the key (legend)
set key at screen 1,0.95

splot 'F_phi_0.plt' w l, 'F_phi_1.plt' w l, 'F_phi_2.plt' w l, 'F_phi_3.plt' w l

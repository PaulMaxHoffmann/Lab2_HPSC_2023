set term png
set output "plot.png"
set xlabel "X-axis label"
set ylabel "Y-axis label"
set zlabel "Z-axis label"
set title "Diffusion Equation"
splot 'F_phi_0.plt' w l, 'F_phi_1.plt' w l, 'F_phi_2.plt' w l, 'F_phi_3.plt' w l

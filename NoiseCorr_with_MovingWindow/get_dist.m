clear ; clc

sta = loadtxt('./sta/sta_obs.txt')


distance_ab = distance(coordinate_a(1),coordinate_a(2),coordinate_b(i,1),coordinate_b(i,2))/180*pi*6371;
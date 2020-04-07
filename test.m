clear all, close all, clc

[Cnodes,Celems,fib]=benchmark_ellipse_linear();

figure
plot3(Cnodes(:,1),Cnodes(:,2),Cnodes(:,3),'+')
grid on
axis equal
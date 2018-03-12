% Metapost TeX macro
% Title: /media/fernando/FER1G-EXT4/FORMULATION/UnBornITPrograms/SIMULACEO/Hu0Omegam/modelo2.dia
% Creator: Dia v0.97.2
% CreationDate: Thu Aug  6 23:54:50 2015
% For: fernando


beginfig(1);
verbatimtex
%&latex
\documentclass{minimal}
\begin{document}
etex
% Define macro for horizontal centering.
vardef hcentered primary P =
  P shifted -(xpart center P, 0)
enddef;
% Define macro for right justification.
vardef rjust primary P =
  P shifted -(xpart (lrcorner P - llcorner P), 0)
enddef;
  % picture(5.950000,-11.050000)(32.890000,-8.267687)
  x = 1.000000cm; y = -1.000000cm;

  t = 1.000000;

% set_line_color 0.000000, 0.000000, 0.000000
% set_line_color 1.000000, 1.000000, 1.000000
  path p;
  p = (6.000000x,9.000000y)--(6.000000x,11.000000y)--(9.110000x,11.000000y)--(9.110000x,9.000000y)--cycle;
  fill p withcolor (1.000000,1.000000,1.000000);
% set_linewidth 0.100000
% set_line_color 0.000000, 0.000000, 0.000000
  draw (6.000000x,9.000000y)--(6.000000x,11.000000y)--(9.110000x,11.000000y)--(9.110000x,9.000000y)--cycle
    withpen pencircle scaled 0.1000x;
% set_line_color 0.000000, 0.000000, 0.000000
  draw hcentered btex {\usefont{OT1}{phv}{m}{n} $U_0$} etex scaled 1.64t shifted (7.555000x,10.195000y);
  path p;
  p = (21.000000x,9.000000y)--(21.000000x,11.000000y)--(27.050000x,11.000000y)--(27.050000x,9.000000y)--cycle;
  fill p withcolor (1.000000,1.000000,1.000000);
% set_linewidth 0.100000
% set_line_color 0.000000, 0.000000, 0.000000
  draw (21.000000x,9.000000y)--(21.000000x,11.000000y)--(27.050000x,11.000000y)--(27.050000x,9.000000y)--cycle
    withpen pencircle scaled 0.1000x;
% set_line_color 0.000000, 0.000000, 0.000000
  draw hcentered btex {\usefont{OT1}{phv}{m}{n} $f(\Omega_M)$} etex scaled 1.64t shifted (24.025000x,10.195000y);
% set_line_color 0.000000, 0.000000, 0.000000
  draw btex {\usefont{OT1}{phv}{m}{n} $U_0$} etex scaled 1.95841t shifted (10.000000x,9.000000y);
% set_linewidth 0.100000
% set_line_color 0.000000, 0.000000, 0.000000
  draw (27.097107x,10.000000y)--(29.388197x,10.000000y)
    withpen pencircle scaled 0.1000x;
% set_linewidth 0.100000
% set_line_color 0.000000, 0.000000, 0.000000
  draw (29.388197x,10.250000y)--(29.888197x,10.000000y)--(29.388197x,9.750000y)--cycle
    withpen pencircle scaled 0.1000x;
% set_line_color 0.000000, 0.000000, 0.000000
  draw btex {\usefont{OT1}{phv}{m}{n} $\hat{U}_0$} etex scaled 1.95841t shifted (28.000000x,9.000000y);
  path p;
  p = (12.991250x,9.000000y)--(12.991250x,11.000000y)--(17.008750x,11.000000y)--(17.008750x,9.000000y)--cycle;
  fill p withcolor (1.000000,1.000000,1.000000);
% set_linewidth 0.100000
% set_line_color 0.000000, 0.000000, 0.000000
  draw (12.991250x,9.000000y)--(12.991250x,11.000000y)--(17.008750x,11.000000y)--(17.008750x,9.000000y)--cycle
    withpen pencircle scaled 0.1000x;
% set_line_color 0.000000, 0.000000, 0.000000
  draw hcentered btex {\usefont{OT1}{phv}{m}{n} $g(U_0)$} etex scaled 1.64t shifted (15.000000x,10.195000y);
% set_linewidth 0.100000
% set_line_color 0.000000, 0.000000, 0.000000
  draw (9.160419x,10.000000y)--(12.329734x,10.000000y)
    withpen pencircle scaled 0.1000x;
% set_linewidth 0.100000
% set_line_color 0.000000, 0.000000, 0.000000
  draw (12.329734x,10.250000y)--(12.829734x,10.000000y)--(12.329734x,9.750000y)--cycle
    withpen pencircle scaled 0.1000x;
% set_linewidth 0.100000
% set_line_color 0.000000, 0.000000, 0.000000
  draw (17.058498x,10.000000y)--(20.341700x,10.000000y)
    withpen pencircle scaled 0.1000x;
% set_linewidth 0.100000
% set_line_color 0.000000, 0.000000, 0.000000
  draw (20.341700x,10.250000y)--(20.841700x,10.000000y)--(20.341700x,9.750000y)--cycle
    withpen pencircle scaled 0.1000x;
% set_line_color 0.000000, 0.000000, 0.000000
  draw btex {\usefont{OT1}{phv}{m}{n} $\Omega_M$} etex scaled 1.95841t shifted (18.000000x,9.000000y);
endfig;
end;

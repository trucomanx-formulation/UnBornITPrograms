% Metapost TeX macro
% Title: /media/fernando/FER1G-EXT4/FORMULATION/UnBornITPrograms/SIMULACEO/Hu0Omegam/modelo.dia
% Creator: Dia v0.97.2
% CreationDate: Thu Aug  6 22:59:40 2015
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
  % picture(2.395000,-18.050000)(34.912312,-0.267688)
  x = 1.000000cm; y = -1.000000cm;

  t = 1.000000;

% set_line_color 0.000000, 0.000000, 0.000000
% set_line_color 1.000000, 1.000000, 1.000000
  path p;
  p = (2.445000x,9.000000y)--(2.445000x,10.900000y)--(5.555000x,10.900000y)--(5.555000x,9.000000y)--cycle;
  fill p withcolor (1.000000,1.000000,1.000000);
% set_linewidth 0.100000
% set_line_color 0.000000, 0.000000, 0.000000
  draw (2.445000x,9.000000y)--(2.445000x,10.900000y)--(5.555000x,10.900000y)--(5.555000x,9.000000y)--cycle
    withpen pencircle scaled 0.1000x;
% set_line_color 0.000000, 0.000000, 0.000000
  draw hcentered btex {\usefont{OT1}{phv}{m}{n} $U_0$} etex scaled 1.64t shifted (4.000000x,10.145000y);
% set_line_color 0.000000, 0.000000, 0.000000
  draw btex {\usefont{OT1}{phv}{m}{n} $U_1$} etex scaled 1.95841t shifted (19.000000x,2.000000y);
% set_line_color 0.000000, 0.000000, 0.000000
  draw btex {\usefont{OT1}{phv}{m}{n} $U_2$} etex scaled 1.95841t shifted (19.000000x,6.950000y);
% set_line_color 0.000000, 0.000000, 0.000000
  draw btex {\usefont{OT1}{phv}{m}{n} $U_m$} etex scaled 1.95841t shifted (19.000000x,16.000000y);
% set_linewidth 0.100000
% set_line_color 0.000000, 0.000000, 0.000000
  draw (5.604614x,9.966046y)--(9.000000x,10.000000y)--(9.000000x,3.050000y)--(15.388212x,3.004370y)
    withpen pencircle scaled 0.1000x;
% set_linewidth 0.100000
% set_line_color 0.000000, 0.000000, 0.000000
  draw (15.389998x,3.254364y)--(15.888199x,3.000799y)--(15.386427x,2.754376y)--cycle
    withpen pencircle scaled 0.1000x;
% set_linewidth 0.100000
% set_line_color 0.000000, 0.000000, 0.000000
  draw (5.604614x,9.966046y)--(9.000000x,10.000000y)--(9.000000x,8.000000y)--(15.388212x,7.954370y)
    withpen pencircle scaled 0.1000x;
% set_linewidth 0.100000
% set_line_color 0.000000, 0.000000, 0.000000
  draw (15.389998x,8.204364y)--(15.888199x,7.950799y)--(15.386427x,7.704376y)--cycle
    withpen pencircle scaled 0.1000x;
% set_linewidth 0.100000
% set_line_color 0.000000, 0.000000, 0.000000
  draw (5.604614x,9.966046y)--(9.000000x,10.000000y)--(9.000000x,17.000000y)--(15.388212x,16.954370y)
    withpen pencircle scaled 0.1000x;
% set_linewidth 0.100000
% set_line_color 0.000000, 0.000000, 0.000000
  draw (15.389998x,17.204364y)--(15.888199x,16.950799y)--(15.386427x,16.704376y)--cycle
    withpen pencircle scaled 0.1000x;
% set_line_color 0.000000, 0.000000, 0.000000
  draw btex {\usefont{OT1}{phv}{m}{n} ...} etex scaled 3.61597t shifted (17.000000x,11.950000y);
  path p;
  p = (22.000000x,2.000000y)--(22.000000x,18.000000y)--(26.000000x,18.000000y)--(26.000000x,2.000000y)--cycle;
  fill p withcolor (1.000000,1.000000,1.000000);
% set_linewidth 0.100000
% set_line_color 0.000000, 0.000000, 0.000000
  draw (22.000000x,2.000000y)--(22.000000x,18.000000y)--(26.000000x,18.000000y)--(26.000000x,2.000000y)--cycle
    withpen pencircle scaled 0.1000x;
% set_line_color 0.000000, 0.000000, 0.000000
  draw hcentered btex {\usefont{OT1}{phv}{m}{n} Decoder} etex scaled 1.64t shifted (24.000000x,10.195000y);
% set_linewidth 0.100000
% set_line_color 0.000000, 0.000000, 0.000000
  draw (18.000000x,3.000000y)--(21.388197x,3.000000y)
    withpen pencircle scaled 0.1000x;
% set_linewidth 0.100000
% set_line_color 0.000000, 0.000000, 0.000000
  draw (21.388197x,3.250000y)--(21.888197x,3.000000y)--(21.388197x,2.750000y)--cycle
    withpen pencircle scaled 0.1000x;
% set_linewidth 0.100000
% set_line_color 0.000000, 0.000000, 0.000000
  draw (18.000000x,7.950000y)--(21.388244x,7.992353y)
    withpen pencircle scaled 0.1000x;
% set_linewidth 0.100000
% set_line_color 0.000000, 0.000000, 0.000000
  draw (21.385120x,8.242334y)--(21.888205x,7.998603y)--(21.391369x,7.742373y)--cycle
    withpen pencircle scaled 0.1000x;
% set_line_color 0.000000, 0.000000, 0.000000
  draw btex {\usefont{OT1}{phv}{m}{n} $U_0$} etex scaled 1.95841t shifted (6.000000x,9.000000y);
% set_linewidth 0.100000
% set_line_color 0.000000, 0.000000, 0.000000
  draw (26.041016x,10.000000y)--(28.388197x,10.000000y)
    withpen pencircle scaled 0.1000x;
% set_linewidth 0.100000
% set_line_color 0.000000, 0.000000, 0.000000
  draw (28.388197x,10.250000y)--(28.888197x,10.000000y)--(28.388197x,9.750000y)--cycle
    withpen pencircle scaled 0.1000x;
% set_line_color 0.000000, 0.000000, 0.000000
  draw btex {\usefont{OT1}{phv}{m}{n} $\hat{U}_0$} etex scaled 1.95841t shifted (30.000000x,10.000000y);
% set_linewidth 0.100000
% set_linewidth 0.100000
  path p;
  p = (18.000000x,3.000000y)..(17.000000x,4.000000y)..(16.000000x,3.000000y)..(17.000000x,2.000000y)..cycle;
  fill p withcolor (1.000000,1.000000,1.000000);
% set_line_color 0.000000, 0.000000, 0.000000
  draw (18.000000x,3.000000y)..(17.000000x,4.000000y)..(16.000000x,3.000000y)..(17.000000x,2.000000y)..cycle
    withpen pencircle scaled 0.1000x;
% set_line_color 0.000000, 0.000000, 0.000000
  draw (16.292900x,2.292900y)--(17.707100x,3.707100y)
    withpen pencircle scaled 0.1000x;
% set_line_color 0.000000, 0.000000, 0.000000
  draw (16.292900x,3.707100y)--(17.707100x,2.292900y)
    withpen pencircle scaled 0.1000x;
% set_linewidth 0.100000
% set_line_color 0.000000, 0.000000, 0.000000
  draw (17.000000x,1.000000y)--(17.000000x,1.388197y)
    withpen pencircle scaled 0.1000x;
% set_linewidth 0.100000
% set_line_color 0.000000, 0.000000, 0.000000
  draw (16.750000x,1.388197y)--(17.000000x,1.888197y)--(17.250000x,1.388197y)--cycle
    withpen pencircle scaled 0.1000x;
% set_line_color 0.000000, 0.000000, 0.000000
  draw btex {\usefont{OT1}{phv}{m}{n} $E_1$} etex scaled 1.95841t shifted (17.000000x,1.000000y);
% set_linewidth 0.100000
% set_linewidth 0.100000
  path p;
  p = (18.000000x,7.950000y)..(17.000000x,8.950000y)..(16.000000x,7.950000y)..(17.000000x,6.950000y)..cycle;
  fill p withcolor (1.000000,1.000000,1.000000);
% set_line_color 0.000000, 0.000000, 0.000000
  draw (18.000000x,7.950000y)..(17.000000x,8.950000y)..(16.000000x,7.950000y)..(17.000000x,6.950000y)..cycle
    withpen pencircle scaled 0.1000x;
% set_line_color 0.000000, 0.000000, 0.000000
  draw (16.292900x,7.242900y)--(17.707100x,8.657100y)
    withpen pencircle scaled 0.1000x;
% set_line_color 0.000000, 0.000000, 0.000000
  draw (16.292900x,8.657100y)--(17.707100x,7.242900y)
    withpen pencircle scaled 0.1000x;
% set_linewidth 0.100000
% set_line_color 0.000000, 0.000000, 0.000000
  draw (17.000000x,5.950000y)--(17.000000x,6.338197y)
    withpen pencircle scaled 0.1000x;
% set_linewidth 0.100000
% set_line_color 0.000000, 0.000000, 0.000000
  draw (16.750000x,6.338197y)--(17.000000x,6.838197y)--(17.250000x,6.338197y)--cycle
    withpen pencircle scaled 0.1000x;
% set_line_color 0.000000, 0.000000, 0.000000
  draw btex {\usefont{OT1}{phv}{m}{n} $E_2$} etex scaled 1.95841t shifted (17.000000x,5.950000y);
% set_linewidth 0.100000
% set_line_color 0.000000, 0.000000, 0.000000
  draw (18.000000x,16.950000y)--(21.388244x,16.992353y)
    withpen pencircle scaled 0.1000x;
% set_linewidth 0.100000
% set_line_color 0.000000, 0.000000, 0.000000
  draw (21.385120x,17.242334y)--(21.888205x,16.998603y)--(21.391369x,16.742373y)--cycle
    withpen pencircle scaled 0.1000x;
% set_linewidth 0.100000
% set_linewidth 0.100000
  path p;
  p = (18.000000x,16.950000y)..(17.000000x,17.950000y)..(16.000000x,16.950000y)..(17.000000x,15.950000y)..cycle;
  fill p withcolor (1.000000,1.000000,1.000000);
% set_line_color 0.000000, 0.000000, 0.000000
  draw (18.000000x,16.950000y)..(17.000000x,17.950000y)..(16.000000x,16.950000y)..(17.000000x,15.950000y)..cycle
    withpen pencircle scaled 0.1000x;
% set_line_color 0.000000, 0.000000, 0.000000
  draw (16.292900x,16.242900y)--(17.707100x,17.657100y)
    withpen pencircle scaled 0.1000x;
% set_line_color 0.000000, 0.000000, 0.000000
  draw (16.292900x,17.657100y)--(17.707100x,16.242900y)
    withpen pencircle scaled 0.1000x;
% set_linewidth 0.100000
% set_line_color 0.000000, 0.000000, 0.000000
  draw (17.000000x,14.950000y)--(17.000000x,15.338197y)
    withpen pencircle scaled 0.1000x;
% set_linewidth 0.100000
% set_line_color 0.000000, 0.000000, 0.000000
  draw (16.750000x,15.338197y)--(17.000000x,15.838197y)--(17.250000x,15.338197y)--cycle
    withpen pencircle scaled 0.1000x;
% set_line_color 0.000000, 0.000000, 0.000000
  draw btex {\usefont{OT1}{phv}{m}{n} $E_m$} etex scaled 1.95841t shifted (17.000000x,14.950000y);
% set_line_color 0.000000, 0.000000, 0.000000
  draw btex {\usefont{OT1}{phv}{m}{n} $P(e_1=1)=P_s$} etex scaled 1.95841t shifted (10.000000x,2.000000y);
% set_line_color 0.000000, 0.000000, 0.000000
  draw btex {\usefont{OT1}{phv}{m}{n} $P(e_2=1)=P_s$} etex scaled 1.95841t shifted (10.000000x,7.000000y);
% set_line_color 0.000000, 0.000000, 0.000000
  draw btex {\usefont{OT1}{phv}{m}{n} $P(e_m=1)=P_s$} etex scaled 1.95841t shifted (10.000000x,16.000000y);
% set_line_color 0.000000, 0.000000, 0.000000
  draw btex {\usefont{OT1}{phv}{m}{n} $P(u_0=1)=0.5$} etex scaled 1.95841t shifted (3.000000x,12.000000y);
endfig;
end;

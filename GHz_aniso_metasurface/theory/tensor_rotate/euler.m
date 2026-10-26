clear all; close all; clc;


al = 0;
be = pi/2; 
ga = 0;

Mal = [cos(al), sin(al), 0;
      -sin(al), cos(al), 0;
      0,0,1];

Mbe = [1,0,0;
    0,cos(be), sin(be);
    0, -sin(be),cos(be)];

Mga = [cos(ga), sin(ga), 0;
      -sin(ga), cos(ga), 0;
      0,0,1];

Mga*Mbe*Mal


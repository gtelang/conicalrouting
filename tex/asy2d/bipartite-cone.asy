settings.outformat="pdf";
size(150);
import geometry;

real dx = 0.40;
point O =(0,0);
point A =(1,1.8);
point B =(1.7,1);

real cvx=0.1;
point Ka=(1-cvx)*O + cvx*A;
point Kb=(1-cvx)*O + cvx*B;


segment OA=line(O,A);
segment OB=line(O,B);
segment KA=line(Ka,A);
segment KB=line(Kb,B);

// Draw the cone arms and origin O along with appropriate labels
draw(OA);
draw(OB); 

// Draw points along segment KA separated by `dx`
point [] kapts = {};
point [] kbpts = {};

for (int i=0; i < 5; ++i){
  kapts[i] = curpoint(KA,(i)*dx);
}

// Ditto for KB
for (int i=0; i < 5; ++i){
  kbpts[i] = curpoint(KB,(i)*dx);
}



// Draw  segments joining each point on OA to every point on OB
for (int i=0 ; i<kapts.length ; ++i){
  for (int j=0 ; j<kbpts.length ; ++j)
    {
        draw(kapts[i]--kbpts[j],heavygreen+0.6) ;
    }
}



draw(O--Ka,red+1);
draw(O--Kb,red+1);


real d = 0.1;
draw(baseline("$\Delta x$"),(kapts[1] + kapts[1]*I*d) -- (kapts[2] + kapts[1]*d*I),black,Bars,Arrows(TeXHead), PenMargins,align=NW);
draw(baseline("$\Delta x$"),(kbpts[1] - kbpts[1]*I*d) -- (kbpts[2] - kbpts[1]*d*I),black,Bars,Arrows(TeXHead),PenMargins,align=SE);

for (int i=0 ; i<kapts.length ; ++i){
  dot(kapts[i],blue);
}
for (int i=0 ; i<kbpts.length ; ++i){
  dot(kbpts[i],blue);
}


dot(Ka,blue);
dot(Kb,blue);



label(scale(0.7)*"$\varepsilon$",O,1*SE+E);
label(scale(0.7)*"$\varepsilon$",O,1*NW+N);

label("$A$",A, N);
label("$B$",B, E);

dot(O); label("$O$", O, SW);


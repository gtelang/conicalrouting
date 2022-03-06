size(150);
import geometry;

real dx = 0.40;
point O =(0,0);
point A =(1,1.8);
point B =(1.5,1);
segment OA=line(O,A);
segment OB=line(O,B);

// Draw the cone arms and origin O along with appropriate labels
draw(OA);
draw(OB); dot(O); label("$O$", O, SW);

// Draw points along segment OA separated by `dx`
point [] oapts = {};
point [] obpts = {};

for (int i=0; i < 5; ++i){
  oapts[i] = curpoint(OA,(i+1)*dx);
}

// Ditto for OB
for (int i=0; i < 4; ++i){
  obpts[i] = curpoint(OB,(i+1)*dx);
}



// Draw red segments joining each point on OA to every point on OB
for (int i=0 ; i<oapts.length ; ++i){
  for (int j=0 ; j<obpts.length ; ++j)
    {
        draw(oapts[i]--obpts[j],heavygreen+0.6) ;
    }
}

real d = 0.1;
draw(baseline("$\Delta x$"),(oapts[1] + oapts[1]*I*d) -- (oapts[2] + oapts[1]*d*I),black,Bars,Arrows(TeXHead), PenMargins,align=NW);
draw(baseline("$\Delta x$"),(obpts[1] - obpts[1]*I*d) -- (obpts[2] - obpts[1]*d*I),black,Bars,Arrows(TeXHead),PenMargins,align=SE);

for (int i=0 ; i<oapts.length ; ++i){
  dot(oapts[i]);
}
for (int i=0 ; i<obpts.length ; ++i){
  dot(obpts[i]);
}

label("$A$",A, N);
label("$B$",B, E);

// Code for drawing stalactites and stalagmites
size(300);
import geometry;

// Draw start and target points
point s = (0,0);
point t = (3,1);
dot(s); dot(t);

// Draw surrounding box
draw(box(s,t));

// Draw stalagmites (mites go up)
point   [] mtips  = {(0.5,0.7), (1.2,0.5),  (1.8,0.7)};
segment [] mites = {};

for (int i=0 ; i< mtips.length; ++i){
  mites[i] = line((mtips[i].x,0), mtips[i]);
  draw(mites[i]);
}

// Draw stalactites (tites down)
point   [] ttips = {(0.8,0.2), (1.5,0.3),  (2.5,0.2)};
segment [] tites = {};

for (int i=0 ; i< ttips.length; ++i){
  tites[i] = line((ttips[i].x,1), ttips[i]);
  draw(tites[i]);
}

label("$s$", s, SW );
label("$t$", t, NE );

// // draw starting segment
// draw(s--mtips[0], dashed+purple);

// // draw segments joining tips of tites and mites
// for (int i=0 ; i<mtips.length ; ++i){
//   draw(mtips[i]--ttips[i], dashed+red);
// }

// for (int i=0 ; i<ttips.length-1 ; ++i){
//   draw(ttips[i]--mtips[i+1], dashed+darkgreen);
// }

// //draw ending segment
// draw(ttips[ttips.length-1]--t, dashed+purple);

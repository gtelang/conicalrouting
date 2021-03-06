size(250);
import geometry;

real scale = 1.0;
point B=(0,-1);
point O=(0,0);
segment OB = line(O,B);
int N = 8;
real dtheta = pi/N;
for (int i=0 ; i<N ; ++i){
      real  theta = i*dtheta;
      point h1 = scale*(cos(theta), sin(theta));
      point h1p = scale*(cos(theta - dtheta), sin(theta - dtheta)) ;
      point h2 = -h1; 
      segment hseg = segment(h1,h2);

      if (i==1){
	draw(hseg,(dotted+blue+1.3), arrow=Arrows(TeXHead));
      }
      
      else if (i==2){
	draw(hseg,(blue+1.2), arrow=Arrows(TeXHead));

	markangle(L=scale(0.56)*"$\Delta \theta_i$",  h1p, O, h1, BeginArcArrow(TeXHead));
	markangle(L=scale(0.56)*"$\Delta \theta_i$",  -h1p, O, -h1, BeginArcArrow(TeXHead));
	label(scale(0.7)*"Rocking line ", h1,NE);
      }
      else{
	draw(hseg,dotted+gray+0.99,Arrows(TeXHead)) ;
      }


      // fill a cone and its antipode cone with yellow
      path s1=buildcycle(h1--O,O--h1p,h1p--h1); 
      path s2=buildcycle((-h1)--O,O--(-h1p),(-h1p)--(-h1));

      fill(s1, (((i==4) || (i==5))?green:yellow)+opacity(0.25));
      fill(s2, (((i==4) || (i==5))?white:yellow)+opacity(0.25));
}



point P1=-scale*(cos(3*dtheta),sin(3*dtheta));
point P2=-scale*(cos(5*dtheta),sin(5*dtheta));
fill(P1--O--B--cycle,white);
fill(P2--O--B--cycle,white);

draw(OB,black+1.5);

draw( (P1--O)  ,red+1.5 , MidArrow(HookHead,4.0));
draw( (O--P2)  ,red+1.5, MidArrow(HookHead,4.0) );
//label(scale(0.5)*"Shortest path edges", B, SE+(1,-1));

arrow(scale(0.7)*"Shortest path edge", P1,(-1,-1), arrow=Arrow(DefaultHead,1.5));
arrow(scale(0.7)*"Shortest path edge", P2,(1,-1), arrow=Arrow(DefaultHead,1.5));
label(scale(0.5)*"Segment Obstacle", B,S);

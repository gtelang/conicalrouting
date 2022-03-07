size(250);
settings.outformat="pdf";
import geometry;
import fontsize;

//srand(21);
real T=6;
point start=(0,0);
point target = (T,0);

path p = start;
int L=1;

point b = (1,-L);
point t = (1,L);
segment fobs = line(b,t);

segment fpart = line(b,t);
for (int i=0 ; i<=T-2; ++i){
      segment seg = shift(i)*fpart;
      draw(seg,blue+2+opacity(0.5));
      int j=0;
      while (-L+j <= L){
	dot( b + (i,j), black+4    );
	if (i==0){
	  draw(  start -- (b + (i,j)), black+opacity(0.5), MidArrow(TeXHead));
	}
 	if (i==T-2){
	  draw(  (b + (i,j)) -- target, black+opacity(0.5), MidArrow(TeXHead));
	}

	j=j+1;
	
      }


      
      for(int j1=(-L) ; j1 <=L ; ++j1) {
	for (int j2=(-L) ; j2<=L  ; ++j2){
	        if (i>=1){
		  draw(  (i,j2) -- (i+1,j1), black+opacity(0.5), Arrow(TeXHead, position=0.60));}
	  }
      }


	if (i>=1){
	int yr = rand()%(2*L) - L;
	p = p--(i,yr);
     }


}

int yr = rand()%(2*L) - L;
p = p--(T-1,yr)--target;




draw(p,heavygreen+2, Arrow(DefaultHead, position=0.55,size=6.5));
draw(p,heavygreen+2, Arrow(DefaultHead, position=1.65,size=6.5));
draw(p,heavygreen+2, Arrow(DefaultHead, position=2.65,size=6.5));
draw(p,heavygreen+2, Arrow(DefaultHead, position=3.65,size=6.5));
draw(p,heavygreen+2, Arrow(DefaultHead, position=4.65,size=6.5));
draw(p,heavygreen+2, Arrow(DefaultHead, position=5.65,size=6.5));

dot(start,red+6);
dot(target,red+6);

label( scale(0.7)*rotate(-25)*"$w_1+w_2 \geq \ell$"  , (start+b)/2,  5.5*NE, fontsize(8)   );
path qstart=brace((0.05,-0.05), (2,-1),.1); 
path q1=brace((2.9,-1.2), (1,-1.2),.1); label("$w_2+w_3 \geq \ell$", (2,-1.5), fontsize(6));
path q2=brace((2,1.2), (4,1.2),.1);label("$w_3 +w_4\geq \ell$", (3,1.5), fontsize(6));
path q3=brace((5,-1.2), (3.1,-1.2),.1);label("$w_4+w_5 \geq \ell$", (4,-1.5), fontsize(6));
path qtarget=brace((6,0), (4,-1),.1); label(scale(0.7)*rotate(25)*"$w_5 + w_6 \geq \ell$", (5.4,-0.2), fontsize(8));

draw(qstart,black+1bp+opacity(0.8));
draw(qtarget,black+1bp+opacity(0.8));
draw(q1,black+1bp);
draw(q2,black+1bp);
draw(q3,black+1bp);
label("$\ell  = 1.3$ \qquad $W=2$", (3,-2.0), fontsize(9));
label(scale(0.5)*"(lookahead parameter) \qquad (window size)", (2.9,-2.3), fontsize(9));

label(scale(0.4)*minipage("Shortest Path with \newline weight constraints",4cm),(-0.5,-1));

arrow((start+b)/2,(-0.5,-0.5),black+0.3, ArcArrow);

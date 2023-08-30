#######################################
#
# BioShapes
# Maintainer: L. Mada
#
# https://github.com/discoleo/BioShapes
#
# Continuation of:
# 1. Bachelor Thesis (2022-2023)
# Candidate: Adrian Cotoc
# Faculty of Mathematics and Informatics, UVT
#
# Coordinator:
#   Prof. Daniela Zaharie
#   Dr. med. Leonard Mada (Syonic SRL)
#   in collaboration with Syonic SRL
# GitHub: https://github.com/Adi131313/BioShapes
#
# 2. Bachelor Thesis: Darian Voda (2021-2022)


##########################
#### Helper Functions ####
#### Demo / Examples  ####


### Various BioShapes
#' @export
examples.bioshapes = function(col = list("#48B000", 1, 1,
		c("blue", "red"), c("purple", "orange")),
		lwd=2, y.txt = c(6, 0), axt=c(1,2)) {
	if(length(lwd) == 1) lwd = rep(lwd, 5);
	# Plot
	plot.base(xlim=c(-1,10), ylim=c(-1,10), axt=axt);

	### Row 1:

	### Ex 1: Liposome
	lst = liposomes(c(30, 17), r=0.15, phi=c(0, pi/34), d=0.1, center=c(0.5, 8));
	lines(lst, lwd=lwd, col=col[[1]]);
	text(0.5, y.txt[1], "Liposome");


	### Ex 2: Brush-Border Cell
	p1 = c(3, 6.5)
	cell = cell.BrushBorder(p1, w=2, h=2.5, A=1/2, lwd = lwd[2]);
	lines(cell, col=col[[2]]);
	text(4, y.txt[1], "Brush-Border Cell");


	### Ex 3: Smooth Muscles / Connective Tissue
	# Cell 1:
	fill.nc = "#0080FF";
	lst = cell.SmoothMuscle(c(6, 9), c(8, 9), r = 0.3, fill="#F04832",
		col.nc = NA, fill.nc=fill.nc);
	lines(lst, lwd=lwd[3], col=col[[3]]);
	# Cell 2:
	# - someone in the lab sipped some toluidine blue in the haematoxylin;
	lst = cell.SmoothMuscle(c(6, 9), c(6.5, 8), r = 0.4, fill="#D03264",
		col.nc = NA, fill.nc=fill.nc);
	lines(lst, lwd=lwd[3], col=col[[3]]);
	text(7.5, y.txt[1], "Smooth Muscles");


	### Row 2:

	### Ex 4: Vertical Helix
	# Note: explicit function for DNA also available!
	p1 = c(2.5, 0.5); p2 = c(p1[1], 5);
	lst1 = helix(p1, p2);
	lst2 = helix(p1, p2, phi=-pi/2);
	lines(lst1, col=col[[4]][1], lwd=lwd[4])
	lines(lst2, col=col[[4]][2], lwd=lwd[4])
	text(2.5, y.txt[2], "Helix/DNA");

	### Ex 5: Vertical Spirals
	p1 = c(5.5, 1.5); p2 = c(p1[1], 4); dx = c(2.25, 0);
	lst1 = spirals(p1, p2)
	lst2 = spirals(p2 + dx, p1 + dx)
	lines(lst1, col=col[[5]][1], lwd=lwd[5])
	lines(lst2, col=col[[5]][2], lwd=lwd[5])
	text(6.5, y.txt[2], "Spirals/Coils");
}


#### Diagram of a Liposome ####
# d = Dimensions of ArrowHead;
#' @export
diagramLiposome = function(lbl = c("Outer lipid layer",
		"Inner lipid layer", "Lipid bilayer"),
		title = "Liposome", lwd=2, d=-0.4, n = c(30, 17), col="#48B000",
		cex.title = 1.5, xy.title = c(0, -6.5), new = TRUE) {

  if(new){
    plot.base(xlim=c(-10,10), ylim=c(-10,10))
  }
  if(length(lbl) == 1) lbl = rep(lbl, 3);
  # Liposome
  lst = liposomes(n, r=0.5, phi=c(0, pi/(2*n[[2]])), d=0.2)
  lines(lst)

  # Title
  if( ! is.null(title)) text(xy.title[1], xy.title[2], title, cex=cex.title);

  # Left arrow
  # TODO: fix (-d)!
  a1 = arrowSimple(x=c(-2.7,-5), y=c(-4.6,-8), d=-d, lwd=lwd);
  text(-5.5, -9, lbl[[1]])

  # Right arrow
  a2 = arrowSimple(x=c(1.4, 5), y=c(-2.4,-7), d=d, lwd=lwd);
  text(5, -8, lbl[[2]])

  # Upper arrow
  a3 = arrowSimple(x=c(0.08, 1), y=c(3.5,8), d=d, lwd=lwd);
  text(1, 9, lbl[[3]])
}


#### Liposome Measurements ####
#' @export
measureLiposome = function(lbl = c("D = 60 nm", "d=50"), center=c(0,0), add=FALSE,
		lwd.arrow=2, lwd=1.5, title="Liposome", xy.title = c(0,-6.5), cex.title=1.5, ...) {
  if( ! add) plot.base(xlim=c(-10,10), ylim=c(-10,10));

  # Liposome
  lst = liposomes(c(30, 17), r=0.5, phi=c(0, pi/34), d=0.2, center=center);
  lines(lst, lwd=lwd, ...);
  text(xy.title[1], xy.title[2], title, cex=cex.title);

  # Measurement arrow D1
  x = center[1] + c(6,6); y = center[2] + c(-5, 5);
  measure(x, y, lwd=lwd.arrow);
  xy = c(7, 0) + center;
  text(xy[1], xy[2], srt=-90, lbl[1]);

  # Measurement arrow d2
  x = center[1] + c(-2,2); y = center[2] + c(0, 0);
  measure(x, y, lwd=lwd.arrow, d=c(-0.5, 0.5));
  xy = c(0,1) + center;
  text(xy[1], xy[2], lbl[2]);
}

### Example Enzyme ####
#' @export
enzymeReaction = function(x = c(2,5), y = c(1,1),
		lbl = c("A", "B", "Enzyme", "Inhibitor"),
		col = c("black", "black", "black", "red", "red"),
		dx=1, dy=c(0.1, 0.1, 0.5), dI= - c(2, 0.75, 2.4), dH=0.5,
		lwd=c(1, 2), scale=1) {
  if(length(y) == 1) y = c(y, y);
  slope = slope(x, y);
  l1 = shiftLine(x, y, d = dy[[1]], slope=slope, scale=scale);
  l2 = shiftLine(rev(x), rev(y), d = - dy[[2]], slope=slope, scale=scale);
  #
  arrowSimple(l1$x, l1$y, d = -dH, d.head=c(0, 0.5), col=col[[1]], lwd=lwd[[1]]);
  arrowSimple(l2$x, l2$y, d = dH, d.head=c(0, - 0.5), col=col[[1]], lwd=lwd[[1]]);
  text(x[1] - dx, y[1], lbl[[1]], col=col[[2]]);
  text(x[2] + dx, y[2], lbl[[2]], col=col[[2]]);
  # Enzyme
  midx = sum(x)/2;
  midy = sum(y)/2;
  mid  = shiftLine(midx, midy, d = dy[[3]], slope=slope);
  text(mid$x, mid$y, lbl[[3]], col=col[[3]]);
  # Inhibitor
  if(length(lbl) > 3) {
    slopeT = -1/slope;
    pI = shiftPoint(c(midx, midy), d=dI, slope=slopeT, scale=scale);
    arrowT(pI[1:2, 1], pI[1:2, 2], col=col[[4]], lwd=lwd[[2]]);
    text(pI[3,1], pI[3,2], lbl[[4]], col=col[[5]]);
  }
}

### Chemistry ###

#' @export
examples.SpiroGons = function(which = 0, R = 4, ngon = c(5,7,0,5)) {
  if(is.na(match(which, 0:4))) stop("Wrong id!");
  if(which == 0) {
    cat(c("Note:\n",
          "Use examples.SpiroGons(which) to run the individual plots.\n"));
    par.old = par(mfrow = c(2,2));
  }
  lim = R + 1; lim = c(-lim, lim); ngon[ngon == 0] = 6;
  plot0 = function() plot.base(xlim=lim, ylim=lim, axt=NULL, asp=1);
  # Auto: In;
  if(which == 0 || which == 1) {
    ng = ngon[1]; n = 16;
    plot0();
    gg = circle.spiro(n=n, ngon=ng, R=R)
    lines(gg)
    gg = circle.spiro(n = n + 8, ngon=ng, R = 2.5, type="out")
    lines(gg)
  }
  # Out
  if(which == 0 || which == 2) {
    ng = ngon[2]; n = 14; # "incidental" coupling of edge;
    plot0();
    gg = circle.spiro(n=n, ngon=ng, R=R)
    lines(gg)
    gg = circle.spiro(n = n + 8, ngon=ng, R = 2.5, type="out")
    lines(gg)
  }
  # Clockwise: Even (ngon = 6)
  if(which == 0 || which == 3) {
    ng = ngon[3];
    plot0();
    n = 15
    gg = circle.spiro(n=n, R=R, ngon=ng, r.adj = 0.015)
    lines(gg)
    n = 10
    gg = circle.spiro(n=n, R = 2, ngon=ng, r.adj = 0.02)
    lines(gg)
  }
  # Clockwise: Odd
  if(which == 0 || which == 4) {
    plot0();
    n = 31; ng = ngon[4];
    gg = circle.spiro(n=n, R=R, ngon=ng, type = "real-clock")
    lines(gg)
    # through Contact points:
    plot.circle(R, col="green", lty=2)
    gg = circle.spiro(n = n - 8, R = R - 1.2, ngon=ng, type = "clock", r.adj = -0.02)
    gg = as.bioshape(lapply(c(seq(10), 21:23), function(id) gg[[id]]));
    lines(gg)
    gg = circle.spiro(n = n - 8, R = R - 1.2, ngon=ng, type = "in")
    gg = as.bioshape(lapply(12:19, function(id) gg[[id]]));
    lines(gg)
    # Centers of polygons:
    plot.circle(R - 1.2, col="red", lty=3)
    gg = circle.spiro(n = n - 10, R = R - 2.5, ngon=ng, type = "in")
    lines(gg)
    plot.circle(R - 2.5, col="red", lty=3)
  }
  if(which == 0) par.old = par(par.old);
}

#' @export
examples.GonsOnCircle = function(R = 3, ngons = c(3,4,5,6)) {
  par.old = par(mfrow = c(2,2))
  lim = R + 1; lim = c(-lim, lim);
  # Plot 1:
  n = ngons[1];
  plot.base(xlim=lim, ylim=lim, axt=NULL)
  ng = circle.ngon(R=R, n=n)
  lines(ng, fill=3)
  ng = circle.ngon(R = R - 1, n=n, clockwise = FALSE)
  lines(ng, fill=4)
  # Plot 2:
  n = ngons[2];
  plot.base(xlim=lim, ylim=lim, axt=NULL)
  ng = circle.ngon(R=R, n=n)
  lines(ng)
  ng = circle.ngon(R = R - 1, n=n, phi = pi/20)
  lines(ng)
  # Plot 3:
  n = ngons[3];
  plot.base(xlim=lim, ylim=lim, axt=NULL)
  ng = circle.ngon(R=R, n=n)
  lines(ng, fill=5)
  # Plot 4:
  n = ngons[4];
  plot.base(xlim=lim, ylim=lim, axt=NULL)
  ng = circle.ngon(R=R, n=n)
  lines(ng)

  par(par.old);
  invisible();
}


### Basic Shapes

### Star shape polygon
#' @export
example.star = function(n = 5, fill = "red"){
  plot.base()
  star = star(n, R = c(3, 1), center = c(5, 5), fill = fill)
  lines(star)
}


#' @export
example.curves = function(R = 5, nr = 20, axt = c(1, 2)) {
  par.old = par(mfrow = c(2,2));

  ### Circular Helix
  lim = c(-R, R) + c(-1,1);
  plot.base(xlim = lim, ylim = lim)
  xy = helix.rad(R=R, n=nr, N=512)
  lines(xy)
  xy = helix.rad(R=R, n=nr, phi=pi/2, N=512)
  lines(xy, col="red")


  ### Flower 1:
  lim = c(-R, R) + c(-2,2);
  plot.base(xlim = lim, ylim = lim)
  xy = helix.rad(R=R, n=nr-9, phi=pi)
  lines(xy)
  xy = helix.rad(R = R - 1, n=nr-9, phi=pi)
  lines(xy, col="red")


  ### Flower 2:
  lim = c(-R, R)*2 + c(-1,1);
  plot.base(xlim = lim, ylim = lim)
  xy = helix.rad(R=R, r=R, n=round(nr/2), N=256)
  lines(xy, col="red")


  ### Other
  lim = c(-R, R) + c(-2,2);
  plot.base(xlim = lim, ylim = lim)
  abline(h=0, col="green", lty=2)
  for(r in c(3,4,6,9)) {
    xy = helix.rad(R = R - 2/r, n=nr-10, phi=pi, r = 2/r);
    lines(xy, col="red");
  }
  par(par.old);
  invisible();
}

### DNA: Double Helix
#' @export
example.dna = function(n = 3, n.lin = 6, phi = c(0, pi),
		lwd = 2, lwd.lines = lwd, col = c("red", "green"), axt = NULL){
	h2 = dna.new(c(1,8), c(1,4), phi = phi,
		n=n, n.lin=n.lin, lwd=lwd, lwd.lines=lwd.lines, col=col);
	plot.base(axt=axt);
	lines(h2);
	invisible(h2);
}

#' @export
example.dna.tests = function(n = 3, phi = c(pi, -pi/2), A = 1,
		x0 = c(0, 8), y0 = c(0, 0), dy = 3,
		dphi = c(0, pi/6, pi/3, pi/2), col = c("red", "green"),
		xlim = NULL, ylim = NULL, axt = c(1,2)) {
	LEN = length(dphi);
	if(length(dy) == 1) {
		dy = seq(0, LEN - 1) * dy;
	}
	if(is.null(xlim)) {
		sgn  = if(x0[1] < x0[2]) c(-1, 1) else c(1, -1);
		xlim = x0 + 2*sgn;
	}
	if(is.null(ylim)) ylim = c(-2, 10);
	plot.base(xlim=xlim, ylim=ylim, axt=axt);
	for(i in seq(LEN)) {
		lines(dna.new(x0, y0 + dy[i], A=A, col=col, n=n, phi = phi + dphi[i]));
	}
}

### Ducts
#' @export
example.ducts = function(n = c(15, 10)) {
  par.old = par(mfrow = c(1,2));

  ### Duct 1: Composite Plot
  lim = c(-20, 20);
  plot.base(xlim=lim, ylim=lim);
  lines(duct(n[1], c(15,18), phi=pi/2/n[1], fill = "#0048C0"))
  lines(duct(n[1], c(10,13), phi=c(0, pi/n[1]), nc.r = NULL,
             fill = "#8080F0"), lwd = 5);
  lines(duct(n[1], c(5,8)))
  abline(h=0, col="green")

  ### Duct 2:
  R = c(3, 5)
  tmp = duct(n = n[2], R = R, center = c(1,3))
  plot.base(xlim = c (-10, 10), ylim = c(-10, 10))
  lines(tmp)

  tmp = ngon.circle(4, N = n[2], R = R[2], center = c(1,3))
  lines(tmp, col = "blue")

  par(par.old);
  invisible();
}

#' @export
example.complexDuct = function(n = 8, lim = c(-10, 10)){
  # warning: just if n is even
  center = c(0, 0)
  radius = c(7, 5, 2);
  h.scale = 1
  R.scale = c(1, 1)
  dr = 0.25

  plot.base(xlim=lim, ylim=lim);
  tmp = duct.complex(n=n);
  lines(tmp);
}


###################

### Examples: Arcs

### Design of Lens: Simple Arcs
#' @export
example.arcs = function(th = pi/3) {
  par.old = par(mfrow = c(2,2));

  ### Half of Middle Arc:
  r = 3; dr = 0.3;
  plot.base()
  plot.circle.arc(r + dr, c(3,3), c(pi - th, pi), lty=2)
  plot.circle.arc(r     , c(3,3), c(pi - th, pi + th), lty=2)
  plot.circle.arc(r - dr, c(3,3), c(pi, pi + th), lty=2)

  ### Curvatures: only with lens()
  cc = c(4,3)
  plot.base()
  plot.circle.arc(3, cc, c(pi - th, pi + th), lty=2)
  plot.circle.arc(4, cc, c(pi - th, pi + th), lty=2)
  plot.circle.arc(5, cc, c(pi - th, pi + th), lty=2)
  
  ### Mirrored:
  plot.base()
  phi = c(2*pi - th, th);
  plot.circle.arc(3, c(3,3), phi, lwd=1, col="#6432B0")
  plot.circle.arc(3, c(3 - 0.2,3), phi + pi, lwd=1, col="#6432B0")
  # exact = r * sin(pi/3) + 3;
  abline(h = 6, lty=2, col="green")
  
  ### Ex 4:
  r = 5;
  plot.base()
  for(ri in seq(1, r, by=1)) {
	plot.circle.arc(ri, c(3,3), c(0, 2*th), lwd=1, col="#6432B0")
  }

  par(par.old);
  invisible();
}


### Arcs by Distance
# dL = distance between the 2 halves;
# col.line = colour of middle line;
#' @export
example.ArcsByDist = function(d = c(0.5, 1, 1.5, 2), dL = 0.5,
		col = c("green", "red"), col.line = "blue", lwd=1, debug=FALSE) {
	plot.arcs = function(x, y, dd) {
		if( ! is.null(col.line)) lines(x, y, col = col.line);
		for(d in dd) {
			sg = if(d < 0) -1 else 1;
			tmp = circle.ArcByDist(x, y, d = d, col=col[1], lwd=lwd);
			lines(tmp);
			tmp = circle.ArcByDist(x, y, d = d + 3*sg, col=col[2], lwd=lwd);
			lines(tmp);
		}
	}
	if(is.logical(debug)) {
		if(debug == FALSE) debug = 0;
	}
	par.old = par(mfrow=c(2,2));
	# V
	plot.base()
	x = c(4, 4); y = c(0, 6);
	if(debug < 0) y = 6 - y;
	plot.arcs(x, y, d);
	plot.arcs(x + dL, y, - d);
	# H
	plot.base()
	x = c(2, 5); y = c(4, 4);
	if(debug < 0) x = 5 - x;
	plot.arcs(x, y, d);
	plot.arcs(x, y - dL, - d);
	# Oblique:
	plot.base()
	x = c(2, 5); y = c(0, 6);
	if(debug < 0) { x = 5 - x; y = 6 - y; }
	plot.arcs(x, y, d);
	xy = shift.ortho(x, y, d = - dL);
	plot.arcs(xy$x, xy$y, - d);
	# Oblique:
	plot.base()
	x = c(5, 2); y = c(0, 6);
	if(debug < 0) { x = 5 - x; y = 6 - y; }
	plot.arcs(x, y, d);
	xy = shift.ortho(x, y, d = - dL);
	plot.arcs(xy$x, xy$y, - d);
	
	par(par.old);
	invisible();
}

### Examples of convex lenses
#' @export
example.lens = function(pos = c(0, 1/2, 1), fill = "#6480D0"){

  par.old = par(mfrow = c(2,2));
  R = 5;
  lens = lens(R = R, x = c(1, 2), y = c(0, 4))
  plot.base()
  lines(lens)

  ### Example 1:
  R = 5;
  plot.base()
  lens = lens(R = R, x = c(1, 2), y = c(0, 4))
  lines(lens)
  #
  lens = lens(R = R, x = c(5, 3), y = c(0, 5))
  lines(lens, col="Red")
  # negative R: semi-concave Lens
  lens = lens(R = c(4,-7), x = c(5, 0), y = c(1, 0) + 6)
  lines(lens, col="#329624")

  ### Example 2: Group of Lenses
  # pos = c(0, poz, 1)
  h = c(2, 3, 4)
  scale.R = c(1, 1.5, 2)
  x = c(0, 6); y = c(3, 2);
  # fill: does NOT work with concave lenses;

  lst = lens.group(x=x, y=y, h=h, pos=pos, l.scale = scale.R, fill=fill)
  plot.base()
  lines(lst)

  ### Example: Lens Group
  # pos = c(0, poz, 1)
  h = c(2, 1.2, 1.5)
  scale.R = c(1, 1.5, 2)
  x = c(0, 6); y = c(0, 4);
  # fill: does NOT work with concave lenses;
  #
  lst = lens.group(x=x, y=y, h=h, pos=pos, l.scale = scale.R, fill=fill)
  plot.base()
  lines(lst)
  lines(x, y, lty=2, lwd=2, col="green")

  par(par.old);
  invisible();
}


### Braces
#' @export
example.braces = function(){
  plot.base(ylim = c(-10,10))
  lines(braces.curly(c(0,0)), lwd=3)
  lines(braces.curly(c(4,0), left=FALSE), lwd=3)
}

#####################

### Virus

#' @export
example.virus = function(N = 12, R = 3, lwd = 10){

  par.old = par(mfrow = c(1,3));

  # Virus with polygons spikes
  virus = virus(N, R = R, center = c(5, 5), lwd = lwd)
  plot.base()
  lines(virus)

  # Virus with circle spikes
  virus = virus(N, R = R, center = c(5, 5), ngon.spike = 0, lwd = lwd)
  plot.base()
  lines(virus)

  # Virus with both polygons and circles
  plot.base(xlim = c(-6, 6), ylim = c(-6, 6))
  tmp = virus2()
  lines(tmp)

  par(par.old);
  invisible();
}



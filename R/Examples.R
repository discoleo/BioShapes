#######################################
#
# Title: BioShapes
# Maintainer: L. Mada
#
# https://github.com/discoleo/BioShapes
#
# Continuation of:
# 1. Bachelor Thesis: Adrian Cotoc (2022-2023)
# Faculty of Mathematics and Informatics, UVT
# Coordinator:
#   Prof. Daniela Zaharie
#   Dr. med. Leonard Mada (Syonic SRL)
#   in collaboration with Syonic SRL
#   GitHub: https://github.com/Adi131313/BioShapes
#
# 2. Bachelor Thesis: Darian Voda (2021-2022)


##########################
#### Helper Functions ####
#### Demo / Examples  ####

#' #export
test = function(...) {
	UseMethod("test");
}
#' #export
test.lines = function(...) {
	UseMethod("test.lines");
}

######################

### Line Intersections

# x, y   = 1st line;
# xB, yB = 2nd line;
# OR x = matrix encoding both lines;
#' @export
test.lines.simple = function(x = c(1, 5), y = c(1, 9), xB = c(2,7), yB = c(7, 3),
		lty = c(1, 2), col = c("black", "blue", "green", "red"),
		ylim = c(0, 12), add = FALSE) {
	if(inherits(x, "matrix")) {
		y = x[,2];
		xB = x[c(3,4)]; xA = x[c(1,2)];
		yB = y[c(3,4)]; yA = y[c(1,2)];
	} else { xA = x; yA = y; }
	#
	if( ! add) plot.base(ylim=ylim);
	lines(xA, yA, col = col[1], lty = lty[[1]]);
	lines(xB, yB, col = col[2], lty = lty[[2]]);
	p = intersect.lines(xA, yA, xB, yB);
	# NONE of the segments:
	if(any(is.na(p$x))) {
		text(sum(xA, xB)/4, sum(yA, yB)/4, "NO", col = "red");
	}
	# Point intersects outside one of the segments
	isInt = is.intersect.lines(p);
	idCol = if(isInt) 3 else 4;
	points(p$x, p$y, col = col[idCol]);
	#
	invisible(p);
}
# Helper:
#' #export
test.lines.list = function(xy, lty = c(1, 2), col = c("black", "blue", "green")) {
	plotf = function(xA, xB, yA, yB) {
		plot.base(ylim = c(0, 12))
		lines(xA, yA, col = col[1], lty = lty[[1]]);
		lines(xB, yB, col = col[2], lty = lty[[2]]);
		p = intersect.lines(xA, yA, xB, yB);
		if(any(is.na(p$x))) {
			text(sum(xA, xB)/4, sum(yA, yB)/4, "NO", col = "red");
		}
		points(p$x, p$y, col = col[3]);
	}
	lapply(xy, function(xy) {
		plotf(xy$xA, xy$xB, xy$yA, xy$yB);
	})
	invisible();
}
#' @export
test.lines.special = function(lty = c(1, 2), col = c("black", "blue", "green")) {
	par.old = par(mfrow = c(2,2))

	### Test
	lst = list(
		list(xA = c(1,6),yA = c(1,8),
			xB = c(1,7), yB = c(4,2)),
		### Special Cases:
		# Overlap
		list(xA = c(1,3), yA = c(1,5),
			xB = c(2,6), yB = c(3,11)),
		# NO
		list(xA = c(1,3), yA = c(1,5),
			xB = c(4,6), yB = c(7,11)),
		# Overlap
		list(xA = c(1,6), yA = c(1,11),
			xB = c(3,5), yB = c(5,9))
	);
	
	test.lines.list(lst);
	
	par(par.old);
	invisible();
}

#####################

### Various BioShapes
#' @export
example.bioshapes = function(n.lipo = c(30, 17), d.lsc = 0, d.lipo = 0.1,
		col = list("#48B000", c("#FFFF80", "#9648C0"), 1,
			c("blue", "red"), c("#649664", "orange")),
		lwd = 2, y.txt = c(6, 0), axt = c(1,2)) {
	if(length(lwd) == 1) lwd = rep(lwd, 5);
	# Plot
	plot.base(xlim=c(-1,10), ylim=c(-1,10), axt=axt);

	### Row 1:

	### Ex 1: Liposome
	lst = liposomes(n.lipo, r = 0.15, phi = c(0, pi/34), center = c(0.5, 8),
		d = d.lipo, d.lsc = d.lsc);
	lines(lst, lwd=lwd, col=col[[1]]);
	text(0.5, y.txt[1], "Liposome");


	### Ex 2: Brush-Border Cell
	p1 = c(3, 6.5)
	cell = cell.BrushBorder(p1, w=2, h=2.5, lwd = lwd[2],
		col = 1, fill = col[[2]][[1]], fill.nc = col[[2]][[2]]);
	lines(cell);
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

	### Ex 4a: dsDNA Helix
	x = c(0.5, 0.5); y = c(0.5, 5);
	lst1 = dna.new(x, y, lwd=lwd[4]);
	lines(lst1)

	### Ex 4b: Vertical Helix
	# Note: explicit function for DNA also available!
	# (see Ex. 4a)
	p1 = c(3, 0.5); p2 = c(p1[1], 5);
	lst1 = helix(p1, p2, col=col[[4]][1], lwd=lwd[4]);
	lst2 = helix(p1, p2, col=col[[4]][2], lwd=lwd[4], phi=-pi/2);
	lines(lst1)
	lines(lst2)
	text(1.5, y.txt[2], "Helix/DNA");

	### Ex 5: Vertical Spirals
	p1 = c(6, 1.5); p2 = c(p1[1], 4); dx = c(2.25, 0);
	lst1 = spirals(p1, p2)
	lst2 = spirals(p2 + dx, p1 + dx)
	lines(lst1, col=col[[5]][1], lwd=lwd[5])
	lines(lst2, col=col[[5]][2], lwd=lwd[5])
	text(7.2, y.txt[2], "Spirals/Coils");
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
  a1 = arrowSimple(x=c(-2.7,-5), y=c(-4.6,-8), d=d, lwd=lwd);
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


### Example Enzymes
#' @export
enzymeReaction = function(x = c(2,5), y = c(1,1),
		lbl = c("A", "B", "Enzyme", "Inhibitor"),
		col = c("blue", "black", "black", "red", "red"),
		dx=1, dy=c(0.1, -0.1, 0.5), dI= - c(2, 0.75, 2.4), dH=0.5,
		lwd=c(1, 2), scale=1, new = FALSE) {
  if(length(y) == 1) y = c(y, y);
  slope = slope(x, y);
  if(new) plot.base();
  arrowDHalf(x, y, d = dy[1:2], col=col[[1]], lwd=lwd[[1]], plot = TRUE);
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

###############

### Ensemble of Shapes

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


##############

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


#############

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


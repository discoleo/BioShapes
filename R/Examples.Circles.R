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


### Examples:
# - Circles;
# - Ellipses;
# - Arcs: e.g. Lenses;


################

################
### Ellipses ###

### Test Ellipses
# - Range & Tangents to given points;
#' @export
test.ellipse.tan = function(x = c(0.5), r = c(1,3), phi = 0, center = c(0,0),
		dx = c(-2, 2), N = 64, col.xlim = "#FF6496", lbl = "Tangents") {
	x0 = range.ellipse.x(r=r, phi=phi, center=center);
	xx = seq(x0[1], x0[2], length.out = N);
	y = solve.ellipse.all(xx, r=r, phi=phi, center=center);
	plot.base(ylim = range(y) + c(-3, 3));
	lines(xx, y[1,])
	lines(xx, y[2,])
	
	sol = sapply(x, function(x) {
		solve.ellipse.all(x=x, r=r, phi=phi, center=center);
	});
	idNA = which(is.na(sol[1,]));
	if(length(idNA) > 0) warning("NAs for solutions: ", paste0(idNA, collapse = ", "));
	# Slope:
	sl = sapply(seq_along(x), function(id) {
		slope.ellipse(x[id], y = sol[, id], r=r, phi=phi, center=center);
	});
	for(id in seq_along(x)) {
		yi = sol[, id]; si = sl[, id];
		# TODO: both L;
		lines.slope(c(x[id], yi[1]), slope = si[1], L = 2*max(dx), col = "green");
		lines.slope(c(x[id], yi[2]), slope = si[2], L = 2*max(dx), col = "green");
		points(rep(x[id], 2), yi, col = "red");
	}
	abline(v = x0, col=col.xlim);
	if( ! is.null(lbl)) {
		text(3, 0, labels = lbl);
	}
}

#' @export
test.ellipse.intersect = function(x0 = c(5.52, 5), y0 = c(3.5, 3),
		phi = pi - pi/3, r = c(3,2), center = c(4,4), verbose = TRUE) {
	plot.base()
	plot.ellipse(r = r, center = center, th = phi, phi = c(0, 2*pi));
	lines.slope(center, tan(phi), L = 8, col = "green");
	lines.slope(center, - 1 / tan(phi), L = 7, col = "blue");
	
	len = length(x0)
	for(i in seq_along(x0)) {
		xy0 = c(x0[i], y0[i]);
		xy  = intersect.ellipse.ParallelLine(r = r, xy = xy0, center=center, phi=phi)
		mid = apply(xy, 2, mean);
		if(verbose) {
			xym = rbind(x = c(xy[,1], mid[1]), y = c(xy[,2], mid[2]));
			colnames(xym) = NULL;
			print(xym);
		}
		lines.slope(mid, tan(phi), L = 10, col = "pink");
		points(c(xy[,1], mid[1]), c(xy[,2], mid[2]), col="red");
	}
	# Center:
	points(center[1], center[2], col="red");
}

#' @export
test.ellipse.bySlope = function(
		d = c(1.5, 1, 1, 1.5, 1),
		rr = c(3, 3, 1/2, 2, 2)) {
	ellf = function(d = 1, rr = 2, col = NULL) {
		lst = solve.ellipse.bySlope(c(2, 7), c(5, 3), d=d, rr=rr)
		lst$phi = c(0, 2*pi);
		if( ! is.null(col)) lst$col = col;
		do.call(plot.ellipse, lst);
	}
	#
	len = length(d);
	col = list("purple", "#D09600", "#3296E0", "#64D032", "black");
	col = rep(col, ceiling(len / length(col)));
	plot.base();
	for(i in seq_along(d)) {
		ellf(d = d[[i]], rr = rr[[i]], col = col[[i]]);
	}
	points(c(2,7), c(5,3), col = "red")
}

######################

######################
### Examples: Arcs ###

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


### Convex Lenses
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

###################

### Derived Shapes:
### Based on Ellipses

#' @export
test.cylinder = function(lty.back = 2, col = "red") {
	plot.base()
	# "/\"
	lines(cylinder(c(1, 3), c(2, 6), lty.back=lty.back, col=col));
	lines(cylinder(c(6, 4), c(2, 6), lty.back=lty.back, col=col));
	# H: not best visual examples
	# (just necessary test)
	lines(cylinder(c(5, 0), c(8, 8), lty.back=lty.back, col=col));
	lines(cylinder(c(2, 7), c(0, 0), lty.back=lty.back, col=col));
}

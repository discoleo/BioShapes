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


###############

###############
### Circles ###

### Hashed Circles
# lty = line type for the hashed-lines: circle 1 vs circles 2:4;
#' @export
test.circle.hash = function(phi = c(pi/5, pi/3, pi/2 + 0.2, pi - 0.2),
		n = list(11, 13, 13, c(13, 11)), r = 2,
		lty = c(2, 1), scale = c(1,2)) {
	par.old = par(mfrow = c(2,2));
	
	### Base:
	plot.circle.hash(n[[1]], phi = phi[1], r=r, scale = scale[1],
		lty = lty[1], col = "red");
	
	### Scale
	scale = scale[2]; lty = lty[2];
	plot.circle.hash(n[[2]], phi = phi[2], r=r, scale=scale, lty=lty);
	#
	plot.circle.hash(n[[3]], phi = phi[3], r=r, scale=scale, lty=lty);
	#
	phi = phi[4]; n = n[[4]];
	len = length(n);
	plot.circle.hash(n[1], phi = phi, r=r, scale=scale, lty=lty);
	if(len > 1) {
		for(id in seq(2, len)) {
			phi_i = phi + pi / id;
			lines(circle.hash(n[id], phi = phi_i, center = c(4,4),
				r=r, scale=scale, lty=lty, col = "red"));
	}}
	
	par(par.old)
}

### Plot Objects formed from circles;
# - convenience function;
# - xy  = chain of circles;
# - pin = hack to set par(pin) = mean(...);
#   Note: asp = 1 is the better approach (and is set automatically);
#' @export
test.FilledCircle = function(xy, r=NULL, R=NULL, lim=NULL, line=TRUE,
                            col="#B0B032", col.line="green", add=FALSE, pin = FALSE, ...) {
	if(missing(xy)) {
		# Basic object:
		cat("Note: Generating object with chained circles!\n");
		xy = circles.OnCircle(11, r = 5);
	}
	if(is.null(r)) {
		r = attr(xy, "r");
		if(is.null(r)) stop("Missing r!");
	} else {
		attr(xy, "r") = r;
	}
  ### New Plot
  if( ! add) {
    x = xy$x; y = xy$y;
    mid = attr(xy, "center");
    if(is.null(lim)) {
      R0  = attr(xy, "R");
      lim = R0 + r + 1;
      lim = c(-lim, lim);
    } else if(length(lim) == 1) {
      lim = c(-lim, lim);
      mid = c(0, 0); # remove center offset;
    } else {
      mid = c(0, 0); # remove center offset;
    }
    plot.base(xlim = lim + mid[1], ylim = lim + mid[2], asp = 1);
  }
  if(pin){
    pin = mean(par("pin")) + 0.25;
    par.old = par(pin = c(pin, pin));
    lines.circles(xy, R=R, line=line, fill=col, col.line=col.line, ...)
    par(par.old);
  }
  else {
    lines.circles(xy, R=R, line=line, fill=col, col.line=col.line, ...)
  }
}


### Quasi-Uniform Points inside Circle
#' @export
test.circle.uniform.text = function(n = 200, phi = 0, d = 0.5) {
	xy = uniform.circle(n[1], phi = phi[1], d = d[1]);
	plot(xy, type="n");
	text(xy, labels=seq(0, n[1]));
	len = length(n);
	if(len > 1) {
		if(length(phi) == 1) phi = rep(phi, len);
		if(length(d) == 1) d = rep(d, len);
		for(id in seq(2, length(n))) {
			ni = n[id];
			xy = uniform.circle(ni, phi = phi[id], d = d[id]);
			text(xy, labels=seq(0, ni), col = "red");
		}
	}
}


#####################

#####################
### Circle Chains ###

# phi = rotation;
#' @export
test.Chains = function(n = 15, phi = pi/n) {
	if(length(phi) == 1) phi = rep(phi, 4);
	par.old = par(mfrow = c(2,2))
	
	### Closed Circles
	# - n & r known, R unknown;
	r = 1
	xy = circles.OnCircle(n, r, phi=phi[[1]]);
	test.FilledCircle(xy);
	text(0, 0, "R = unknown", col = "#32B096");
	text(xy$x, xy$y, seq(n), col = "red");
	
	### Radius of Big Circle: known
	# - small circles: on the big circle;
	# - n & R known, r unknown;
	R = 7
	xy = circles.OnFixedCircle(n, r=R, phi=phi[[2]]);
	test.FilledCircle(xy, R=R);
	text(0, 0, "R = known", col = "#32B096");
	
	### Inside Circle
	# Tangent & Inside given Circle
	# Outer Circle: known
	R = 15
	xy = circles.InFixedCircle(n, R, phi=phi[[3]]);
	test.FilledCircle(xy, R=R)
	
	### Outside Circle
	# Tangent & Outside given Circle
	R = 6
	xy = circles.OutsideFixedCircle(n, R, phi=phi[[4]]);
	test.FilledCircle(xy, R=R);
	
	par(par.old);
}

# phi = rotation;
#' @export
test.ChainsCombined = function(n = 15, phi = pi/n) {
	if(length(phi) == 1) { phi = rep(phi, 4); phi[3] = 0; }
	par.old = par(mfrow = c(2,2));
	
	### Inside Circle / Combined
	# - Circle 1: R unknown
	r = 1
	xy = circles.OnCircle(n, r, phi=phi[[1]]);
	test.FilledCircle(xy)
	# - Circle 2: based on R of Circle 1;
	# Outer Circle: known
	# Inner Circle: unknown
	R = attr(xy, "R") - r;
	xy = circles.InFixedCircle(n, r=R, phi=phi[[1]]);
	test.FilledCircle(xy, add=TRUE, line=FALSE);
	
	### Tangent to outer/inner Chain of Circles
	R = 3.25;
	fill = c("#B0B032", "#90B048");
	lst1 = circles.TanToChainShape(n, R, center = c(4,4), phi=phi[[2]],
		col = NA, fill=fill);
	lst2 = circles.TanToChainShape(n, R, center = c(4,4), phi=phi[[2]],
		col = NA, fill=fill, type = "Outer");
	plot.base()
	lines(lst1)
	lines(lst2)
	
	### Outside Circle: Shifted Center
	# 2 Circles
	R = 6
	mid1 = c(-R, 0); mid2 = mid1 + c(2*R, 0);
	xy1 = circles.InFixedCircle(n, r=R, center=mid1, phi=phi[[3]]);
	xy2 = circles.InFixedCircle(n, r=R, center=mid2, phi=phi[[4]]);
	test.FilledCircle(xy1, R=R, lim = 2*R + 1);
	test.FilledCircle(xy2, R=R, add=TRUE);
	
	# TODO: Example 4;
	
	par(par.old);
}


################
### Ellipses ###

### Test Ellipses
# - Range & Tangents to given points;
#' @export
test.ellipse.tan = function(x = c(2.5, 1.5), r = c(1,3), phi = 0, center = c(2,0),
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
	idNA  = which(is.na(sol[1,]));
	hasNA = length(idNA) > 0;
	if(hasNA) warning("NAs for solutions: ", paste0(idNA, collapse = ", "));
	# Arc:
	toArc = function(x, y) {
		arc.ellipse(rep(x, each=2), as.vector(y), r=r, center=center, theta = phi);
	}
	arc.phi = if(hasNA) {
		toArc(x[- idNA], sol[, - idNA]);
	} else toArc(x, sol);
	arc.phi = matrix(arc.phi, nrow=2);
	lapply(seq(1, length.out = ncol(arc.phi)), function(nc) {
		arc.phi = arc.phi[, nc];
		if(arc.phi[2] - arc.phi[1] >= pi) {
			tmp = arc.phi[1]; arc.phi[1] = arc.phi[2]; arc.phi[2] = tmp;
		}
		plot.ellipse(r=r, center=center, theta = phi, phi = arc.phi,
			lwd = 2, col = "#6432E0");
	})
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
	abline(v = x, lty = 2, col = "#A0D032");
	if( ! is.null(lbl)) {
		text(3, 0, labels = lbl);
	}
}

### Tangents with given slope
# - all tangents are drawn for each of the ellipses;
#' @export
test.ellipse.whereTan = function(slope = c(0, 1/3, 2, Inf), r = c(1,2),
		phi = c(0, pi/3, pi/2, 3*pi/5), center = c(0,0),
		paired = FALSE, verbose = FALSE, dx = c(-2,2)) {
	xlim = max(abs(r)) + 0.5;
	ylim = max(abs(r)) + 0.5;
	xlim = center[1] + c(-xlim, xlim);
	ylim = center[2] + c(-ylim, ylim);
	#
	L = 2*max(dx);
	doTan = function(slope, phi) {
		xy = solve.ellipse.xytan(slope=slope, r=r, phi=phi, center=center);
		if(verbose) print(xy);
		lapply(seq(nrow(xy)), function(id) {
				lines.slope(c(xy[id, 1], xy[id, 2]), slope=slope, L=L, col = "green");
		});
		points(xy[,1], xy[,2], col="red");
	}
	len = length(phi);
	par.old = if(len == 1) par(mfrow = c(1,1)) else par(mfrow = c(2,2));
	if(paired) {
		for(i in seq_along(phi)) {
			plot.base(xlim=xlim, ylim=ylim, asp=1);
			plot.ellipse(r=r, center=center, th = phi[i], phi = c(0, 2*pi));
			doTan(slope = slope[i], phi = phi[i]);
		}
	} else for(phi_i in phi) {
		plot.base(xlim=xlim, ylim=ylim, asp=1);
		plot.ellipse(r=r, center=center, th = phi_i, phi = c(0, 2*pi));
		for(sl in slope) {
			doTan(slope = sl, phi = phi_i);
		}
	}
	
	par(par.old);
	invisible();
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


#################

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

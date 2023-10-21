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


### Helper Functions: Circle-Chains

# - Tangent circles forming a large circle;

circles = function(n, r, center = c(0,0), phi = 0, ...) {
	UseMethod("circles");
}
circles.math = function(n, r, ...) {
	UseMethod("circles.math");
}


### Large Circle of Unknown Radius
# - n small circles each of radius r;
#' @export
circles.OnCircle = function(n, r, center = c(0,0), phi=0) {
  R  = r / sin(pi/n);
  xy = points.circle(n, r=R, center=center, phi=phi);
  attr(xy, "R") = R;
  attr(xy, "r") = r;
  return(xy);
}

### Outside a large circle of given radius
#' @export
circles.OutsideFixedCircle = function(n, r, center = c(0,0), phi=0) {
  r1 = r / (1/sin(pi/n) - 1);
  R1 = r + r1;
  xy = points.circle(n, r=R1, center=center, phi=phi);
  attr(xy, "R") = R1;
  attr(xy, "r") = r1;
  return(xy);
}

### Forming large circle of given radius
#' @export
circles.OnFixedCircle = function(n, r, center = c(0,0), phi=0) {
  r1 = r * sin(pi/n);
  xy = points.circle(n, r=r, center=center, phi=phi);
  # Outer R: reuse same attribute name ???
  attr(xy, "R") = r1 + r;
  attr(xy, "r") = r1;
  return(xy);
}

### Inside a large circle of given radius
#' @export
circles.InFixedCircle = function(n, r, center = c(0,0), phi=0) {
  r1 = r / (1 + 1/sin(pi/n));
  R  = r - r1;
  xy = points.circle(n, r=R, center=center, phi=phi);
  attr(xy, "R") = R;
  attr(xy, "r") = r1;
  return(xy);
}


### Circles Tangent to a Chain of Circles
# n = number of circles in chain;
# R = radius of outer/inner chain;
# r = radius of circles in outer chain;
#' @export
circles.math.TanToChain = function(n, R, r, type = c("Inner", "Outer")) {
	# Note: assumes that n circles of radius r form a chain of radius R;
	type = match.arg(type);
	fR = R / r; fR2 = fR^2 - 1;
	dR = R^2 - r^2;
	b1 = r + fR * sqrt(dR);
	dd = b1^2 - dR * fR2;
	dd = sqrt(dd);
	if(type == "Outer") dd = - dd;
	r2 = (b1 - dd) / fR2;
	R2 = fR * r2;
	return(list(R = R2, r = r2));
}
#' @export
circles.TanToChain = function(n, r, center = c(0,0), phi = 0,
		type = c("Inner", "Outer")) {
	# TODO: type of baseline-chain;
	c1 = circles.OnFixedCircle(n=n, r=r, center=center, phi=phi);
	R1 = r; r1 = attr(c1, "r");
	R2 = circles.math.TanToChain(n=n, R=R1, r=r1, type=type);
	c2 = pointsCircle(n, r = R2$R, center=center, phi = phi + pi/n);
	attr(c2, "R") = R2$R;
	attr(c2, "r") = R2$r;
	lst = list(C1 = c1, C2 = c2);
	return(lst);
}

# - Generates the actual shape;
#' @export
circles.TanToChainShape = function(n, r, center = c(0,0), phi = 0,
		col = NULL, fill = NULL, type = c("Inner", "Outer")) {
	lst = circles.TanToChain(n=n, r=r, center=center, phi=phi, type=type);
	cc1 = cbind(lst$C1$x, lst$C1$y);
	cc2 = cbind(lst$C2$x, lst$C2$y);
	c1 = list(r = attr(lst$C1, "r"), center = cc1, phi = c(0, 2*pi));
	c2 = list(r = attr(lst$C2, "r"), center = cc2, phi = c(0, 2*pi));
	if( ! is.null(col)) {
		c1$col = col[[1]];
		c2$col = if(length(col) > 1) col[[2]] else col[[1]];
	}
	if( ! is.null(fill)) {
		c1$fill = fill[[1]];
		c2$fill = if(length(fill) > 1) fill[[2]] else fill[[1]];
	}
	lst = list(as.circle(c1), as.circle(c2));
	invisible(as.bioshape(lst));
}


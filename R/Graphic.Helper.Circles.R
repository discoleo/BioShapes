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


### Helper Functions to Generate Circles


### Plot Objects formed from circles;
# - convenience function;
#' @export
testFilledCircle = function(xy, r=NULL, R=NULL, lim=NULL, line=TRUE,
                            col="#B0B032", col.line="green", add=FALSE, pin = FALSE, ...) {
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
    plot(x, y, xlim = lim + mid[1], ylim = lim + mid[2], asp = 1);
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


#####################

### Helper Functions

### Points on a Circle
# r = radius;
# phi = rotation (counter-clockwise);
#' @export
pointsCircle = function(n, r, center = c(0,0), phi=0) {
  x = r * cos(seq(0, n-1) * 2*pi/n + phi) + center[1];
  y = r * sin(seq(0, n-1) * 2*pi/n + phi) + center[2];
  lst = list(x=x, y=y);
  attr(lst, "R") = r;
  attr(lst, "center") = center;
  return(lst);
}


################

### Constructors

#' @export
circle = function(...) {
	UseMethod("circle");
	# NO default method
}


### Circle through 3 Points
#' @export
circle.p3 = function(x, y) {
	if(length(x) != 3 || length(y) != 3) {
		stop("Circle needs 3 points!");
	}
	mid = center.p3(x, y);
	r   = dist.xy(c(x[1], mid[1]), c(y[1], mid[2]));
	lst = list(r=r, center = mid);
	class(lst) = c("circle", "list");
	return(lst);
}

### Tangent & 2 Points:
# - tangent in p1 at line of given slope,
#   and passing also through p2;
#' @export
circle.p2s = function(p1, p2, slope) {
	# y1 = - x1/slope + y0;
	# (x1-x2)*xc + (y1-y2)*yc = (x1^2 - x2^2 + y1^2 - y2^2)/2;
	x1 = p1[1]; y1 = p1[2];
	x2 = p2[1]; y2 = p2[2];
	dy = y1 - y2;
	sinv = 1 / slope; # as (+) value;
	xdiv = (x1 - x2) - sinv*dy;
	xc = (x1^2 - x2^2 + dy*(y1 + y2))/2 - dy*(y1 + x1*sinv);
	xc = xc / xdiv;
	yc = y1 - sinv*(xc - x1);
	r  = dist.xy(c(x1,xc), c(y1,yc));
	lst = list(r=r, center = c(xc, yc));
	class(lst) = c("circle", "list");
	return(lst);
}


### Center of Triangle
#' @export
center.p3 = function(x, y) {
	# Special Cases: slope == Inf or 0;
	if(x[1] == x[3]) {
		if(x[1] == x[2]) stop("Points are collinear!");
		x = x[c(1,3,2)]; y = y[c(1,3,2)];
	} else if(x[2] == x[3]) { x = x[c(2,3,1)]; y = y[c(2,3,1)]; }
	#
	if(x[1] == x[2]) {
		# V Line
		mid12.x = x[1];
		mid12.y = (y[1] + y[2])/2;
		mid13.x = (x[1] + x[3])/2;
		mid13.y = (y[1] + y[3])/2;
		slope = slope(x[c(1,3)], y[c(1,3)]);
		if(slope == 0) {
			return(c(mid13.x, mid12.y));
		} else {
			xc = mid13.x + (mid13.y - mid12.y) * slope;
			return(c(xc, mid12.y));
		}
	}
	if(y[1] == y[3]) {
		if(y[1] == y[2]) stop("Points are collinear!");
		x = x[c(1,3,2)]; y = y[c(1,3,2)];
	}
	if(y[1] == y[2]) {
		# H Line
		mid12.x = (x[1] + x[2])/2;
		mid12.y = y[1];
		mid13.x = (x[1] + x[3])/2;
		mid13.y = (y[1] + y[3])/2;
		slope = slope(x[c(1,3)], y[c(1,3)]);
		yc = mid13.y - (mid12.x - mid13.x)/slope;
		return(c(mid12.x, yc));
	}
	# General Case:
	slope2  = slope(x[c(1,2)], y[c(1,2)]);
	slope3  = slope(x[c(1,3)], y[c(1,3)]);
	# Collinearity check;
	if(slope2 == slope3) stop("Points are collinear!");
	mid12.x = (x[1] + x[2])/2;
	mid12.y = (y[1] + y[2])/2;
	mid13.x = (x[1] + x[3])/2;
	mid13.y = (y[1] + y[3])/2;
	# mid12.y - yc = xc/slope2 - mid12.x/slope2;
	# mid13.y - yc = xc/slope3 - mid13.x/slope3;
	xc = (mid12.y - mid13.y) + mid12.x/slope2 - mid13.x/slope3;
	xc = xc / (1/slope2 - 1/slope3);
	yc = mid12.y + (mid12.x - xc) / slope2;
	return(c(xc, yc));
}


#####################

### Chains of Circles
# - Tangent circles forming a large circle;

### Large Circle of Unknown Radius
# - n small circles each of radius r;
#' @export
circlesOnCircle = function(n, r, center = c(0,0), phi=0) {
  R  = r / sin(pi/n);
  xy = pointsCircle(n, r=R, center=center, phi=phi);
  attr(xy, "R") = R;
  attr(xy, "r") = r;
  return(xy);
}

### Outside a large circle of given radius
#' @export
circlesOutsideFixedCircle = function(n, r, center = c(0,0), phi=0) {
  r1 = r / (1/sin(pi/n) - 1);
  R1 = r + r1;
  xy = pointsCircle(n, r=R1, center=center, phi=phi);
  attr(xy, "R") = R1;
  attr(xy, "r") = r1;
  return(xy);
}

### Forming large circle of given radius
#' @export
circlesOnFixedCircle = function(n, r, center = c(0,0), phi=0) {
  r1 = r * sin(pi/n);
  xy = pointsCircle(n, r=r, center=center, phi=phi);
  attr(xy, "R") = r1 + r; # reuse same attribute name ???
  attr(xy, "r") = r1;
  return(xy);
}

### Inside a large circle of given radius
#' @export
circlesInFixedCircle = function(n, r, center = c(0,0), phi=0) {
  r1 = r / (1 + 1/sin(pi/n));
  R  = r - r1;
  xy = pointsCircle(n, r=R, center=center, phi=phi);
  attr(xy, "R") = R;
  attr(xy, "r") = r1;
  return(xy);
}

###############

### Arcs

# Generates the points on an arc of circle.
#' @export
circle.arc = function(r=1, phi, center = c(0,0), N = 64) {
  id = seq(0, N) / N;
  dp = phi[2] - phi[1];
  dw = id * dp + phi[1];
  x  = r * cos(dw) + center[1];
  y  = r * sin(dw) + center[2];
  xy = list(x=x, y=y);
  return(xy);
}


### Circle arc through p1 & p2
# d = distance from the p1-p2 segment;
# x, y = corresponding x & y coordinates of p1 and p2;
#' @export
circle.ArcByDist = function(x, y, d, col=NULL, fill=NULL, lwd=1, tol=1E-8) {
	L  = dist.xy(x, y);
	sg = sign(d); d = abs(d);
	qq = which.quadrant(x, y);
	dd = L - 2*d;
	center.x = (x[1] + x[2]) / 2;
	center.y = (y[1] + y[2]) / 2;
	center = c(center.x, center.y);
	slope  = slope(x, y);
	if(abs(dd) < tol) {
		r  = L / 2;
		th = pi / 2;
	} else {
		r  = L^2/(8*d) + d/2;
		di = r - d;
		di = sg * di;
		if(qq == 2) {
			th = atan2( - L/2, - di);
		} else if(qq == 4) {
			th = atan2( - L/2, - di);
		} else {
			th = atan2(L/2, di);
		}
		center = shift.ortho(center.x, center.y, d = - di, slope = slope);
		center = unlist(center[1, 1:2]);
	}
	phi = pi/2 + atan(slope) - th;
	# print(c(th=th, phi=phi)/pi);
	if(sg >= 0) {
		phi = c(phi, phi + 2*th);
	} else {
		phi = c(phi + 2*th, phi);
	}
	phi = as.radians0(phi);
	arc = list(center=center, r=r, phi=phi, lwd=lwd);
	if( ! is.null(col))  arc$col = col;
	if( ! is.null(fill)) arc$fill = fill;
	class(arc) = c("circle.arc", "list");
	arc = list(arc);
	return(as.bioshape(arc));
}


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


### Helper Functions to Generate Circles


#' @export
as.circle = function(x) {
	if( ! inherits(x, "circle")) {
		class(x) = c("circle", class(x));
	}
	return(x);
}


#####################

### Helper Functions

### Points on a Circle
# aka Regular n-Gon;
# r = radius;
# phi = rotation (counter-clockwise);
#' @export
pointsCircle = function(n, r, center = c(0,0), phi = 0, clock = FALSE) {
	pp = points.circle(n=n, r=r, center=center, phi=phi, clock=clock);
	return(pp);
}
#' @export
points.circle = function(n, r, center = c(0,0), phi = 0, clock = FALSE) {
  x = r * cos(seq(0, n-1) * 2*pi/n + phi) + center[1];
  y = r * sin(seq(0, n-1) * 2*pi/n + phi) + center[2];
  if(clock) {
	x = rev(x); y = rev(y);
  }
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


### Ring/Annulus
# w   = width of annulus;
# lwd = width of border;
#' @export
circle.ring = function(xy, r = 1, w = 0.25, col = c("#D0D048", "#969612"), lwd = 4, w.scale = 68) {
	d = w.scale * w;
	lst = list(r=r, center = xy, lwd = d, col = col[1]);
	lst = list(An = as.circle(lst));
	brd = list(r = r + w, center = xy, lwd=lwd, col = col[2]);
	lst$Bout = as.circle(brd);
	brd = list(r = r - w, center = xy, lwd=lwd, col = col[2]);
	lst$Bin = as.circle(brd);
	invisible(as.bioshape(lst));
}

###############
###############

############
### Arcs ###

#' @export
arc = function(x, y, center = c(0,0), ...) {
	UseMethod("arc");
}

### Arc: Start & Stop-Angles
# Note:
# - does NOT need r, as it does NOT verify
#   that the points are on a circle of given radius!
#' @export
arc.circle = function(x, y, center = c(0, 0)) {
	dx = x - center[1];
	dy = y - center[2];
	atan2(dy, dx);
}

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

#' @export
circle.ArcByDiam = function(x, y, top = TRUE) {
		r  = dist.xy(x, y) / 2;
		cx = (x[1] + x[2]) / 2;
		cy = (y[1] + y[2]) / 2;
		phi = atan2(y[2] - y[1], x[2] - x[1]);
		if( ! top) phi = as.radians0(phi + pi);
		lst = list(r=r, center = c(cx, cy), phi = c(phi, phi + pi));
		class(lst) = c("circle.arc", "list");
		return(lst);
}
#' @export
circle.ArcByDiam.xy = function(p1, p2, top = TRUE) {
	circle.ArcByDiam(c(p1[1], p2[1]), c(p1[2], p2[2]), top=top);
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


##################

### Fill Circle

### Uniform inside Circle
# Note: NOT random;
#' @export
uniform.circle = function(n, r = 1, center = c(0,0), phi = 0, d = 0.5, d.sep = 5) {
	id = seq(0, n) + d;
	rr = sqrt(id/n) * r;
	th = pi * (1 + sqrt(d.sep)) * id + phi;
	xy = cbind(rr * cos(th), rr * sin(th));
	xy[,1] = xy[,1] + center[1];
	xy[,2] = xy[,2] + center[2];
	return(xy);
}


### Hashed Circle
# n = number of lines;
# phi = counter-clockwise rotation starting at pi/2;
#' @export
circle.hash = function(n, center = c(0, 0), r = 1, phi = 0, scale = 1,
		lwd = 1, lty = 1, col = NULL, neps = 1/(2*n)) {
	n0 = c(-1, 1) + if(length(neps) == 1) c(neps, - neps) else neps;;
	tx = seq(n0[1], n0[2], length.out = n);
	ty = sqrt(1 - tx^2);
	tx = r * tx; ty = r * ty;
	cs = cos(phi); sn = sin(phi);
	x1 = tx * cs + ty * sn + center[1];
	y1 = (tx * sn - ty * cs) * scale + center[2];
	x2 = tx * cs - ty * sn + center[1];
	y2 = (tx * sn + ty * cs) * scale + center[2];
	hasCol = ! is.null(col);
	if(hasCol && length(col) == 1) col = rep(col, n);
	xy = lapply(seq(n), function(id) {
		lst = list(x = c(x1[id], x2[id]), y = c(y1[id], y2[id]),
			lwd=lwd, lty=lty);
		if(hasCol) lst$col = col[id];
		return(lst);
	});
	xy = as.bioshape(xy);
	return(xy);
}


### Simulate Emboss
# - but larger lwd.bg causes line to extend beyond circle boundary:
#   is NOT a good replacement for the blur/diffusion operation;
#' @export
circle.hashEmboss = function(n, center = c(0, 0), r = 1, phi = 0, scale = 1,
		d = - lwd/(3*n), lwd = 1, lwd.bg = 1.5 * lwd, lty = 1,
		col = NULL, col.bg = "#A0A0A0") {
	# cBG = circle.hash(n=n, center=center, r=r, phi=phi, scale=scale,
	#	lwd = lwd.bg, col = col.bg, neps = c(1/2 + 1/3, - 1/2 + 1/3) / n);
	lst = circle.hash(n=n, center=center, r=r, phi=phi, scale=scale,
		lwd=lwd, lty=lty, col=col);
	slope = tan(phi + pi/2);
	cBG = lapply(seq(n), function(id) {
		ln = lst[[id]];
		pp = shift.ortho(ln$x, ln$y, slope=slope, d=d, scale=scale, id.offset=id);
	})
	cBG = do.call(rbind, cBG);
	cBG = list(cBG, lwd = lwd.bg, col = col.bg);
	cBG = as.bioshape(cBG);
	lst = list(BG = cBG, Lines = lst);
	return(as.bioshape(lst));
}


#########################

### Ellipse

# - returns solution closest to "y.guess";
#' @export
solve.ellipse.closest = function(x, y.guess = NULL, r, phi = 0, center = c(0,0), tol = 1E-10) {
	cc = as.coeff.ellipse(r, phi=phi);
	x0 = center[1]; y0 = center[2];
	xx = x - x0;
	y  = solve.ellipse.v0(xx, coef = cc, tol=tol);
	y  = y + y0;
	# NO Guess:
	if(is.null(y.guess) || length(y) == 1) return(y);
	# Solution closest to "y.guess"
	# Note: solution must be valid!
	dG = abs(y - y.guess);
	y  = ifelse(dG[1] <= dG[2], y[1], y[2]);
}

#' @export
solve.ellipse.all = function(x, r, phi = 0, center = c(0,0), tol = 1E-10) {
	x0 = center[1]; y0 = center[2];
	xx = x - x0;;
	cc = as.coeff.ellipse(r, phi=phi);
	y = sapply(xx, function(xx) {
		y = solve.ellipse.v0(xx, coef = cc, tol=tol);
		if(length(y) == 1) y = c(y, y);
		return(y);
	});
	y = y + y0;
	return(y);
}

### Basic Solver
#' @export
solve.ellipse.v0 = function(x, coef, tol = 1E-10) {
	cc = coef;
	cc[2] = cc[2]*x;
	cc[3] = cc[3]*x^2 - 1;
	Delta = cc[2]^2 - 4*cc[1]*cc[3];
	if(abs(Delta) <= tol) {
		y = - cc[2] / (2*cc[1]);
		return(y);
	}
	if(Delta < 0) return(NA);
	Delta = sqrt(Delta);
	y = (- cc[2] + c(-1, 1)*Delta) / (2*cc[1]);
	return(y);
}

### Arc Angles
# - computes Start & Stop-Angles;
#' @export
arc.ellipse = function(x, y, r = c(1,2), center = c(0,0), theta = 0) {
	x0 = x - center[1];
	y0 = y - center[2];
	sn = sin(theta); cs = cos(theta);
	B  = matrix(c(r[1]*cs, r[1]*sn, - r[2]*sn, r[2]*cs), nrow=2);
	trig.phi = solve(B, rbind(x0, y0));
	atan2(trig.phi[2,], trig.phi[1,]);
}


### Ellipse:
# - passing through (x1, y1) & (x2, y2);
# - with slope = slope(x, y);
# - center = mid(xy) + d;
# - rr = ratio of 2 radii;
#' @export
solve.ellipse.bySlope = function(x, y, d = 1, rr = 2) {
	sl = slope(x, y);
	mid.x = (x[1] + x[2]) / 2;
	mid.y = (y[1] + y[2]) / 2;
	xy = shift.ortho(mid.x, mid.y, slope = sl, d=d);
	xy = unlist(xy[1, c(1,2)]);
	cc = as.coef.ellipse.slope(x=x, y=y, slope = sl, center = xy);
	cc = cc[1,];
	a2 = cc[1] + cc[2] * rr^2;
	a = sqrt(a2);
	b = a / rr;
	return(list(r = c(a, b), center = xy, th = atan(sl)));
}
#' @export
as.coef.ellipse.slope = function(x, y, slope, center = c(0,0)) {
	sq = sqrt(slope^2 + 1);
	sn = 1 / sqrt(1 + 1/slope^2);
	if(slope < 0) sn = - sn;
	cs = 1 / sq;
	x = x - center[1];
	y = y - center[2];
	coeff = cbind(
		(cs*x + sn*y)^2,
		(sn*x - cs*y)^2);
	return(coeff);
}

### Coefficients
# Eq: c[1]*(y - y0)^2 + c[2]*(x - x0)*(y - y0) + c[3]*(x - x0)^2 - 1 = 0;
# - Note: c[] does NOT include: - 1;
#' @export
as.coeff.ellipse = function(r, phi = 0) {
	sn = sin(phi); cs = cos(phi);
	a = r[1]; b = r[2];
	cc = c((sn^2 / a^2 + cs^2 / b^2),
		(1/a^2 - 1/b^2)*sin(2*phi), (cs^2 / a^2 + sn^2 / b^2));
	return(cc);
}


### Ellipse: x-Range
#' @export
range.ellipse.x = function(r, phi = 0, center = c(0,0)) {
	cc = as.coeff.ellipse(r, phi=phi);
	# 2*c1*y + c2*x = 0
	# (4*c1*c3 - c2^2)*x^2 - 4*c1 # = 0
	x = 4*cc[1] / (4*cc[1]*cc[3] - cc[2]^2);
	# TODO: check explicit sin/cos formula;
	if(x < 0) warning("Invalid solution!");
	x = sqrt(x);
	x = center[1] + c(-x, x);
	return(x);
}

### Slope of Tangents to Ellipse
# (x, y) = given points on ellipse;
#' @export
slope.ellipse = function(x, y, r, phi = 0, center = c(0,0)) {
	# Note: does NOT check validity of point (x, y);
	cc = as.coeff.ellipse(r, phi=phi);
	xx = x - center[1];
	yy = y - center[2];
	sl = - (2*cc[3]*xx + cc[2]*yy) / (2*cc[1]*yy + cc[2]*xx);
	return(sl);
}

### Points at which Tangent has given slope
#' @export
solve.ellipse.xtan = function(slope, r, phi = 0, center = c(0,0), coeff = NULL) {
	cc = if( ! is.null(coeff)) coeff else as.coeff.ellipse(r, phi=phi);
	# (2*sl*cc[1] + cc[2])*yy = - (2*cc[3] + sl*cc[2])*xx
	if(abs(slope) == Inf) {
		dv = 4*cc[3]*cc[1]^2 - cc[1]*cc[2]^2;
		c4 = 2*cc[1] * sign(slope);
	} else {
		c4 = (2*slope*cc[1] + cc[2]);
		c5 = (2*cc[3] + slope*cc[2]);
		dv = cc[1]*c5^2 + cc[3]*c4^2 - cc[2]*c4*c5;
	}
	if(dv < 0) return(c(NA, NA));
	x = c4 / sqrt(dv);
	x = c(-x, x) + center[1];
	return(x);
}
#' @export
solve.ellipse.xytan = function(slope, r, phi = 0, center = c(0,0)) {
	cc = as.coeff.ellipse(r, phi=phi);
	xx = solve.ellipse.xtan(slope=slope, coeff = cc, center = c(0,0));
	yy = if(abs(slope) == Inf) { - cc[2] * xx / (2*cc[1]); }
		else if((div <- 2*slope*cc[1] + cc[2]) == 0) {
			# print("Div 0!")
			c(-1, 1) / sqrt(cc[1]); # TODO: check!
		} else - (2*cc[3] + slope*cc[2]) * xx / div;
	x  = xx + center[1];
	y  = yy + center[2];
	return(cbind(x, y));
}


### Ellipse tangent to Line
# - Passing through (x1, y1) and (x2, y2);
# - Tangent to Lines;

# TODO


### Intersection: Ellipse w Line
# Line passing through xy with same slope as the ellipse;
# r = radii of ellipse;
#' @export
intersect.ellipse.ParallelLine = function(r, xy, center = c(0,0), phi = 0) {
	sn = sin(phi); cs = cos(phi);
	sl = sn / cs; a = r[1]; b = r[2];
	xc = center[1]; yc = center[2];
	sab = sn^2 / a^2 + cs^2 / b^2;
	dab = 1/a^2 - 1/b^2;
	ssl = sl*(xc - xy[1]) - yc + xy[2];
	sn2 = 2*sn*cs;
	# cc = c((cs^2 / a^2 + sn^2 / b^2) * x^2,
	#	(1/a^2 - 1/b^2)*sin(2*phi), (sn^2 / a^2 + cs^2 / b^2) * y^2);
	cc = c(cs^2 / a^2 + sn^2 / b^2 + sl^2*sab + sl*dab*sn2,
		dab*ssl*sn2 + 2*sl*sab*ssl, sab*ssl^2 - 1);
	Dt = cc[2]^2 - 4*cc[1]*cc[3];
	if(Dt < 0) return(NA);
	Dt = sqrt(Dt);
	x = center[1] - (cc[2] + c(Dt, - Dt)) / (2*cc[1]);
	y = sl*(x - xy[1]) + xy[2];
	return(cbind(x, y));
}

#####################
#####################

### Cylinder

#' @export
cylinder = function(x, y, w = 1, rr = 0.5, lwd = NULL, col = NULL, lty.back = 2) {
	lst = box.capEllipse(x, y, h = w, y.rel = rr, col=col);
	lst$Cap2$phi = c(0, 2*pi);
	if( ! is.null(lty.back)) {
		lst$Cap3 = lst$Cap1;
		lst$Cap3$phi = rev(lst$Cap3$phi);
		lst$Cap3$lty = lty.back;
	}
	return(lst)
}


#' @export
cylinder.bySlope = function(xy, slope = Inf, w = 1, h = 4*w, rr = 0.5,
		lwd = NULL, col = NULL, lty.back = 2) {
	xyE = shift.point(xy, slope = slope, d = h);
	x = c(xy[1], xyE[1]);
	y = c(xy[2], xyE[2]);
	lst = cylinder(x=x, y=y, w=w, rr=rr, lwd=lwd, col=col, lty.back=lty.back);
	return(lst);
}

# n = number of filling circles;
#' @export
cylinder.circleFill = function(x, y, w = 1, rr = 0.5, n = 15,
		lwd = NULL, lwd.fill = lwd, col = NULL, col.circles = col, lty.back = NULL) {
	if(n < 0) stop("Invalid number of circles!");
	lst = cylinder(x=x, y=y, w=w, rr=rr, lwd=lwd, col=col, lty.back=lty.back);
	tt = seq(n) / (n+1);
	cx = (1-tt) * x[1] + tt*x[2];
	cy = (1-tt) * y[1] + tt*y[2];
	tmp = lst$Cap1;
	crc = lapply(seq(n), function(id) {
		tmp$center = c(cx[id], cy[id]);
		tmp$lwd = lwd.fill;
		if( ! is.null(col.circles)) tmp$col = col.circles;
		return(tmp);
	});
	crc = as.bioshape(crc);
	lst = c(Fill = list(crc), lst);
	invisible(as.bioshape(lst));
}

### Multipple Cyclinders
# - parallel cylinders:
#   slope[1] = direction of cylinders;
#   slope[2] = direction of ensemble;
#' @export
cylinder.tubes = function(xy, n, w = 1, h = 4*w, d = 0.25, slope = c(Inf, 0), col = NULL, ...) {
	if(n < 1) stop("Invalid number of domains!");
	if(n == 1) return(cylinder.bySlope(xy, w=w, h=h, slope = slope[1], ...));
	#
	wd  = w + d;
	dH  = (n - 1) * wd;
	xyE = shift.point(xy, slope = slope[2], d = dH);
	xyr = shift.point(xy,  slope = slope[1], d = h);
	xys = shift.point(xyE, slope = slope[1], d = h);
	if(length(col) == 1) col = rep(col, n);
	n = n - 1;
	lst = lapply(seq(0, n), function(id) {
		tt  = wd * id / dH;
		xyS = (1 - tt)*xy  + tt*xyE;
		xyT = (1 - tt)*xyr + tt*xys;
		cylinder(c(xyS[1], xyT[1]), c(xyS[2], xyT[2]), w=w, col = col[id + 1], ...);
	});
	invisible(as.bioshape(lst));
}



### Ellipse

### Constructor: Trivial
#' @export
ellipse = function(x, y, r, theta = 0, phi = c(0, 2*pi), lwd = NULL, col = NULL, fill = NULL) {
	if(missing(y)) {
		if(! inherits(x, "matrix")) x = matrix(x, ncol = 2);
		xy = xy.coords(x);
	} else xy = xy.coords(x, y);
	len = length(xy$x);
	if(len == 1) {
		center = c(xy$x, xy$y);
	} else {
		center = cbind(xy$x, xy$y);
	}
	lst = list(center=center, r=r, theta=theta, phi=phi);
	if( ! is.null(fill)) lst$fill = fill;
	if( ! is.null( col)) lst$col  = col;
	if( ! is.null( lwd)) lst$lwd  = lwd;
	lst = as.ellipse(lst);
	return(lst);
}

### Ellipse on Helix Wheel
# - Helix Wheel: in Graphic.Chemistry.Helix.R;
#' @export
ellipse.HelixWheel = function(id, data, labels = NULL,
			phi = c(0, pi), col = 1, lwd = 2,
			r.scale = 3.25, dy.scale = 1, cex = NULL) {
	xy = data$C$center[id, ];
	cx = mean(xy[, 1]);
	cy = mean(xy[, 2]);
	rr = data$C$r[[1]];
	iL = nrow(xy);
	th = atan2(xy[iL,2] - xy[1,2], xy[iL,1] - xy[1,1]);
	d  = dist.xy(xy[,1], xy[,2]) / 2;
	ee = ellipse(cx, cy, r = c(d + rr, rr*r.scale),
			theta = th, phi=phi, col=col, lwd=lwd);
	lst = list(E = ee);
	# Label:
	if( ! is.null(labels)) {
		thpi = th / pi;
		sg = if(thpi <= -3/4) - dy.scale else dy.scale;
		# print(th / pi * 180);
		txt = list(x = cx, y = cy + 1.5625 * sg*rr*r.scale,
			labels=labels, col=col);
		if( ! is.null(cex)) txt$cex = cex;
		lst$Lbl = txt;
	}
	return(lst);
}

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


###################
#
# Title: BioShapes
#
# Maintainer: L. Mada
#
#
# Continues the work of:
# - Adrian Cotoc (BSc 2022-2023)
# - Darian Voda (BSc 2021-2022)
#
# GitHub: https://github.com/discoleo/BioShapes


### Arrows: Other Types

# - Pins/Tags;
# - Boxes & Labels;


#######################

### Pin/Tag
# (x, y) = from (pin-point, center);
# theta  = angle of pin;
#' @export
pin.center = function(x, y, theta = pi/3, lwd=1, lwd.circle = lwd,
		col = 1, col.circle = col, fill = NULL, debug = FALSE) {
	th = theta / 2;
	d = dist.xy(x, y);
	r = d * sin(th);
	# Circle Arc
	slope = slope(x, y);
	phi = atan2(y[2] - y[1], x[2] - x[1]);
	phi = - pi/2 + (phi - th);
	phi = c(phi, phi - (pi - theta));
	if(debug) print(phi/pi);
	lst = list(r = r, center = c(x[2], y[2]), phi = phi,
		lwd = lwd.circle, col = col.circle);
	if( ! is.null(fill)) lst$fill = fill;
	class(lst) = c("circle.arc", "list");
	# Triangle:
	tt = cos(th)^2;
	xx = tt*x[2] + (1-tt)*x[1];
	yy = tt*y[2] + (1-tt)*y[1];
	hT = d * sin(theta) / 2;
	pp = shift.ortho(c(xx, yy), slope = slope(x, y), d = c(-hT, hT));
	xx = c(pp[1,1], x[1], pp[2,1]);
	yy = c(pp[1,2], y[1], pp[2,2]);
	lstT = list(x = xx, y = yy, lwd=lwd, col=col);
	if( ! is.null(fill)) {
		lstP = lstT;
		lstP$col = NA; # remove border;
		lstP$fill = fill;
		class(lstP) = c("polygon", "list");
		lstT = list(Poly = lstP, Triangle = lstT);
		lstT = as.bioshape(lstT);
	}
	lst0 = list(Triangle = lstT, Arc = lst);
	return(as.bioshape(lst0));
}

######################
######################

### Boxes & Labels


# TODO:
# - correct length by - h;
# - correct/verify arcs when scale != 1;
#' @export
box.cap = function(x, y, h = 1, lwd = 1, col = NULL, fill = NULL, scale = 1) {
	slope = slope(x, y);
	h2 = h/2;
	xy = shift.ortho(x, y, slope=slope, d = c(-h2, h2), scale=scale);
	qd = which.quadrant(x, y);
	lst  = list(Poly = xy, lwd = lwd);
	### Caps
	phi0 = atan(slope / scale);
	if(qd == 2 || qd == 4) phi0 = pi/2 - phi0;
	if(slope == -Inf) phi0 = pi/2;
	phiA = c(phi0 + pi/2, phi0 + 3*pi/2);
	phiB = c(phi0 + 3*pi/2, phi0 + pi/2);
	# 2 | 1
	# 3 | 4
	if(qd == 1 || qd == 2) {
		phi1 = phiA; phi2 = phiB;
	} else {
		phi1 = phiB; phi2 = phiA;
	}
	cap1 = list(r = h2, center = c(x[1], y[1]), phi = phi1, lwd=lwd);
	cap2 = list(r = h2, center = c(x[2], y[2]), phi = phi2, lwd=lwd);
	if( ! is.null(col)) {
		lst$col = col;
		cap1$col = col;
		cap2$col = col;
	}
	if( ! is.null(fill)) {
		cap1$fill = fill;
		cap2$fill = fill;
		# Polygon:
		pp = list(
			x = lst$Poly$x[c(1,2,4,3)],
			y = lst$Poly$y[c(1,2,4,3)]);
		pp$fill = fill; pp$col = NA;
		class(pp) = c("polygon", "list");
		lst = c(PolyFill = list(pp), lst);
	}
	class(cap1) = c("circle.arc", "list");
	class(cap2) = c("circle.arc", "list");
	cap = as.bioshape(list(Cap1 = cap1, Cap2 = cap2))
	lst = c(lst, Cap = cap);
	return(as.bioshape(lst));
}


#' @export
box.capEllipse = function(x, y, h = 1, y.rel = 0.25, lwd = 1, col = NULL, fill = NULL, scale = 1) {
	slope = slope(x, y);
	h2 = h/2;
	xy = shift.ortho(x, y, slope=slope, d = c(-h2, h2), scale=scale);
	lst  = list(Poly = xy, lwd = lwd);
	### Caps
	isV  = (x[2] == x[1]);
	phi0 = if(isV) atan2(y[2] - y[1], 0)
		else atan2(y[2] - y[1], (x[2] - x[1]) * scale);
	pi2  = pi / 2;
	if(slope == -Inf) phi0 = pi2;
	th = phi0 + pi2;
	phi1 = c(0, pi); phi2 = c(pi, 2*pi);
	# BUG: in shape::getellipse ???
	if(isV && y[2] < y[1]) {
		tmp = phi1; phi1 = phi2; phi2 = tmp;
	}
	r = c(h2, y.rel * h2 / scale);
	cap1 = list(r=r, center = c(x[1], y[1]), th = th, phi = phi1, lwd=lwd, scale=scale);
	cap2 = list(r=r, center = c(x[2], y[2]), th = th, phi = phi2, lwd=lwd, scale=scale);
	if( ! is.null(col)) {
		lst$col = col;
		cap1$col = col;
		cap2$col = col;
	}
	if( ! is.null(fill)) {
		cap1$fill = fill;
		cap2$fill = fill;
		# Polygon:
		pp = list(
			x = lst$Poly$x[c(1,2,4,3)],
			y = lst$Poly$y[c(1,2,4,3)]);
		pp$fill = fill; pp$col = NA;
		class(pp) = c("polygon", "list");
		lst = c(PolyFill = list(pp), lst);
	}
	class(cap1) = c("ellipse", "list");
	class(cap2) = c("ellipse", "list");
	cap = as.bioshape(list(Cap1 = cap1, Cap2 = cap2))
	lst = c(lst, Cap = cap);
	return(as.bioshape(lst));
}


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

# - Circular Arrows;
# - Pins/Tags;
# - Boxes & Labels;


###################

### Circular Arrows
# r = radius of circle;
# w = width of Annulus;
# type = Types of circular arrow;
# - Equal: Head = Symmetric/Isosceles triangle;
# - OutEqual: as Equal, but tip extends beyond phi[2];
# - Gene: End has no arrow;
# - ShGene: simpler, but tip is shifted above midline circle;
#' @export
arrow.circular = function(phi, r = 2, center = c(0,0), w = 0.5,
		type = c("Equal", "Gene", "ShGene", "OutEqual", "Ring", "Ugly"),
		tip.scale = 1, N = NULL, ...) {
	type = match.arg(type);
	w2 = w/2;
	rr = c(r - w2, r + w2);
	#
	c1 = list(r=rr[1], center=center, phi=phi, N=N);
	c1 = as.circle.arc(c1);
	c2 = list(r=rr[2], center=center, phi = rev(phi), N=N);
	c2 = as.circle.arc(c2);
	# Simple Ring segment
	if(type == "Ring") {
		lst = list(Ci = c1, C2 = c2);
		class(lst) = c("polycircle", "list");
		lst = as.bioshape(list(A = lst));
		return(lst);
	}
	# Complex Types:
	isGene = type == "Gene";
	if(isGene) {
		# Experimental: midline tip;
		d = w2 * tip.scale;
		alpha = solve.circular.tip.phi(w2, d=d, r=r, ...);
		id  = 1; # TODO
		az1 = alpha$dphi1[id];
		az2 = alpha$dphi2[id];
		c1$phi[2] = phi[2] - az1;
		c2$phi[1] = phi[2] - az2;
		dphi = 0;
	}
	#
	w2 = w2 * tip.scale;
	isShGene = type == "ShGene";
	if(type == "Ugly)") {
		dphi = 2 * asin(w2 / r);
	} else if( ! isGene) {
		r0 = r;
		r  = sqrt(r^2 + w2^2);
		dphi = atan(w2/r);
		if(type == "Equal") {
			# shifted phi
			phi = phi - dphi;
			c1$phi = phi;
			c2$phi = rev(phi);
		} else if(isShGene) {
			phi[2] = phi[2] - dphi;
			c1$phi[2] = phi[2];
			c2$phi[1] = phi[2];
		}
	}
	if(isShGene || isGene) {
		xyE = NULL;
	} else {
		xE = r*cos(phi[1] + dphi) + center[1];
		yE = r*sin(phi[1] + dphi) + center[2];
		xyE = cbind(xE, yE);
	}
	# rr = rev(rr);
	xT = r*cos(phi[2] + dphi) + center[1];
	yT = r*sin(phi[2] + dphi) + center[2];
	#
	lst = list(Cin = c1, cbind(xT, yT), Cout = c2, xyE);
	class(lst) = c("polycircle", "list");
	lst = as.bioshape(list(A = lst));
	return(lst);
}

# Experimental
solve.circular.tip = function(w, d, r, rm.complex = TRUE,
		as.square = FALSE, debug = TRUE, tol = 1E-8) {
	r2 = r*r;
	w2 = w*w; w4 = w2*w2;
	d2 = d*d; d4 = d2*d2; d6 = d2*d4;
	coeff = c(4*w2*r2*d4, w4*d2 + 4*w2*r2*d2 - 2*w2*d4 - 4*r2*d4 + d6,
		(w4 - 4*w2*d2 - 4*r2*d2 + 3*d4), 3*d2 - 2*w2, 1);
	h2 = polyroot(coeff);
	if(debug) print(h2);
	isCZero = abs(Im(h2)) < tol;
	if(rm.complex) {
		h2 = Re(h2[isCZero]);
		H2 = d2 + h2;
		isReal = H2 >= 0;
		H2 = H2[isReal];
	} else {
		h2[isCZero] = Re(h2[isCZero]);
		H2 = d2 + h2;
	}
	if(as.square) {
		return(H2);
	} else {
		H = sqrt(H2);
		return(H);
	}
}
solve.circular.tip.phi = function(w, d, r, rm.complex = TRUE, debug = TRUE, tol = 1E-8) {
	h2 = solve.circular.tip(w, d, r, as.square = TRUE,
		rm.complex=rm.complex, debug=debug, tol=tol);
	r2 = r*r; rn = r-w; rp = r+w;
	cs1 = (r2 + rn^2 - h2) / (2*r*rn);
	cs2 = (r2 + rp^2 - h2) / (2*r*rp);
	if(rm.complex) {
		isReal = cs1 >= -1 & cs1 <= 1 & cs2 >= -1 & cs2 <= 1;
		if(debug) cat("Is Real: ", isReal, "\n");
		cs1 = cs1[isReal];
		cs2 = cs2[isReal];
	}
	az1 = acos(cs1);
	az2 = acos(cs2);
	if(debug) print(rbind(az1, az2));
	return(list(dphi1 = az1, dphi2 = az2));
}


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
	lst = c(lst, cap);
	return(as.bioshape(lst));
}


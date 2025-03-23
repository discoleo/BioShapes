#######################################
#
# Title: BioShapes
# Maintainer: L. Mada
#
# https://github.com/discoleo/BioShapes


### Grids


### Square Grid

### Circles: Distributed Evenly
#' @export
grid.squareEC = function(x, y, r = 1/2, col = NULL, fill = NULL,
		type = c("square", "interlaced"), d = r/4) {
	type = match.arg(type);
	if(inherits(d, "formula")) {
		d = eval(d[[2]]);
	}
	if(type == "square") {
		cc = grid.squareE.centers(x, y, r=r, d=d);
	} else {
		cc = grid.laced.centers(x, y, r=r, d=d);
	}
	lst = list(center = as.matrix(cc), r=r);
	if(! is.null(col))  lst$col  = col;
	if(! is.null(fill)) lst$fill = fill;
	clShape = if(length(r) == 1) "circle" else "ellipse";
	class(lst) = c(clShape, "list");
	lst = list(G = lst);
	return(as.bioshape(lst));
}

### Circles: Distributed Evenly + Random Effects
#' @export
grid.squareREC = function(x, y, r = 1/2, phi = 0,
		col = NULL, fill = NULL, type = c("interlaced", "square"),
		d = 3/2*r, scale.center = 1/4, scale.r = 1/6, scale.phi = 1) {
	type = match.arg(type);
	if(inherits(d, "formula")) {
		d = eval(d[[2]]);
	}
	if(length(r) == 1) r = c(r,r);
	if(type == "square") {
		cc = grid.squareE.centers(x, y, r=r, d=d);
	} else {
		cc = grid.laced.centers(x, y, r=r, d=d);
	}
	# Random effects:
	n = nrow(cc); cc.sc = r * scale.center;
	cc$x  = cc$x + runif(n, - cc.sc[1], cc.sc[1]);
	cc$y  = cc$y + runif(n, - cc.sc[2], cc.sc[2]);
	r.sc  = r * scale.r;
	cc$rx  = runif(n, r[1] - r.sc[1], r[1] + r.sc[1]);
	cc$ry  = runif(n, r[2] - r.sc[2], r[2] + r.sc[2]);
	phi.sc = pi/2 * scale.phi;
	cc$phi = phi + runif(n, - phi.sc, phi.sc);
	# Object:
	lst = list(center = cc);
	if(! is.null(col))  lst$col  = col;
	if(! is.null(fill)) lst$fill = fill;
	clShape = if(length(r) == 1) "circle" else "ellipse";
	class(lst) = c(clShape, "list");
	lst = list(G = lst);
	return(as.bioshape(lst));
}

### Distributed Evenly
# r = radius of shapes;
# d = distance between individual shapes;
# (x, y) = pair of coordinates of rectangle;
#' @export
grid.squareE.centers = function(x, y, r = 1/2, d = r/4) {
	if(length(r) == 1) r = c(r, r);
	if(length(d) == 1) d = c(d, d);
	div = 2*r + d;
	nx = floor((abs(x[1] - x[2]) + d[1]) / div[1]);
	ny = floor((abs(y[1] - y[2]) + d[2]) / div[2]);
	sx = if(x[2] < x[1]) -1 else 1;
	sy = if(y[2] < y[1]) -1 else 1;
	cx = seq(x[1] + sx*r[1], x[2] - sx*r[1], length.out = nx);
	cy = seq(y[1] + sy*r[2], y[2] - sy*r[2], length.out = ny);
	cc = expand.grid(cx, cy);
	names(cc) = c("x", "y")
	return(cc);
}

# Interlaced Rows:
#' @export
grid.laced.centers = function(x, y, r = 1/2, d = r/2) {
	if(length(r) == 1) r = c(r, r);
	if(length(d) == 1) d = c(d, d);
	div = 2*r + d;
	nx = floor((abs(x[1] - x[2]) + d[1]) / div[1]);
	ny = floor((abs(y[1] - y[2]) + d[2]) / div[2]);
	sx = if(x[2] < x[1]) -1 else 1;
	sy = if(y[2] < y[1]) -1 else 1;
	cx = seq(x[1] + sx*r[1], x[2] - sx*r[1], length.out = nx);
	cy = seq(y[1] + sy*r[2], y[2] - sy*r[2], length.out = ny);
	cc = expand.grid(cx, cy);
	names(cc) = c("x", "y");
	if(ny <= 1) return(cc);
	# Interlaced rows:
	cx = (cx[ - nx] + cx[-1]) / 2;
	cy = (cy[ - ny] + cy[-1]) / 2;
	ci = expand.grid(cx, cy);
	names(ci) = c("x", "y");
	cc = rbind(cc, ci);
	return(cc);
}


### Hexagonal Grid

### Hexagonal Grid: Circular
# n = Number of levels;
# r = Radius of each hexagon;
# r.hexa = actual radius of each hexagon (as drawn);
# phi = Rotation of the grid/hexagons;
#' @export
grid.hexa = function(n, r = 1, center = c(1,1), col = NULL,
		phi = 0, r.hexa = r) {
	r0 = r; r = r.hexa;
	if(n < 0) stop("The number of levels n must be positive!");
	if(n == 0) return(as.bioshape.empty());
	if(is.null(col)) {
		# col = # TODO
	} else if(length(col) == 1) { col = rep(col, n); }
	else if(length(col) < n) {
		warning("Vector of colours is shorter!");
	}
	### Grid:
	L = ngon(6, center=center, r=r, col = col[1], phi = phi + pi/6);
	G = list(L1 = L);
	if(n == 1) return(as.bioshape(G));
	# L2:
	L = ring.hexa(6, center=center, r=r, phi=phi,
		R = sqrt(3) * r0, col = col[2]);
	G = c(G, L2 = as.bioshape(L));
	if(n == 2) return(as.bioshape(rev(G)));
	# L 2:n
	L = lapply(seq(3, n), function(nL) {
		# Simpler version: only 6 outer hexagons;
		R = r0 * sqrt(3) * (nL-1);
		cc = grid.centers.hexaSimple(n = nL-1, r = R, center=center, phi=phi);
		Ln = lapply(seq(nrow(cc)), function(id) {
			ngon(6, center = unlist(cc[id,]), r=r,
				col=col[nL], phi = phi + pi/6);
		});
		return(as.bioshape(Ln));
	});
	names(L) = paste0("L", seq(3, n));
	G = rev(c(G, L));
	return(invisible(as.bioshape(G)));
}

### Ring of Hexagons
# Note: Helper function;
# n = Number of Hexagons;
# p = Number of partitions between 2 outer hexagons;
# R = Radius of ring;
# id.filter = Which hexagons to keep from the n-ring;
#' @export
grid.hexa.ring = function(n, p, R, r = 1, center = c(1,1),
		id.filter = NULL, col = NULL, phi = 0, phi.h = pi/6) {
	cc = grid.hexa.centers(n, p=p, r = R, center=center,
		id.filter=id.filter, phi=phi);
	if(nrow(cc) == 0) return(NULL);
	lst = lapply(seq(nrow(cc)), function(id) {
		ngon(6, center = c(cc$x[id], cc$y[id]), r=r, col=col, phi=phi.h);
	});
	lst = as.bioshape(lst);
	invisible(lst);
}
# Centers of the Hexagon-Ring
#' @export
grid.centers.hexaSimple = function(n, r, center = c(1,1), phi = 0) {
	cc = ngon(6, r=r, center=center, phi=phi);
	cc$x = c(cc$x, cc$x[1]);
	cc$y = c(cc$y, cc$y[1]);
	# Partition Lines between the 6 Points/Centers:
	tt = seq(0, n-1) / n;
	partf = function(id) {
		x1 = cc$x[id]; x2 = cc$x[id+1];
		y1 = cc$y[id]; y2 = cc$y[id+1];
		x  = (1-tt)*x1 + tt*x2;
		y  = (1-tt)*y1 + tt*y2;
		return(data.frame(x=x, y=y));
	}
	cc = lapply(seq(6), partf);
	cc = do.call(rbind, cc);
	return(cc);
}
#' @export
grid.hexa.centers = function(n, p, r, center = c(1,1),
		id.filter = NULL, phi = 0) {
	cc = ngon(n, r=r, center=center, phi=phi);
	if( ! is.null(id.filter)) {
		cc$x = cc$x[id.filter];
		cc$y = cc$y[id.filter];
		n = length(id.filter);
	}
	cc$x = c(cc$x, cc$x[1]);
	cc$y = c(cc$y, cc$y[1]);
	# Partition Line between Centers:
	partf = function(id) {
		x1 = cc$x[id]; x2 = cc$x[id+1];
		y1 = cc$y[id]; y2 = cc$y[id+1];
		tt = seq(0, p-1) / p;
		x  = (1-tt)*x1 + tt*x2;
		y  = (1-tt)*y1 + tt*y2;
		return(data.frame(x=x, y=y));
	}
	cc = lapply(seq(n), partf);
	cc = do.call(rbind, cc);
	return(cc);
}
# Hexagon Ring:
# - used only for Level 2;
#' @export
ring.hexa = function(n, R, r = 1, center = c(1,1),
		col = NULL, id.filter = NULL, phi = 0, phi.h = pi/6) {
	c2 = ngon(n, center = center, r = R, phi=phi);
	if(is.null(id.filter)) id.filter = seq(n);
	h2 = lapply(id.filter, function(id) {
		ngon(6, center = c(c2$x[id], c2$y[id]), r=r,
			col=col, phi = phi.h + phi);
	});
	invisible(h2);
}


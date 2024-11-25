#######################################
#
# Title: BioShapes
# Maintainer: L. Mada
#
# https://github.com/discoleo/BioShapes


### Grids


### Hexagonal Grid: Circular
# n = Number of levels;
# r = Radius of each hexagon;
# phi = Rotation of the grid/hexagons;
#' @export
grid.hexa = function(n, r = 1, center = c(1,1), col = NULL, phi = 0) {
	if(n < 0) stop("The number of levels n must be positive!");
	if(n == 0) return(as.bioshape.empty());
	if(is.null(col)) {
		# col = # TODO
	} else if(length(col) == 1) col = rep(col, n);
	### Grid:
	L = ngon(6, center=center, r=r, col = col[1], phi = phi + pi/6);
	G = list(L1 = L);
	if(n == 1) return(as.bioshape(G));
	# L2:
	L = ring.hexa(6, center=center, r=r, phi=phi,
		R = sqrt(3) * r, col = col[2]);
	G = c(G, L2 = as.bioshape(L));
	if(n == 2) return(as.bioshape(rev(G)));
	# L 2:n
	L = lapply(seq(2, n), function(id) {
		R = r * sqrt(3) * (id-1); nOuter = 6*id;
		id.filter = seq(1, nOuter, by = id);
		grid.hexa.ring(n = nOuter, p = id-1,
			id.filter=id.filter, R=R, r=r, center=center,
			col=col[id], phi=phi, phi.h = phi + pi/6);
	});
	names(L) = paste0("L", seq(2, n));
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


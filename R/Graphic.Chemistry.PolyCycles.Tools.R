#######################################
#
# Title: BioShapes
# Maintainer: L. Mada
#
# https://github.com/discoleo/BioShapes
#
# Continuation of:
# 1. Bachelor Thesis: Adrian Cotoc (2022-2023)
# 2. Bachelor Thesis: Darian Voda (2021-2022)



### Annotate Cycles
#' @export
annotate.polycycle = function(n, R = 4, center = c(0,0), lwd = 1, col = "red",
		as.letters = FALSE, phi = - 2*pi/n) {
	x = R * cos(2*pi*seq(n) / n + phi) + center[1];
	y = R * sin(2*pi*seq(n) / n + phi) + center[2];
	lbl = if(as.letters) LETTERS[seq(n)] else seq(n);
	lst = list(x = x, y = y, labels = lbl, lwd = lwd, col = col);
	class(lst) = c("text", "list");
	lst = list(lst);
	return(as.bioshape(lst));
}

#' @export
polycycle.polyEE = function(n = 14, ngon = 7, R = 4, dB = 0.4, dpi = 0.125, center = c(0,0)) {
	gg = polycycle.polyene(n=n, ngon=ngon, R=R, d = dpi, center=center);
	gg$EE = piBond.polycycle(gg$Polycycle, d = dB, dpi = dpi/2, center=center);
	invisible(gg);
}

#' @export
polycycle.polyHt = function(n = 14, ngon = 7, Ht = "O", col = "red", center = c(0,0),
		R = 4, dB = 0.4, dpi = 0.125, r.adj = 0, adj.Ht = NULL) {
	gg = polycycle.polyene(n=n, ngon=ngon, R=R, d = dpi, center=center);
	gg$EE = piBond.polycycle(gg$Polycycle, d = dB, dpi = dpi/2, center=center);
	# Ht-Labels:
	if(length(Ht) == 1) Ht = rep(Ht, n);
	if(length(col) == 1) col = rep(col, n);
	if(hasAdj <- (! is.null(adj.Ht)) && nrow(adj.Ht) == 1) {
		adj.Ht = matrix(rep(adj.Ht, each = n), ncol=2);
	}
	#
	rr = attr(gg$Polycycle, "r");
	xy = points.circle(n, r = R - rr - dB - 0.075 * R / rr + r.adj, center=center);
	Ht = list(x = xy$x, y = xy$y, labels = Ht, col = col);
	if(hasAdj) Ht$adj = adj.Ht;
	class(Ht) = c("text", "list");
	gg$Ht = Ht;
	invisible(gg);
}

#' @export
piBond.polycycle = function(x, d = 0.4, center = c(0, 0), dpi = 0.0625) {
	if(length(x) == 0) return(NULL);
	p1 = lapply(x, function(gg) c(gg$x[1], gg$y[1]));
	px = lapply(x, function(gg) c(center[1], gg$x[1]));
	py = lapply(x, function(gg) c(center[2], gg$y[1]));
	#
	qd = which.quadrant.polycycle(x);
	d  = ifelse(qd == 1 | qd == 4, -d, d);
	pO = lapply(seq(along = x), function(id) {
		shiftPoint(p1[[id]], px[[id]], py[[id]], d = d[[id]]);
	});
	pO = do.call(rbind, pO);
	p1 = do.call(rbind, p1); colnames(p1) = c("x1", "y1");
	pO = cbind(pO, p1);
	pD = lapply(seq(along = x), function(id) {
			shift.ortho(pO[id, c("x", "x1")], pO[id, c("y", "y1")], d = c(-dpi, dpi),
				id.offset = 2*(id - 1));
		});
	pD = do.call(rbind, pD);
	pD = list(pD);
	invisible(as.bioshape(pD));
}

#' @export
which.quadrant.polycycle = function(x) {
	sapply(x, function(gg) {
		n  = length(gg$x);
		n2 = (n %/% 2) + 1;
		if(n %% 2 == 0) {
			x = gg$x[c(1, n2)];
			y = gg$y[c(1, n2)];
		} else {
			x = c(gg$x[1], (gg$x[n2] + gg$x[n2 + 1]) / 2);
			y = c(gg$y[1], (gg$y[n2] + gg$y[n2 + 1]) / 2);
		}
		which.quadrant(x, y);
	});
}


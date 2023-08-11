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


### Molecular Biology


### Ig-Domain
#' @export
mol.IgDomain = function(x, y, t = 3/4, l = 1, d = 1.5 * l,
		col = 1, col.ig = col, lwd = 1, lwd.ig = lwd) {
	len = length(t);
	opt = list(l=l, d=d, col.ig=col.ig, lwd.ig=lwd.ig);
	if(len > 1) {
		# TODO: verify if "t" is properly sorted!
		if(any(diff(t) < 0)) warning("May not be properly sorted!");
		opt = rep.all.list(opt, len=len);
	}
	mid = cbind(1-t, t) %*% cbind(x, y);
	slope = slope(x, y);
	xy = sapply(seq(len), function(id) {
		l = opt$l[id];
		shiftPoint(mid[id, c(1,2)], slope=slope, d = c(-l, l));
	}); # t(x, x, y, y);
	# Ig-Loops:
	loops = lapply(seq(len), function(id) {
		xy = xy[, id];
		cc = circle.ArcByDist(xy[c(1,2)], xy[c(3,4)], d = opt$d[id],
			col = opt$col[id], lwd = opt$lwd[id]);
		return(cc[[1]]);
	})
	loops = as.bioshape(loops);
	# Base-Line:
	xy = xy[c(1,3,2,4), , drop=FALSE]; dim(xy) = c(2, ncol(xy) * 2);
	xy = cbind(c(x[1], y[1]), xy, c(x[2], y[2]));
	xy = split.AltLines.matrix(t(xy));
	#
	lst = list(Lines = xy, Loops = loops, lwd=lwd, col=col);
	lst = as.bioshape(lst);
	invisible(lst);
}


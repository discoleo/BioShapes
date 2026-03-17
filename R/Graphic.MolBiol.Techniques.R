###################
#
# Title: BioShapes
#
# Maintainer: L. Mada
#
#
# GitHub: https://github.com/discoleo/BioShapes


### Molecular Biology: Techniques

### Gel-Electrophoresis

### Western Blot

#' @export
gel = function(x, y, n = 8, h = 8, h.well = 1.5, d = 1/4,
		col = c("#B2B2B2"), fill = c("#C2C2C2", "#3224E8"),
		tray.fill = 0.75) {
	# Gel Boundaries:
	xy = as.rect.fromLine(x, y, h=h);
	if(length(n) == 0) {
		class(xy) = c("polygon", class(xy));
		xy = list(Gel = xy);
		return(as.bioshape(xy));
	}
	nn = n; n = n[1];
	# Bottom of wells:
	tb = h.well / dist.xy(x, y); tbi = 1 - tb;
	x1 = tb * x[1] + tbi * xy$x[4];
	y1 = tb * y[1] + tbi * xy$y[4];
	x2 = tb * x[2] + tbi * xy$x[3];
	y2 = tb * y[2] + tbi * xy$y[3];
	xB = c(x1, x2); yB = c(y1, y2);
	# Lanes:
	xyT = split.withSep(xy$x[c(4,3)], xy$y[c(4,3)], n=n, d=d);
	xyB = split.withSep(xB, yB, n=n, d=d);
	# Merge:
	n1 = n + 2; nE = 2*n1 + 1; nE = c(nE, nE+1);
	id = seq(1, by = 2, length.out = n1);
	xx = rbind(xyT$x[id], xyT$x[id+1], xyB$x[id+1], xyB$x[id+2]);
	yy = rbind(xyT$y[id], xyT$y[id+1], xyB$y[id+1], xyB$y[id+2]);
	xx = c(x[1], as.vector(xx), xyT$x[nE], x[2]);
	yy = c(y[1], as.vector(yy), xyT$y[nE], y[2]);
	res = list(x = xx, y = yy, col = col[1], fill = fill[1]);
	class(res) = c("polygon", "list");
	res = list(Gel = res);
	if(length(nn) == 1) return(as.bioshape(res));
	### Loaded Gel:
	tf = tray.fill; tf1 = 1 - tf;
	xyTF   = list();
	xyTF$x = tf * xyT$x + tf1 * xyB$x;
	xyTF$y = tf * xyT$y + tf1 * xyB$y;
	trays = lapply(seq(min(nn)), function(id) {
		id2 = 2*id;
		xi = c(xyTF$x[c(id2, id2 + 1)], xyB$x[c(id2 + 1, id2)]);
		yi = c(xyTF$y[c(id2, id2 + 1)], xyB$y[c(id2 + 1, id2)]);
		tmp = list(x = xi, y = yi, col = fill[2], fill = fill[2]);
		class(tmp) = c("polygon", "list");
		return(tmp);
	})
	trays = as.bioshape(trays);
	res$Trays = trays;
	return(as.bioshape(res));
}


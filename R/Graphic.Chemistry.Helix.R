#######################################
#
# Title: BioShapes
# Maintainer: L. Mada
#
# https://github.com/discoleo/BioShapes


### Biochemistry: Helix Structures


#' @export
seq.helix.numeric = function(n, by = 11, helix.length = 18) {
	r0 = n %% helix.length;
	nt = n;
	if(r0 != 0) nt = n + helix.length - r0;
	id = seq(0, by = by, length.out = nt);
	id = (id %% helix.length) + 1;
	if(nt > helix.length) {
		n0 = (nt - 1) %/% helix.length ;
		id.off = seq(0, by = helix.length, length.out = n0 + 1);
		id.off = rep(id.off, each = helix.length);
		id = id + id.off;
	}
	if(r0 != 0) id[id > n] = NA;
	return(id);
}

#' @export
helix.wheel = function(x, r = 2, center = c(0,0), phi = 0,
		col = 1, fill = "yellow", col.line = "#969696",
		clock = TRUE, by = 11, helix.length = 18) {
	id  = seq.helix.numeric(length(x), by=by, helix.length=helix.length);
	len = length(id);
	# Layers: full turns
	nLayers = len / helix.length;
	R1 = r;
	r1 = R1 * sin(pi / helix.length);
	rr = r1 * 2;
	Rn = R1 + rr * seq(0, nLayers);
	x  = x[id];
	lst = lapply(seq(nLayers), function(nL) {
		id0 = nL * helix.length - c(17,0);
		lbl = x[seq(id0[1], id0[2])];
		xyC = points.circle(helix.length, r = Rn[nL], center=center, phi=phi, clock=clock);
		# NA:
		isNA = is.na(lbl);
		if(any(isNA)) {
			isNA = ! isNA;
			lbl = lbl[isNA]; xyC = list(x = xyC$x[isNA], y = xyC$y[isNA]);
		}
		# Text
		lbl = list(x = xyC$x, y = xyC$y, labels = lbl, col=col);
		class(lbl) = c("text", "list");
		# Circles
		xyC = list(center = cbind(xyC$x, xyC$y), r = r1, col = NA, fill=fill);
		xyC = as.circle(xyC);
		lst = list(C = xyC, T = lbl);
		lst = as.bioshape(lst);
		return(lst);
	});
	# Lines
	if( ! is.null(col.line)) {
		# Simple Line: same entry/exit;
		id  = seq.helix.numeric(helix.length, by=by, helix.length=helix.length);
		tmp = points.circle(helix.length, r = Rn[1] - r1, center=center, phi=phi, clock=clock);
		id  = order(id);
		tmp = list(x = tmp$x[id], y = tmp$y[id], col = col.line);
		lst$L = tmp;
	}
	return(as.bioshape(lst));
}




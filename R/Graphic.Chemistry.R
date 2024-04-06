#######################################
#
# Title: BioShapes
# Maintainer: L. Mada
#
# https://github.com/discoleo/BioShapes


### Draw 2D molecules

# Basic algorithm

#' @export
molecule.phi = function(phi, xy = c(0,0), r = 1, phi0 = 0, th = NULL) {
	xy = molecule.coords(phi=phi, xy=xy, r=r, phi0=phi0,
		th=th, verbose = FALSE);
	xy = list(list(x = xy[,1], y = xy[,2]));
	class(xy) = c("lines.list", class(xy));
	return(as.bioshape(list(Mol = xy)));
}
#' @export
molecule.coords = function(phi, xy = c(0,0), r = 1, phi0 = 0,
		th = NULL, verbose = FALSE) {
	if(is.null(th)) {
		phi[1] = phi[1] + phi0;
		th = expand.angle(phi);
	}
	# TODO: phi0 != 0 && missing(phi);
	# xy:
	n = length(th);
	x = rep(0, n); x[1] = xy[1];
	y = rep(0, n); y[1] = xy[2];
	for(i in seq(2, n)) {
		thi = th[i-1];
		x[i] = cos(thi)*r + x[i-1];
		y[i] = sin(thi)*r + y[i-1];
	}
	x[abs(x) < 1E-8] = 0;
	y[abs(y) < 1E-8] = 0;
	xy = data.frame(x=x, y=y);
	if(verbose) {
		# in Degrees:
		# pd = round(phi * 180 / pi, 2);
		xy$th = th;
	}
	return(xy);
}

#' @export
expand.angle = function(x) {
	n = length(x);
	th = x;
	pi2 = 2*pi;
	as2pi = function(x) x - pi2 * round(x / pi2);
	th[1] = as2pi(th[1]);
	for(i in seq(2, n)) {
		sg = if(th[i-1] <= pi) 1 else -1;
		th[i] = as2pi(sg*(th[i-1] - pi) + th[i]);
	}
	return(th);
}

##################
### Transforms ###

### Pi-Bonds
#' @export
as.piBond = function(xy, which, d = 1/8) {
	len = length(xy);
	if(len == 0) return(xy);
	if(length(d) == 1) d = rep(d, len);
	ln = lapply(seq(len), function(id) {
		xy = xy[[id]];
		ln = shift.ortho(xy$x[which], xy$y[which], d = d[id]);
		ln = list(x = ln$x, y = ln$y);
		return(ln)
	});
	class(ln) = c("lines.list", "list");
	invisible(ln);
}

### Simple Ligands
#' @export
ligandArrow = function(x, y, slope=Inf, solid=TRUE, d = 0.75, w = 0.125 * d) {
  pxy = if(missing(y)) x else c(x, y);
  # if(slope < 0) d = -d;
  p  = shiftPoint(pxy, slope=slope, d=d);
  pV = shiftLine(p, slope=slope, d = c(-w, w));
  list(x = c(pxy[1], pV$x), y = c(pxy[2], pV$y));
}

#' @export
ligand = function(x, cyc, i, slope=Inf, solid = TRUE, col = NULL, col.grey = "grey50", ...) {
  ligandBase = function(id) {
    cyc = if(length(cyc) == 1) cyc else cyc[[id]];
    tmp = x[[cyc]];
    pol = ligandArrow(tmp$x[i[[id]]], tmp$y[i[[id]]],
                      slope=slope[[id]], solid=solid[[id]], ...);
	if(is.null(col)) {
		col  = if(solid[[id]]) 1 else col.grey;
		fill = col;
    } else {
		col  = if(length(col) == 1) col else col[[id]];
		fill = NA;
	}
    lst = c(pol, col=col, fill=fill);
    class(lst) = c("polygon", class(lst));
    return(lst);
  }
  len = length(i);
  pol = lapply(seq(len), ligandBase);
  class(pol) = c("bioshape", class(pol))
  return(pol);
}


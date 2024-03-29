#######################################
#
# Title: BioShapes
# Maintainer: L. Mada
#
# https://github.com/discoleo/BioShapes
#
# Continuation of:
# 1. Bachelor Thesis: Adrian Cotoc (2022-2023)
# Faculty of Mathematics and Informatics, UVT
# Coordinator:
#   Prof. Daniela Zaharie
#   Dr. med. Leonard Mada (Syonic SRL)
#   in collaboration with Syonic SRL
#   GitHub: https://github.com/Adi131313/BioShapes
#
# 2. Bachelor Thesis: Darian Voda (2021-2022)


### Poly-cyclic chemical molecules

### Parser
# x = String encoding structure;
#   - numbers = dimension of cycle;
#   - separators = type of junction between cycles;
#' @export
parseCycles = function(x, r=1, d2=0.0625) {
  reg = "(?<=[0-9])(?=[^0-9])|(?<=[^0-9])(?=[0-9])";
  cyc = strsplit(x, reg, perl=TRUE);
  cyc = cyc[[1]];
  #
  l = list();
  x0 = 0; y0 = 0; r0 = r;
  r  = 0; d  = 0;
  phi0 = 0; phi = 0; s = "";
  isNum = grepl("^[0-9]", cyc);
  nL = rep(NA, length(cyc));
  nL[isNum] = as.numeric(cyc[isNum]);
  for(id in seq(length(cyc))) {
    n = nL[id];
    if( ! is.na(n)) {
      halpha = pi/n;
      r = r0 / 2 / sin(halpha);
      d = r * cos(halpha);
      # Correction to rotation:
      dp = if(phi0 == 0) d else (d + r + d / cos(phi + phi0));
      phi = if(n %% 2 == 1) 0 else halpha;
      lPp = pointsCircle(n, r=r, center=c(x0 + dp, y0), phi = phi + phi0);
      lPp$x = c(lPp$x, lPp$x[1]);
      lPp$y = c(lPp$y, lPp$y[1]);
      # Reset phi0
      if(phi0 != 0 && round(phi0 - pi, 8) == 0) {
        phi0 = 0;
      }
      l = c(l, list(lPp));
      next;
    }
    s = cyc[id];
    if(s == "|") {
      # TODO: "5<5|4"
      x0 = x0 + d + (if(nL[id - 1] %% 2 == 1) r else d);
      next;
    }
    if(s == "<") {
      x0 = x0 + d + r;
      phi0 = phi0 + pi;
      next;
    }
    if(s == "-" || s == "=") {
      # TODO: if(phi != 0)
      x0 = x0 + d + r;
      x1 = x0 + r0;
      l1 = list(x = c(x0, x1), y = c(y0, y0));
      if(s == "=") {
        l1 = shiftLine(l1$x, l1$y, d = c(-d2, d2));
		l1 = as.bioshape(list(DB = l1));
        # l1 = split(l2[, 1:2], l2$id);
      } else l1 = list(l1);
      l = c(l, l1);
      x0 = x1;
      phi0 = phi0 + pi;
      next;
    }
    if(s == "\\") {
      # TODO: now only for n = 6;
      x0 = x0 + d;
      y0 = y0 + r + r0/2;
      next;
    }
  }
  class(l) = c("chemistry", class(l));
  return(l)
}


### Polycyclic Cycles
#' @export
polycycle.cyc = function(n = 14, ngon = 7, R = 4, center = c(0,0)) {
	if(ngon < 3) stop("Invalid cycle!");
	# ngon = 4: does not work either;
	gg = circle.spiro(n=n, ngon=ngon, R=R, center=center, type = "in");
	if(ngon == 3) return(gg);
	bb = c(ngon, ngon - 1);
	gg = as.mean.xy(gg, bb, c(2,3));
	return(invisible(gg));
}


### Polyene
#' @export
polycycle.polyene = function(n = 14, ngon = 7, R = 4, d = 0.125, center = c(0,0)) {
	gg = polycycle.cyc(n=n, ngon=ngon, R=R, center=center);
	gg = list(Polycycle = gg);
	### Double Bonds
	whichq = function(id) {
		sapply(gg$Polycycle, function(gg) {
			which.quadrant(gg$x[id], gg$y[id]);
		});
	}
	id = c(2,3);
	qd = whichq(id);
	d2 = ifelse(qd == 1 | qd == 2, d, - d);
	ln = as.piBond(gg$Polycycle, id, d = d2);
	gg$Pi = as.bioshape(ln);
	# Out:
	if(ngon >= 7 && ngon %% 2 == 1) {
		ii = seq(4, ngon - 3, by = 2);
		n2 = n / 2; n4 = n %/% 4;
		ln = lapply(ii, function(i) {
			id = c(i, i+1);
			qd = whichq(id);
			d2 = ifelse(qd == 1 | qd == 2, d, - d);
			ln = as.piBond(gg$Polycycle, id, d = d2);
		})
		gg$PiOut = as.bioshape(ln);
	}
	invisible(as.bioshape(gg));
}

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

##################

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


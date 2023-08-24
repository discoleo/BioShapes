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


### Ig-like Domains
# t = positions of Ig-domains;
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
	qdr = which.quadrant(x, y);
	slope = slope(x, y);
	xy = sapply(seq(len), function(id) {
		l = opt$l[id];
		if(qdr == 2) l = - l;
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


### Ig Monomer

# TODO:
# - [started] fix all orientations;
# - Styles;


### Ig Heavy Chains: Basic Backbone
#' @export
mol.IgH = function(xy, height = 6, t.Hinge = 2/5, d.HH = 1/2, theta = pi/3, phi = pi/2) {
	th2 = theta / 2;
	slope = tan(phi);
	isUp  = ((phi >= 0) && (phi <= pi/2)) ||
			(phi >= 3*pi/2);
	dsgn  = if(isUp) height else - height;
	isRUp = (phi >= th2) && (phi - th2 <= pi);
	isLUp = (phi + th2 >= 0) && (phi + th2 <= pi);
	# Bottom & Top: Start & End
	pS = shift.ortho(xy, slope=slope, d = c(-d.HH, d.HH)/2);
	pE = shift.points.df(pS, slope=slope, d = dsgn, simplify = FALSE);
	### Variable Region (Hinge):
	# Start:
	as.m = function(id) as.matrix(rbind(pS[id, c("x", "y")], pE[[id]]))
	mid1 = c(t.Hinge, 1 - t.Hinge) %*% as.m(1);
	mid2 = c(t.Hinge, 1 - t.Hinge) %*% as.m(2);
	d2 = height * (1 - t.Hinge);
	dd = d2 * tan(theta/2); # d(Tip VRegion, Midline)
	# Top:
	pT1 = shift.ortho(pE[[1]], slope=slope, d = - dd);
	pT2 = shift.ortho(pE[[2]], slope=slope, d = dd);
	dHR = if(isRUp) - d.HH else   d.HH;
	dHL = if(isLUp)   d.HH else - d.HH;
	LL1 = shift.ortho.df(rbind(mid1, pT1[, c("x", "y")]), d = dHR);
	LL2 = shift.ortho.df(rbind(mid2, pT2[, c("x", "y")]), d = dHL);
	lst = list(pS=pS, pE=pE, mid1=mid1, mid2=mid2,
		pT1=pT1, pT2=pT2, LL1=LL1, LL2=LL2, qq = c(isRUp, isLUp));
	return(lst);
}

# xy  = base of Ig-molecule;
# phi = slope of molecule;
# theta = angle between the 2 chains;
#' @export
mol.Ig = function(xy, height = 6, t.Hinge = 2/5, d.HH = 1/2, d.rel = 1/4,
		n = c(2,2), theta = pi/3, phi = pi/2) {
	phi = as.radians0(phi);
	th2 = theta / 2;
	Ig = mol.IgH(xy=xy, height=height,  t.Hinge=t.Hinge, d.HH=d.HH,
		theta=theta, phi=phi);
	mid1 = Ig$mid1; mid2 = Ig$mid2;
	d2 = height * (1 - t.Hinge);
	lL = d2 / cos(th2); # length of Light Chain
	# Ig-Domains:
	if(length(n) == 1) n = c(n, 0);
	lst = list();
	if(n[1] > 0) {
		nS  = n[1] + 1; l = lL / nS * d.rel;
		tL  = seq(n[1]) / nS;
		LC1 = mol.IgDomain(Ig$LL1$x, Ig$LL1$y, t = tL, l = l, d = - 1.5 * l);
		LC2 = mol.IgDomain(Ig$LL2$x, Ig$LL2$y, t = tL, l = l, d = 1.5 * l);
		LC  = as.bioshape(list(LC1=LC1, LC2=LC2));
		# HV:
		HV1x = c(mid1[1], Ig$pT1$x[1]); HV1y = c(mid1[2], Ig$pT1$y[1]);
		HV2x = c(mid2[1], Ig$pT2$x[1]); HV2y = c(mid2[2], Ig$pT2$y[1]);
		HV1 = mol.IgDomain(HV1x, HV1y, t = tL, l = l, d = 1.5 * l);
		HV2 = mol.IgDomain(HV2x, HV2y, t = tL, l = l, d = - 1.5 * l);
		HV  = as.bioshape(list(HV1=HV1, HV2=HV2));
		lst = c(lst, HV = HV, LC = LC);
	}
	# Base: HC
	HC1x = c(Ig$pS$x[1], mid1[1]); HC1y = c(Ig$pS$y[1], mid1[2]);
	HC2x = c(Ig$pS$x[2], mid2[1]); HC2y = c(Ig$pS$y[2], mid2[2]);
	if(n[2] > 0) {
		nS  = n[2] + 1; l = height * (1 - t.Hinge) / nS * d.rel;
		tL  = seq(n[2]) / nS;
		HC1 = mol.IgDomain(HC1x, HC1y, t = tL, l = l, d = - 1.5 * l);
		HC2 = mol.IgDomain(HC2x, HC2y, t = tL, l = l, d = 1.5 * l);
		HC  = as.bioshape(list(HC1=HC1, HC2=HC2));
	} else {
		HC = as.bioshape(list(
			HC1 = list(x = HC1x, y = HC1y),
			HC2 = list(x = HC2x, y = HC2y)) );
	}
	lst = c(lst, Heavy = HC);
	#
	return(as.bioshape(lst));
}


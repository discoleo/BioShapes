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

### Golgi

# 1 Tubule:
#' export
golgi.p1 = function(r = 1, d = 0.25, xy = c(0,0), phi = 0, theta = pi/3) {
	r = c(r, r + d);
	phi = phi + c(0, theta);
	C1.x = r*cos(phi[1]) + xy[1];
	C1.y = r*sin(phi[1]) + xy[2];
	C2.x = r*cos(phi[2]) + xy[1];
	C2.y = r*sin(phi[2]) + xy[2];
	top = c(FALSE, TRUE); # TODO: proper;
	cap1 = circle.ArcByDiam(C1.x, C1.y, top = top[1]);
	cap2 = circle.ArcByDiam(C2.x, C2.y, top = top[2]);
	lst = list(C1 = as.circle.arc(list(r = r[1], center = xy, phi = phi)));
	lst$C2 = as.circle.arc(list(r = r[2], center = xy, phi = phi));
	lst$Cap1 = cap1;
	lst$cap2 = cap2;
	invisible(as.bioshape(lst));
}

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
		if(qdr == 2 || qdr == 3) l = - l;
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


### Ig: Basic Backbone
#' @export
mol.IgBB = function(xy, height = 6, t.Hinge = 2/5, t.LC = c(0, 1), d.HH = 1/2,
		d.LH = d.HH, theta = pi/3, phi = pi/2, debug = FALSE) {
	th2 = theta / 2;
	slope = tan(phi);
	qq = which.quadrant.phi(phi);
	dsgn = if(qq == 1 || qq == 4) height else - height;
	# Bottom & Top: Start & End
	pS = shift.ortho(xy, slope=slope, d = c(-d.HH, d.HH)/2);
	pE = shift.points.df(pS, slope=slope, d = dsgn, simplify = FALSE);
	### Variable Region (Hinge):
	# Start:
	as.m = function(id) as.matrix(rbind(pS[id, c("x", "y")], pE[[id]]))
	midL = c(t.Hinge, 1 - t.Hinge) %*% as.m(1);
	midR = c(t.Hinge, 1 - t.Hinge) %*% as.m(2);
	d2 = height * (1 - t.Hinge);
	dd = d2 * tan(theta/2); # d(Tip VRegion, Midline)
	# Top:
	pTL = shift.ortho(pE[[1]], slope=slope, d = - dd);
	pTR = shift.ortho(pE[[2]], slope=slope, d = dd);
	# LC = Light Chains
	qR = which.quadrant.phi(phi - th2);
	qL = which.quadrant.phi(phi + th2);
	if(debug) cat("Quadrant: ", c(qR, qL), "\n");
	sgnR = if(qR == 2) -1 else  1;
	sgnL = if(qL == 1)  1 else -1;
	dHR = sgnR * d.LH;
	dHL = sgnL * d.LH;
	LCL = shift.ortho.df(rbind(midL, pTL[, c("x", "y")]), d = dHL);
	LCR = shift.ortho.df(rbind(midR, pTR[, c("x", "y")]), d = dHR);
	#
	tt = rbind(c(1 - t.LC[1], t.LC[1]), c(1 - t.LC[2], t.LC[2]));
	shiftL = function(xy) {
		xy$x = tt %*% xy$x;
		xy$y = tt %*% xy$y;
		return(xy);
	}
	LCL = shiftL(LCL); LCR = shiftL(LCR);
	#
	lst = list(pS=pS, pE=pE, mid1=midL, mid2=midR,
		pT1=pTL, pT2=pTR, LL1=LCL, LL2=LCR, qq = list(qq1=qL, qq2=qR));
	return(lst);
}

# n = number of Ig-domains in the V-region & HC-region;
# xy  = base of Ig-molecule;
# phi = slope of molecule;
# theta = angle between the 2 chains;
# d.HH  = distance between the 2 Heavy Chains;
# d.LH  = distance between the H & L Chains;
# d.rel = relative length of the Ig-Domains;
#' @export
mol.Ig = function(xy, height = 6, t.Hinge = 2/5, d.HH = 1/2, d.LH = d.HH, d.rel = 1/4,
		n = c(2,2), t.LC = c(0, 1), phi = pi/2, theta = pi/3) {
	phi = as.radians0(phi);
	th2 = theta / 2;
	Ig = mol.IgBB(xy=xy, height=height, t.Hinge=t.Hinge, t.LC=t.LC, d.HH=d.HH, d.LH=d.LH,
		theta=theta, phi=phi);
	mid1 = Ig$mid1; mid2 = Ig$mid2;
	d2   = height * (1 - t.Hinge);
	lenL = d2 / cos(th2); # length of Light Chain
	# Ig-Domains:
	if(length(n) == 1) n = c(n, 0);
	lst = list();
	if(n[1] == 0) {
		HV1x = c(Ig$pT1$x, mid1[1]); HV1y = c(Ig$pT1$y, mid1[2]);
		HV2x = c(Ig$pT2$x, mid2[1]); HV2y = c(Ig$pT2$y, mid2[2]);
		VV = as.bioshape(list(
			LV1 = Ig$LL1, LV2 = Ig$LL2,
			HV1 = list(x = HV1x, y = HV1y),
			HV2 = list(x = HV2x, y = HV2y)) );
		lst = c(lst, VV = list(VV));
	} else if(n[1] > 0) {
		nS  = n[1] + 1; l = lenL / nS * d.rel;
		tL  = seq(n[1]) / nS;
		qq  = Ig$qq;
		sg1 = if(qq$qq1 == 1) -1 else 1;
		sg2 = if(qq$qq2 == 2) -1 else 1;
		# d.IgDomain.rel = 1.5; hardcoded!
		LC1 = mol.IgDomain(Ig$LL1$x, Ig$LL1$y, t = tL, l = l, d = - sg1 * 1.5 * l);
		LC2 = mol.IgDomain(Ig$LL2$x, Ig$LL2$y, t = tL, l = l, d = sg2 * 1.5 * l);
		LC  = as.bioshape(list(LC1=LC1, LC2=LC2));
		# HV:
		HV1x = c(mid1[1], Ig$pT1$x[1]); HV1y = c(mid1[2], Ig$pT1$y[1]);
		HV2x = c(mid2[1], Ig$pT2$x[1]); HV2y = c(mid2[2], Ig$pT2$y[1]);
		HV1 = mol.IgDomain(HV1x, HV1y, t = tL, l = l, d = sg1 * 1.5 * l);
		HV2 = mol.IgDomain(HV2x, HV2y, t = tL, l = l, d = - sg2 * 1.5 * l);
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

###########

###########
### TLR ###

### Toll-Like Receptors

#' @export
mol.tlr = function(x, y, r = 2, n = 28, phi = - pi/6,
		lwd = c(1,4), col = "#328012", fill = "#64C232",
		type = c("Dimer", "Monomer"),
		d.rel = - 0.25, r.ll.rel = c(0.375, 0.5)) {
	type = match.arg(type);
	dphi = atan2(y[2] - y[1], x[2] - x[1]) - pi/2;
	molBase = function(phi) {
		# Leucine-Rich Domain:
		phi.cc = phi[1];
		cc = c(
			x[2] - r*cos(phi.cc),
			y[2] - r*sin(phi.cc));
		ll = mol.tlr.lld(cc, r=r, n=n, phi=phi, r.ll.rel=r.ll.rel,
			lwd = lwd[[1]], col = col[[1]], fill = fill[[1]]);
		# Trans-Membrane:
		# TODO: better concept;
		mb = list(x = x, y = y, lwd = lwd[[2]], col = col[[1]]);
		# Full:
		lst = list(LL = ll, Mbr = mb);
	}
	lst = as.bioshape(molBase(phi + dphi));
	if(type == "Dimer") {
		if(dphi >= - pi/2) d.rel = - d.rel;
		d  = dist.xy(x, y) * d.rel;
		xy = shift.ortho(x, y, d=d);
		x  = xy[,1]; y = xy[,2];
		lst = list(TLR1 = lst);
		phi = dphi - phi;
		phi = c(pi + phi, phi);
		if(length(col)  > 1) col  = col[[2]];
		if(length(fill) > 1) fill = fill[[2]];
		lst$TLR2 = as.bioshape(molBase(phi));
	}
	return(as.bioshape(lst));
}

# Leucin-rich Domain
#' @export
mol.tlr.lld = function(xy, r = 2, n = 28, phi = pi/3, r.ll.rel = c(0.375, 0.5),
		lwd = 1, col = "#328012", fill = "#64C232") {
	if(length(phi) == 1) phi = c(phi, phi + pi);
	dp = diff(phi);
	rL = r * r.ll.rel[1]; # Long axis
	rS = r * r.ll.rel[2] * 2 * sin(dp / (2*n)); # Short axis
	rE = c(rS, rL);
	th = seq(0, n, length.out = n) * dp / n;
	th = th + phi[1];
	lst = lapply(seq(n), function(id) {
		th = th[id];
		x  = r * cos(th) + xy[1];
		y  = r * sin(th) + xy[2];
		lst = list(r = rE, center = c(x, y), phi = c(0, 2*pi), th = th + pi/2,
			lwd=lwd, col=col, fill=fill);
		return(as.ellipse(lst));
	});
	return(as.bioshape(lst));
}


####################
####################

# e.g. Transmembrane receptors;
#' @export
protein.domains = function(xy, n, w = 1, h = 4*w, d = 0.25, slope = c(Inf, 0),
		lwd.con = 2, ...) {
	lst = cylinder.tubes(xy=xy, n=n, w=w, h=h, d=d, slope=slope, ...);
	applyCon = function(i0, nm, top) {
		lapply(seq(i0, n - 1, by=2), function(id) {
			xy1 = lst[[id]][[nm]]$center;
			xy2 = lst[[id + 1]][[nm]]$center;
			lst = circle.ArcByDiam.xy(xy1, xy2, top=top);
			lst$lwd = lwd.con;
			return(lst);
		});
	}
	conT = applyCon(1, "Cap2", top = TRUE);
	conB = applyCon(2, "Cap1", top = FALSE);
	lst = c(lst, conT, conB);
	invisible(as.bioshape(lst));
}


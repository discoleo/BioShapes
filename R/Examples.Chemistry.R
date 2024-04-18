


### Cholesterol

#' @export
cholesterol = function(origin = c(0,0), dPi = 0.125) {
	lst = cholesterol.bb(origin = origin);
	lst$Me1 = ligand(lst, 1, 1, slope = Inf);
	lst$Me2 = ligand(lst, 3, 1, slope = Inf);
	lst$Pi  = as.piBond(lst[2], which = c(4,5), d = - dPi)[[1]];
	# OH
	lstOH = ligand(lst, 1, i = 4, slope = 1/sqrt(3), d = -0.75)[[1]];
	txtOH = list(x = lstOH$x[2], y = lstOH$y[2], labels = "HO", adj = c(1.2, 0.8));
	class(txtOH) = c("text", "list");
	lstOH = list(OHb = lstOH, OHl = txtOH);
	class(lstOH) = c("bioshape", "list");
	lst$OH = lstOH;
	# TODO: C17-side chain;
	invisible(lst);
}
#' @export
cholesterol.bb = function(origin = c(0,0)) {
	lst = parseCycles("6|6\\6|5");
	if(any(origin != 0)) {
		lst = shift.listNR(lst, by = origin);
	}
	return(lst);
}

#' @export
example.cholesterol = function(origin = c(0,0)) {
	par.old = par(mfrow = c(2,1))
	
	### Cholesterol backbone
	plot.base(ylim = c(-2,4) + origin[2]);
	lines(cholesterol.bb(origin = origin));
	
	### Cholesterol
	plot.base(ylim = c(-2,4) + origin[2]);
	lines(cholesterol(origin = origin));
	
	par(par.old);
}


######################

### Polycycles

#' @export
test.polycycle = function(n = c(14, 17, 11, 19), ngon = 7,
		Ht = c("O", "S", NA, ""), col = c("red", "#FFE832", "grey", "grey"), R = 4) {
	if(length(R) == 1) R = rep(R, 4);
	if(length(ngon) == 1) ngon = rep(ngon, 4);
	#
	par.old = par(mfrow = c(2,2));
	for(i in seq(4)) {
		lim = R[i] + 1; lim = c(- lim, lim);
		plot.base(xlim = lim, ylim = lim, axt = NULL, asp = 1);
		if(is.na(Ht[i])) {
			lines(polycycle.polyene(n[i], ngon[i]));
		} else if(nchar(Ht[i]) == 0) {
			lines(polycycle.polyEE(n[i], ngon[i]));
		} else {
			lines(polycycle.polyHt(n[i], ngon[i], Ht = Ht[i], col = col[i]));
		}
		lines(annotate.polycycle(n[i]))
	}
	par(par.old);
	invisible();
}

##################

### Vitamin A
# TODO: misses the double bonds;
#' @export
vitamin.A = function(x = 0, y = 0, r = 1, phi0 = 0, adj = NULL, rev.OH = FALSE) {
	# TODO: parse without overlaps;
	th = c(5/4, -1,1,1,1,
		-1+1/4,3,-1/2,3, 2 + 1/4,1,3, # END Cyclo 6;
		1,1,-1, 2,3,2,-1, 1,-1,2,3,2,-1, 1,-1);
	xy = molecule.coords(2*pi/3 * th, xy = c(x, y), r=r, phi0=phi0);
	# Proper Structure:
	idMe = list(c(7,6), c(9,6), c(1,2), c(17,16), c(23,22), c(27,26));
	Me = lapply(idMe, function(id) {
		list(x = xy[id, 1],
			y = xy[id, 2]);
	});
	# BackBone:
	idBB = c(6:2, 11, 14:16, 19:22, 25:26);
	BB = list(
		x = xy[idBB, 1],
		y = xy[idBB, 2]);
	xy = as.bioshape(list(list(BB=BB), SC = as.bioshape(Me)));
	xy$SC$Cyclo = list(
		x = xy[[1]]$BB$x[c(1,6)],
		y = xy[[1]]$BB$y[c(1,6)]);
	txt = if(rev.OH) "HO" else "OH";
	OH = list(x = xy$SC[[6]]$x[1], y = xy$SC[[6]]$y[1], labels = txt);
	if( ! is.null(adj)) OH$adj = adj;
	class(OH) = c("text", class(OH));
	xy$SC$OH = OH;
	return(xy);
}

#' @export
test.vitamin.A = function() {
	layout(matrix(c(1,2,1,3), 2));
	par.old = par(mar = c(5, 4, 1, 2) + 0.1);
	#
	xy = vitamin.A();
	plot.molecule(xy, xlim = c(-3, 9), main = NULL);
	#
	xy = vitamin.A(phi0 = pi * 1/5, adj = c(0, 0.5));
	plot.molecule(xy, xlim = c(-4, 8), cex = 0.9);
	xy = vitamin.A(phi0 = pi * 6/5, adj = c(1, 0.5), rev = TRUE);
	plot.molecule(xy, xlim = c(-8, 4), cex = 0.9);
	#
	par(par.old);
	invisible();
}


###########################

### Test: Newman Projection

#' @export
test.proj.newman = function(phi.add = c(0,0), ligand = "F,H,H|OH,\\CH[3],H") {
	plot.base()
	# Row 1:
	lines(proj.newman(ligand, center = c(1,6), phi = c(0.7, 1) + phi.add))
	lines(proj.newman(ligand, center = c(6,6), phi = c(1.6, 1.3) + phi.add))
	# Row 2:
	lines(proj.newman(ligand, center = c(1,1), phi = c(0.5, 0.8) + phi.add))
	lines(proj.newman(ligand, center = c(6,1) + phi.add))
}


###########################

### Schiffer–Edmundson Helical Wheel Diagrams
# Barroso, C, et al. (2020)
# The Diverse Piscidin Repertoire of the European Sea Bass (Dicentrarchus labrax):
# Molecular Characterization and Antimicrobial Activities.
# https://doi.org/10.3390/ijms21134613
#
#' @export
example.helix.piscidin = function(r = 2, as.position = FALSE, add.labels = TRUE,
		fill = "yellow", col.lines = "#969696", col.arrow = "red",
		col.labels = c("#D06420", "blue", "red"),
		lwd.arrow = 2, dy = NULL, cex = NULL, cex.title = 1.5) {
	par.old = par(mfrow = c(1,2));
	center = c(4,4);
	lim = c(0, 8);
	if(length(r) == 1) r = c(r,r);
	
	lines.ellipseDomain = function(id, data, labels = NULL,
			phi = c(0, pi), r.scale = 3.25, col = 1, lwd = 2) {
		xy = data$C$center[id, ];
		cx = mean(xy[, 1]);
		cy = mean(xy[, 2]);
		rr = data$C$r[[1]];
		th = atan2(diff(xy[,2]), diff(xy[,1]));
		d  = dist.xy(xy[,1], xy[,2]) / 2;
		lines(ellipse(cx, cy, r = c(d + rr, rr*r.scale),
			theta = th, phi=phi, col=col, lwd=lwd));
		if( ! is.null(labels)) {
			thpi = th / pi;
			sg = if(thpi <= -3/4) - 1.25 else 1.25; # print(th / pi * 180);
			if(is.null(cex)) cex = 1;
			text(cx, cy + 1.25 * sg*rr*r.scale, labels=labels, col=col, cex=cex);
		}
	}
	
	# Piscidin 1
	# FFHHIFRGIVHVGKTIHRLVTG
	x = strsplit("FFHHIFRGIVHVGKTIHRLVTG", "")[[1]];
	if(as.position) x = seq(length(x));
	h1 = helix.wheel(x, r = r[1], center=center, cex=cex,
		fill=fill, col.lines=col.lines, lwd.arrow=lwd.arrow, col.arrow=col.arrow);
	# Plot:
	plot.base(xlim = lim, ylim = lim, axt = NULL);
	lines(h1);
	if(add.labels) {
		lines.ellipseDomain(c(1,2), data = h1$HL2, labels = "Hydrophobic", col = col.labels[[1]]);
		lines.ellipseDomain(c(11,16), data = h1$HL1, labels = "Cationic", col = col.labels[[2]]);
	}
	# Title:
	tmp.dy = if(is.null(dy)) 2*r[1] else dy[[1]];
	text(center[1], center[2] + tmp.dy, labels = "Piscidin 1", cex = cex.title);

	# Piscidin 2
	# MKCATLFFVLSMVVLMAEPGEG FLGRFFRRTQAILRGARQGWRAHKAVSRYRDRYIPETDNNQEQP YNQR
	x = "FLGRFFRRTQAILRGARQGWRAHKAVSRYRDRYIPETDNNQEQP"
	x = strsplit(x, "")[[1]];
	if(as.position) x = seq(length(x));
	h2 = helix.wheel(x, r = r[2], center=center, cex=cex,
		fill=fill, col.lines=col.lines, lwd.arrow=lwd.arrow, col.arrow=col.arrow);
	# Plot:
	rr = h2$HL3$C$r[[1]];
	xlim = range(h2$HL3$C$center[, 1]) + c(-rr, rr);
	if(xlim[1] > lim[1]) xlim[1] = lim[1];
	if(xlim[2] < lim[2]) xlim[2] = lim[2];
	plot.base(xlim = xlim, ylim = lim, axt = NULL);
	lines(h2);
	# Title:
	tmp.dy = if(is.null(dy)) 2*r[2] + 1 else dy[[2]];
	text(center[1], center[2] + tmp.dy, labels = "Piscidin 2", cex = cex.title);
	
	par(par.old);
	invisible(list(H1 = h1, H2 = h2));
}

#######################################
#
# Title: BioShapes
# Maintainer: L. Mada
#
# https://github.com/discoleo/BioShapes


### Genomics Tools


#' @export
genes.mutations = function(x, y, lbl = NULL, which = 1,
		nrow = c(3,3), fill = "#FFA096", cex.text = 1,
		d = 0.575, d.row = 5.5*d, d.col = 0.5,
		dna.join = c(0.25, 0.75)) {
	if(is.null(lbl)) {
		lbl = c("Wild-type", "Heterozygous", "Homozygous", "",
			"Compound\nHeterozygous\n(in trans)",
			"Compound\nHeterozygous\n(in cis)");
		if(which == 1) lbl = lbl[c(1,NA,3,NA,5)];
		if(which == 2) lbl = lbl[c(NA,2,NA,NA,NA,6)];
	}
	lst = genes.boxes(x=x, y=y, nrow=nrow, lbl=lbl,
		cex.text = cex.text, d = d, d.row = d.row, d.col = d.col,
		dna.join = dna.join);
	lbl = lst$Labels;
	lst$Labels = NULL;
	lst$Gene4  = NULL; # remove Top/Col 2;
	# Mutations:
	tM = list(
		list(ID = "Gene2", G1 = 0.75, G2 = NA),
		list(ID = "Gene3", G1 = 0.75, G2 = 0.75),
		list(ID = "Gene5", G1 = 0.75, G2 = 0.25),
		list(ID = "Gene6", G1 = c(0.25, 0.75), G2 = NA) );
	nmG = if(which == 1) "G1" else "G2";
	lsP = list();
	for(id in seq_along(tM)) {
		mM = tM[[id]];
		tt = mM[[nmG]]; tt = tt[! is.na(tt)];
		if(length(tt) > 0) {
			xy = lst[[mM$ID]]$Gene;
			pM = lapply(tt, function(ti) {
				tt = c(ti - 0.05, ti + 0.05); # Hard-coded
				x1 = xy$x[1] * (1-tt) + xy$x[2] * tt;
				x2 = xy$x[4] * (1-tt) + xy$x[3] * tt;
				y1 = xy$y[1] * (1-tt) + xy$y[2] * tt;
				y2 = xy$y[4] * (1-tt) + xy$y[3] * tt;
				r  = list(x = c(x1, x2[2:1]), y = c(y1, y2[2:1]), fill=fill);
				class(r) = c("polygon", "list");
				return(r);
			});
			lsP = c(lsP, pM);
		}
	}
	if(length(lsP) > 0) {
		lsP = as.bioshape(lsP);
		lst$M = lsP;
	}
	# Labels: add back;
	lst$Labels = lbl;
	return(lst);
}


#' @export
genes.boxes = function(x, y, nrow = c(3,2), lbl = NULL,
		cex.text = 1, d = 1, d.row = 3*d, d.col = 0.5,
		dna.join = c(0.25, 0.75)) {
	slope = slope(x, y);
	# Gene: Base Unit
	as.gene = function(x, y) {
		t0 = rev(dna.join);
		t1 = 1 - t0;
		xE = t0*x[1] + t1*x[2];
		yE = t0*y[1] + t1*y[2];
		# Gene:
		xyG = shift.ortho(xE, yE, d = c(-d, d), slope=slope);
		xyG = list(x = xyG$x[c(1,2,4,3)], y = xyG$y[c(1,2,4,3)]);
		class(xyG) = c("polygon", "list");
		# DNA Strand / Backbone:
		xyB = as.bioshape(list(
			L1 = list(x = c(x[1], xE[1]), y = c(y[1], yE[1])),
			L2 = list(x = c(x[2], xE[2]), y = c(y[2], yE[2]))
		));
		lst = list(Strand = xyB, Gene = xyG);
		return(as.bioshape(lst));
	}
	# Panels:
	nc = length(nrow);
	nr = max(nrow);
	hasRows = (nr > 1);
	hasCols = (nc > 1);
	# TODO:
	# - [basic version] Split into columns & Labels;
	# - Mutations;
	if(! hasRows) {
		return(as.gene(x, y));
	} else {
		if(nr %% 2 == 0) {
			d.rows = seq(d.row/2, by = d.row, length.out = nr / 2);
			d.rows = c(rev(d.rows), -d.rows);
		} else {
			d.rows = seq(d.row, by = d.row, length.out = (nr - 1) / 2);
			d.rows = c(rev(d.rows), 0, -d.rows);
		}
		if(hasCols) {
			dTot = dist.xy(x, y);
			dCol = (dTot - (nc-1)*d.col) / nc;
			tSep = d.col / dTot; tUnit = (dCol + d.col) / dTot;
			ti = tUnit * seq(nc - 1);
			ti = as.vector(rbind(ti - tSep, ti));
			xi = x[1] * (1 - ti) + x[2] * ti;
			yi = y[1] * (1 - ti) + y[2] * ti;
			xi = c(x[1], xi, x[2]);
			yi = c(y[1], yi, y[2]);
			lst = list();
			for(idU in seq(0, nc-1)) {
				id  = 2*idU; id = c(id + 1, id + 2);
				xy  = shift.ortho(xi[id], yi[id], slope=slope,
					d = d.rows, simplify = FALSE);
				tmp = lapply(xy, function(xy) {
					as.gene(xy$x, xy$y);
				});
				lst = c(lst, tmp);
			}
		} else {
			xy  = shift.ortho(x, y, slope=slope, d = d.rows, simplify = FALSE);
			lst = lapply(xy, function(xy) {
				as.gene(xy$x, xy$y);
			});
		}
		names(lst) = paste0("Gene", seq_along(lst));
		# Labels:
		hasLbl = any(! is.na(lbl));
		if(hasLbl) {
			cc = sapply(lst, function(obj) {
				center.xy(obj$Gene$x, obj$Gene$y);
			});
			len = length(lst);
			txt = list(x = cc[1,], y = cc[2,], labels = lbl[seq(len)],
				cex = cex.text);
			class(txt) = c("text", "list");
			lst$Labels = txt;
		}
		return(as.bioshape(lst));
	}
}


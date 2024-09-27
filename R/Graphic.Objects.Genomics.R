#######################################
#
# Title: BioShapes
# Maintainer: L. Mada
#
# https://github.com/discoleo/BioShapes


### Genomics Tools


#' @export
genes.mutations = function(x, y, nrow = c(3,2), lbl = NULL,
		cex.text = 1, d = 1, d.row = 3*d, d.col = 0.5,
		dna.join = c(0.25, 0.75)) {
	if(is.null(lbl)) {
		lbl = c("Wild-type", "Heterozygous", "Homozygous",
			"Compound\nHeterozygous\n(in trans)",
			"Compound\nHeterozygous\n(in cis)");
	}
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
	# - Split into columns;
	# - Labels & mutations;
	if(! hasRows) {
		return(as.gene(x, y));
	} else {
		if(nr %% 2 == 0) {
			d.row = seq(d.row/2, by = d.row, length.out = nr / 2);
			d.row = c(rev(d.row), -d.row);
		} else {
			d.row = seq(d.row, by = d.row, length.out = (nr - 1) / 2);
			d.row = c(rev(d.row), 0, -d.row);
		}
		xy  = shift.ortho(x, y, d = d.row, simplify = FALSE);
		lst = lapply(xy, function(xy) {
			as.gene(xy$x, xy$y);
		});
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


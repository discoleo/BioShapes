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


### Chem; Analysis

# - various Tools for structure analysis;


##################

#' @export
proj.newman = function(lbl1, lbl2, phi = c(0, pi/length(lbl1)), r = 1, center = c(4,4),
		pos = NULL, dr = r/4, lwd = c(2, 1)) {
	if(missing(lbl2)) {
		if(length(lbl1) > 1) stop("Missing 2nd group of ligands!");
		lbl1 = strsplit(lbl1, "\\||\n");
		lbl1 = lbl1[[1]];
		if(length(lbl1) != 2) stop("Invalid number of ligands!");
		lbl1 = strsplit(lbl1, ",");
		lbl2 = lbl1[[2]]; lbl1 = lbl1[[1]];
		
	}
	# Circle:
	lst = list(r=r, center=center);
	class(lst) = c("circle", "list");
	# Bonds:
	doDF = function(l1, l2, offset = 0) {
		l1 = as.data.frame(l1);
		l2 = as.data.frame(l2);
		l1$id = seq(1 + offset, length.out = nrow(l1));
		l2$id = seq(1 + offset, length.out = nrow(l2));
		rbind(l1, l2);
	}
	rB = r + dr;
	n1 = length(lbl1);
	n2 = length(lbl2);
	b10 = list(x = rep(center[1], n1), y = rep(center[2], n1));
	b1E = points.circle(n = n1, r = rB, phi = phi[1], center=center);
	b20 = points.circle(n = n2, r = r , phi = phi[2], center=center);
	b2E = points.circle(n = n2, r = rB, phi = phi[2], center=center);
	bb1 = doDF(b10, b1E);
	bb2 = doDF(b20, b2E, offset = n1);
	lst = list(C = lst,
		BB1 = as.bioshape(list(bb1, lwd = lwd[1])),
		BB2 = as.bioshape(list(bb2, lwd = lwd[2])));
	# Labels:
	if(is.null(pos)) {
		# TODO:
	} else {
		lT1 = lapply(seq(n1), function(id) {
			lst = list(x = b1E$x[id], y = b1E$y[id], labels = lbl1[id], pos = pos[[1]][id]);
			class(lst) = c("text", "list");
			return(lst);
		});
		lT2 = lapply(seq(n2), function(id) {
			lst = list(x = b2E$x[id], y = b2E$y[id], labels = lbl2[id], pos = pos[[2]][id]);
			class(lst) = c("text", "list");
			return(lst);
		});
		lT = c(lT1, lT2);
		lst$T = as.bioshape(lT);
	}
	return(as.bioshape(lst));
}

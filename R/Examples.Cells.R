###################
#
# Title: BioShapes
#
# Maintainer: L. Mada
#
#
# Continuation of:
# 1. Bachelor Thesis: Adrian Cotoc (BSc 2022-2023)
# Faculty of Mathematics and Informatics, UVT
#
# Coordinator:
#   Prof. Daniela Zaharie
#   Dr. med. Leonard Mada (Syonic SRL)
#   in collaboration with Syonic SRL
#
# 2. Bachelor Thesis: Darian Voda (BSc 2021-2022)
#
# GitHub: https://github.com/discoleo/BioShapes
# [old]
# GitHub: https://github.com/Adi131313/BioShapes


### Cells: Examples

### Muscle Tissue
#' @export
example.muscle = function(n = 6, fill = "red"){
  tmp = muscle(n = n, fill = fill)
  plot.base()
  # TODO: implement fill in muscle/lens
  lines(tmp, fill=fill)
}


### 3 Blood cells
#' @export
example.bloodCells = function(radius = 2, lwd = 10){
  radius = 2;
  lim = c(-radius, radius) * 4;
  plot.base(xlim = lim, ylim = lim)
  tmp = draw_blood_cell(radius = radius, center = c(-3, -1), lwd = lwd)
  lines(tmp)
  tmp = draw_blood_cell(radius = radius, center = c(3, 1), lwd = lwd)
  lines(tmp)
  tmp = draw_blood_cell(radius = radius, center = c(3, -5), lwd = lwd)
  lines(tmp)
}

#' @export
example.cell.phyllopodia = function(n = 12, scale.base = c(1,1.4,1,1),
		lwd = c(1,2,1,1), phi = c(0, 0, 0, pi*7/12), r = 3) {
	if(length(n) == 1) n = rep(n, 4);
	if(length(r) == 1) r = rep(r, 4);
	if(length(phi) == 1) phi = rep(phi, 4);
	if(length(lwd) == 1) lwd = rep(lwd, 4);
	par.old = par(mfrow = c(2,2));
	
	tmp = cell.podia(r = r[1], n = n[1], center = c(4,4), phi = phi[1],
		scale.base = scale.base[1], lwd = lwd[1])
	plot.base(axt = NULL)
	lines(tmp)
	
	# Length + Base
	tmp = cell.podia(r = r[2], n = n[2], center = c(4,4), phi = phi[2],
		scale.base = scale.base[2], l.scale = seq(0.4,1, length=n[2]),
		lwd = lwd[2], col = "#90C020")
	plot.base(axt = NULL)
	lines(tmp)
	
	# Inward: Strange
	tmp = cell.podia(r = r[3], n = n[3], center = c(4,4), phi = phi[3], l.scale = -0.4,
		scale.base = scale.base[3], lwd = lwd[3])
	plot.base(axt = NULL)
	lines(tmp)
	
	# Fill
	tmp = cell.podia(r = r[4], n = n[4], center = c(4,4), phi = phi[4],
		scale.base = scale.base[4], l.scale = seq(0.4,1, length=n[4]),
		lwd = lwd[4], col = "#90C020", fill = "#F0A0B0")
	plot.base(axt = NULL)
	lines(tmp)
	
	par(par.old);
	invisible();
}


##############
### Neuron ###

### Design of Neuron
#' @export
example.neuronDesign = function(n = 5, r = 2, col.mark = "red") {
  if(length(r) == 1) r = rep(r, 2);
  par.old = par(mfrow = c(1,2));

  phi = 2*pi/n;
  tmp = circles.OnCircle(n, r[1], phi = pi/n);
  test.FilledCircle(tmp, lim = c(-6, 6), pin = FALSE);
  R = attr(tmp, "R");
  tmp = sapply(seq(n), function(k) points(
    r[1] * cos(pi/2 + k * phi) + R * cos(phi * (k-1) + phi/2),
    r[1] * sin(pi/2 + k * phi) + R * sin(phi * (k-1) + phi/2),
	col = col.mark) );

  phi = 0; center = c(2, 3);
  plot.base()
  tmp = neuron(n = n, r = r[2], center = center, phi = phi)
  lines(tmp)

  par(par.old);
  invisible();
}

### Neuron
#' @export
example.neuron = function(phi = 0, n = 5, type = "Tree",
		col = NULL, fill.nucl = NULL, center = c(2, 3)) {
  plot.base()
  tmp = neuron(n = n, center = center, phi = phi, type.syn = type,
	col = col, fill.nucl = fill.nucl);
  lines(tmp);
  invisible(tmp);
}

#' @export
description.neuron = function(phi = 0,
			lbl = c("Axon", "Dendrites", "Nucleus"),
			title = "Neuron",
			lwd = 2, col = "#48B000", col.arrow = "red",
			d.arrow = -0.4, cex.title = 1.5, xy.title = c(3, 8)){
  # TODO: parameters
  neuron = example.neuron(phi = phi);

  # Title
  if( ! is.null(title)) text(xy.title[1], xy.title[2], title, cex=cex.title);

  # Labels
  # TODO: Case phi != 0
  a1 = arrowSimple(x=c(5,6), y=c(5.5, 3.5), d = d.arrow,
	lwd=lwd, col = col.arrow);
  text(5, 6, lbl[[1]])

  a2 = arrowSimple(x=c(1, 2), y=c(0.75, 1.5), d = d.arrow,
	lwd=lwd, col = col.arrow);
  text(1, 0, lbl[[2]])

  return(invisible());
}

### Multiple neurons
#' @export
example.neurons = function(n = 5, col = "red",
		lwd = 1, lwd.axon = lwd, lwd.nucl = 1) {
  ###
  center = c(2, 3)
  lim = c(-15, 15);
  plot.base(xlim = lim, ylim = lim, axt = NULL)

  ### Bottom Row:
  
  phi = -pi/2;
  center2 = center + c(-6, -10)
  tmp = neuron(n = n, center = center2, phi = phi, type = "SArrow",
	col=col, lwd=lwd, lwd.axon = lwd.axon, lwd.nucl = lwd.nucl, synapse = list(l = 2));
  lines(tmp)

  ### Middle Row:
  phi = -pi/1.5;
  center2 = center + c(-1, -2)
  tmp = neuron(n = n, center = center2, phi = phi, type = "Tree",
	col=col, lwd=lwd, lwd.axon = lwd.axon, lwd.nucl = lwd.nucl);
  lines(tmp)

  ###
  phi = 2*pi - pi/3;
  center2 = center + c(-11, -2)
  tmp = neuron(n = n, center = center2, phi = phi, type = "Tree",
	col=col, lwd=lwd, lwd.axon = lwd.axon, lwd.nucl = lwd.nucl);
  lines(tmp)

  ### Top Rows: The Square
  
  # TODO: All synapse types
  phi = 0;
  center = c(-9, 13)
  tmp = neuron(n = n, center = center, phi = phi, type = "Tree",
	col=col, lwd=lwd, lwd.axon = lwd.axon, lwd.nucl = lwd.nucl);
  lines(tmp)
  ###
  phi = pi/2;
  center = c(-10, 4)
  tmp = neuron(n = n, center = center, phi = phi, type = "Solid",
	col=col, lwd=lwd, lwd.axon = lwd.axon, lwd.nucl = lwd.nucl,
	synapse = list(fill = "#E08064"));
  lines(tmp)
  ###
  phi = -pi;
  center = c(0, 4)
  tmp = neuron(n = n, center = center, phi = phi, type = "Detail",
	col=col, lwd=lwd, lwd.axon = lwd.axon, lwd.nucl = lwd.nucl);
  lines(tmp)
  ###
  phi = -pi/2;
  center = c(1, 13.5)
  tmp = neuron(n = n, center = center, phi = phi, type = "Tree",
	col=col, lwd=lwd, lwd.axon = lwd.axon, lwd.nucl = lwd.nucl);
  lines(tmp)
}


###################

### Epithelia

### Brush-Border Epithelium
#' @export
example.Epithelium.BB = function(x = c(0,8), y = c(1,1), n = 8, dy = 6, h = ~2,
		col = 1, col.nc = 1, fill = "#FFFF80", fill.nc = "#9648C0") {
	plot.base();
	# Ep 1:
	ep = epithelium.brush(x=x, y=y, n=n, h=h, col=col, col.nc=col.nc,
		fill=fill, fill.nc=fill.nc);
	lines(ep);
	if(is.null(dy)) return(invisible(ep));
	# Ep 2:
	h = if(inherits(h, "formula")) formula.neg(h) else - h;
	xy  = shift.ortho(x, y, d = dy);
	ep2 = epithelium.brush(x = xy$x, y = xy$y, n=n, h=h, col=col, col.nc=col.nc,
		fill=fill, fill.nc=fill.nc);
	lines(ep2);
	#
	ep = list(Ep1 = ep, Ep2 = ep2)
	invisible(as.bioshape(ep));
}

#####################
#####################

### Tumour

#' @export
test.tumor.mass2 = function(
		d = list(
			c(0.5, 0.43), c(0.5, 0.43),
			c(0.51, 0.43), c(0.53, 0.43)), by = 3,
		lim = c(-3,3), axt = NULL) {
	nc = ceiling(length(d) / 2);
	par.old = par(mfrow = c(2,nc))
	for(id in seq(length(d))) {
		plot.base(xlim=lim, ylim=lim, axt=axt);
		lines(tumor.mix(d = d[[id]], by=by))
	}
	par(par.old)
}


#####################

### Phages

# id.faces = number the faces (only for phage 1);
#' @export
test.phage = function(phi = c(0, pi/5, -pi/3, pi/2, -pi/2),
		r = 1, col = NULL, id.faces = TRUE) {
	fill = rep("#A0A06464", 10);
	fill[c(7,8)] = "#64320064";
	ph0 = phage.base(0, 7, phi = phi[1], r=r, col=col, id.faces=id.faces);
	# Plot:
	plot.base(axt = NULL);
	lines(ph0);
	lines(phage.base(2.5, 4, phi = phi[2], fill=fill, r=r));
	lines(phage.base(4, 1, phi = phi[3], fill=fill, r=r));
	#
	lines(phage.base(4.5, 7, phi = phi[4], fill=fill, r=r));
	lines(phage.base(7.25, 7, phi = phi[5], fill=fill, r=r));
}

###################
#
# Title: BioShapes
#
# Maintainer: L. Mada
#
# Based on: Bachelor Thesis (2022-2023)
# Adrian Cotoc
# Faculty of Mathematics and Informatics, UVT
#
# Coordinator:
#   Prof. Daniela Zaharie
#   Dr. med. Leonard Mada (Syonic SRL)
#
# in collaboration with Syonic SRL
# continous the work of Darian Voda (BSc 2021-2022)
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


##############
### Neuron ###

### Design of Neuron
#' @export
example.neuronDesign = function(n = 5, r = 2, col.mark = "red") {
  if(length(r) == 1) r = rep(r, 2);
  par.old = par(mfrow = c(1,2));

  phi = 2*pi/n;
  tmp = circlesOnCircle(n, r[1], phi = pi/n);
  testFilledCircle(tmp, lim = c(-6, 6), pin = FALSE);
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
  tmp = neuron(n = n, center = center, phi = phi, type = type,
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

  ### First group
  
  phi = -pi/2;
  center2 = center + c(-6, -10)
  tmp = neuron(n = n, center = center2, phi = phi, type = "SArrow",
	col=col, lwd=lwd, lwd.axon = lwd.axon, lwd.nucl = lwd.nucl);
  lines(tmp)

  ###
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

  ### Second group
  
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
	fill.synapse = "#E08064");
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

#######################################
#
# Title: BioShapes
#
# Maintainer: L. Mada
#
# https://github.com/discoleo/BioShapes
#
# Continuation of:
# 1. Bachelor Thesis: Adrian Cotoc (2022-2023)
# Faculty of Mathematics and Informatics, UVT
#
# Coordinator:
#   Prof. Daniela Zaharie
#   Dr. med. Leonard Mada (Syonic SRL)
#   in collaboration with Syonic SRL
# GitHub: https://github.com/Adi131313/BioShapes
#
# 2. Bachelor Thesis: Darian Voda (2021-2022)



### Functions to Generate Neuron

### Synapse Object
#' @export
as.synapse = function(x) {
	class(x) = c("synapse", "bioshape", "list");
	return(x);
}


### Neuron: Main function
#' @export
neuron = function(center = c(0, 0), n = 5, r = 2, phi = 0,
			axon.length = 3 * r, dendrite.length = ~ r/2, r.nucl = ~ (R - r)/2,
			type.syn = c("Solid", "Tree", "Arrow", "SArrow",
				"Detail", "Radial", "T"),
			lwd = 1, lwd.axon = lwd, lwd.dendrite = lwd, lwd.nucl = 1,
			col = 1, col.nucl = 1, fill.nucl = NULL,
			col.dendrite = col, col.axon = col,
			synapse = list()) {
	body = neuron.body(center = center, n = n, r = r, phi = phi,
		lwd = lwd, col = col);
	### Init
	axon.length = axon.length; # force = scale * r;
	R = r;
	r = attr(body, "r");
	phin = 2 * pi/n;
	phi0 = - phin/2;
	pi20 = pi/2; pi32 = 3*pi/2;
	phiD = seq(n) * phin + phi;
	phiR = as.radians0(phiD);
	x0 = r * cos(phiD + pi20) + R * cos(phiD + phi0);
	y0 = r * sin(phiD + pi20) + R * sin(phiD + phi0);
	x0 = x0 + center[1];
	y0 = y0 + center[2];
	sg = ifelse(phiR > pi/2 & phiR <= pi32, - 1, 1);
	### Axon
  if(axon.length != 0) {
    phiA  = phiR[n];
    slope = if(abs(phiA - pi20) < 1E-2 || abs(phiA - pi32) < 1E-2) Inf
    else tan(phiA);
    xy = shiftPoint(c(x0[n], y0[n]), d = sg[n] * axon.length, slope = slope);
    axon = list(x = c(x0[n], xy[1]), y = c(y0[n], xy[2]),
		lwd = lwd.axon, col = col.axon);
    axon = list(axon);
    lenDendites = n - 1;
    if( ! is.null(type.syn)) {
      if(inherits(synapse, "synapse")) {
	    syn = synapse;
	  } else if(is.function(synapse)) {
		syn = synapse(xy, phi = phiA);
	  } else {
		if(is.null(synapse$l)) synapse$l = 2/3 * r;
		if(is.null(synapse$col)) synapse$col = col;
		if(is.null(synapse$fill)) synapse$fill = col;
        synapse = c(list(p = xy, phi = phiA, type = type.syn), synapse);
        syn = do.call("synapse", synapse);
	  }
      axon = c(Axon = axon, Synapse = syn);
    }
  } else {
    axon = NULL;
    lenDendites = n;
  }
	### Dendrites
	eval.formula = function(x) {
		xx = if(inherits(x, "formula")) {
			eval(x[[2]]);
		} else x;
	}
	dendrite.len = eval.formula(dendrite.length);
	dendrite.len = sg * dendrite.len;
	dend = lapply(seq(lenDendites), function(k) {
		tree(c(x0[k], y0[k]), d = dendrite.len[k], slope = tan(phiD[k]),
			n = 2, levels = 2); # TODO
	})
	dend$lwd = lwd.dendrite;
	dend$col = col.dendrite;
	dend = as.bioshape(dend);
	dend = list(dend);
	### Nucleus
	r.nucl = eval.formula(r.nucl);
	if(r.nucl > 0) {
		nucleus = list(r = r.nucl, center = center, lwd = lwd.nucl,
			col = col.nucl, fill = fill.nucl);
		class(nucleus) = c("circle", "list");
		nucleus = list(nucleus);
	} else nucleus = NULL;
	
	### Neuron
	neuron = c(Body = body, axon,
		Dendrites = dend, Nucleus = nucleus);
	class(neuron) = c("neuron", "list");
	return(as.bioshape(neuron));
}

#' @export
tree = function(p, d, slope, n=2, levels=2) {
  xy = shiftPoint(p, d=d, slope = slope);
  xy = list(x = c(p[1], xy[1]), y = c(p[2], xy[2]));
  return(xy);
#
  xy = shiftPoint(p, d=d, slope = slope)
  x = c(p[1], xy[1])
  y = c(p[2], xy[2])
  if (n > 2) {
    control_points = matrix(nrow = n-2, ncol = 2)
    for (i in 2:(n-1)) {
      control_points[i-1,] = (p + xy)/2 + c(-1, 1)*abs((xy-p)[2]/10)*(-1)^i
    }
    bezier_points = bezier(x, y, control_points = control_points)
    x = bezier_points$x
    y = bezier_points$y
  }
  xy = list(x = x, y = y)
  return(xy)
}

#' @export
neuron.body = function(center = c(0, 0), n = 5, r = 3, phi = 0,
		lwd = 1, col = 1) {
  phi0 = phi + pi/n;
  cc = circlesOnFixedCircle(n = n, r = r, center = center, phi = phi0);
  R = r;
  r = attr(cc, "r");
  phin = 2*pi/n;
  id = seq(n);
  a1 = pi/2 + phi + id * phin;
  a2 = pi/2 + phi + (id + (n - 2)/2) * phin;
  lst = lapply(id, function(id) {
    lst = list(
      r = r, center = c(x = cc$x[id], y = cc$y[id]),
      phi = c(a1[id], a2[id]) );
	lst$lwd = lwd;
    lst$col = col;
    class(lst) = c("circle.arc", "list");
    return(lst);
  });
  attr(lst, "r") = r;
  class(lst) = c("bioshape", "list");
  return(lst);
}

# phi = angle of axon (in radians);
# slope = slope of axon, as alternative to angle phi;
# alpha = angle of tree-cone (in degrees);
#' @export
synapse = function(p, phi, type = c("Solid", "Tree", "Arrow", "SArrow",
			"Detail", "Radial", "T"),
			col=1, fill="#808080",
			slope = NULL, l = 1/2, n = 4, alpha = 120,
			helix.scale = 1/12, ...) {
  type = match.arg(type);
  if( ! is.null(slope)) {
    th  = atan(slope);
    phi = alpha * pi / 360;
  } else {
    th  = phi;
    phi = alpha * pi / 360;
  }
  if(type == "Detail") {
    if(is.null(slope)) slope = tan(th);
    sg = th > pi/2 & th <= 3*pi/2;
    d  = if(sg) -l else l;
    cc = shiftPoint(p, d = d, slope=slope);
    th = th + pi + c(- phi, phi); # already 1/2;
    xy = circle.arc(r = l, th, center = cc);
    len = length(xy$x);
    p2 = c(xy$x[1], xy$y[1]);
    p1 = c(xy$x[len], xy$y[len]);
    lst = helix(p1, p2, n = n + 0.5, A = l * helix.scale, phi=0, N=129);
    lst = lst[[1]];
    lst$x = c(xy$x, lst$x);
    lst$y = c(xy$y, lst$y);
    lst$col = col; lst$fill = fill;
    class(lst) = c("polygon", "list")
    lst = list(lst);
    return(as.synapse(lst));
  }
  # Tree Types
  if(type == "T") {
	# skip;
  } else if(phi >= 0) {
    th = th + c(- phi, phi);
  } else {
    th = th + pi + c(- phi, phi);
  }
  if(type == "Tree") {
    px = l * cos(th) + p[1];
    py = l * sin(th) + p[2];
    id = seq(0, 1, length.out = n);
    px = id * px[1] + (1 - id) *  px[2];
    py = id * py[1] + (1 - id) *  py[2];
    xy = data.frame(x = p[1], y = p[2], id = seq(n));
    xy = rbind(xy, data.frame(x = px, y = py, id = seq(n)));
  } else if(type == "Radial") {
    id = seq(th[1], th[2], length.out = n);
    px = l * cos(id) + p[1];
    py = l * sin(id) + p[2];
    xy = data.frame(x = p[1], y = p[2], id = seq(n));
    xy = rbind(xy, data.frame(x = px, y = py, id = seq(n)));
  } else if(type == "T") {
    id = c(-pi, pi)/2 + th;
	l  = l/2; # uses L / 2 !!!
    px = l * cos(id) + p[1];
    py = l * sin(id) + p[2];
	xy = list(x = px, y = py);
  } else if(type == "Arrow" || type == "SArrow") {
    id = th + pi;
	l  = l/2; # uses L / 2 !!!
    px = l * cos(id) + p[1];
    py = l * sin(id) + p[2];
    if(type == "Arrow") {
		xy = list(
			x = c(px[1], p[1], px[2]),
			y = c(py[1], p[2], py[2]), col=col);
	} else {
		xy = list(x = c(p[1], px), y = c(p[2], py), col=col, fill=fill);
		class(xy) = c("polygon", "list");
	}
  } else if(type == "Solid") {
    px = l * cos(th) + p[1];
    py = l * sin(th) + p[2];
    xy = list(x = c(p[1], px), y = c(p[2], py), col=col, fill=fill);
    class(xy) = c("polygon", "list");
  } else stop("Not yet implemented!");
  xy$col = col;
  xy = as.synapse(list(xy));
  return(xy);
}

#' @export
synapse.neuron = function(x, type = c("Solid", "Tree", "Arrow", "SArrow",
			"Detail", "Radial", "T"), col=1, fill="#808080",
			l = 1/2, n = 4, alpha = 120, helix.scale = 1/12, ...) {
	A = x$Axon;
	phi = atan2(A$y[2] - A$y[1], A$x[2] - A$x[1]);
	len = length(A$x);
	syn = synapse(c(A$x[len], A$y[len]), phi=phi, type=type, col=col, fill=fill,
		l=l, n=n, alpha=alpha, helix.scale=helix.scale, ...);
	return(syn);
}

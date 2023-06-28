###################
#
# Bachelor Thesis
#
# Title: BioShapes
#
# Candidate: Adrian Cotoc
# Faculty of Mathematics and Informatics, UVT
#
# Coordinator:
#   Prof. Daniela Zaharie
#   Dr. med. Leonard Mada (Syonic SRL)
#
# in collaboration with Syonic SRL
# continous the work of Darian Voda
#
# GitHub: https://github.com/Adi131313/BioShapes

### Functions to Generate Objects


# Generates a star polygon
#' @export
star = function(n, R = c(2,1), center = c(0,0), lwd=1, phi = pi/2,
                col=1, fill = NULL) {
  c1 = pointsCircle(n, r=R[1], center=center, phi=phi[1]);
  c2 = pointsCircle(n, r=R[2], center=center, phi=phi[1] + pi/n);
  x = as.vector(rbind(c1$x, c2$x));
  y = as.vector(rbind(c1$y, c2$y));
  lst = list(x = x, y = y, col = col, lwd = lwd);
  if(! is.null(fill)){
    lst$fill = fill;
  }
  class(lst) = c("polygon", "list");
  lst = list(lst);
  class(lst) = c("bioshape", "list");
  return(invisible(lst));
}

#' @export
ngon = function(n, r = 1, center = c(0, 0), phi = 0,
		lwd = NULL, col = NULL, fill = NULL) {
  p = pointsCircle(n=n, r=r, center=center, phi=phi);
  if( ! is.null(lwd)) p$lwd = lwd;
  if( ! is.null(col)) p$col = col;
  if( ! is.null(fill)) p$fill = fill;
  class(p) = c("polygon", "list");
  return(p);
}

#' @export
ngon.circle = function(n, N, R = 2, r = 1/2, center = c(0, 0), phi = 0,
		lwd = NULL, col = NULL, fill = NULL){
  cp = pointsCircle(n = N, r = R, center=center, phi=phi);
  pphi = seq(0, N-1) * 2*pi/N + phi;
  pg = lapply(seq(N), function(id) {
	cc = c(cp$x[id], cp$y[id]);
    ngon(n, r = r, center = cc, phi = pphi[id], lwd=lwd, col=col, fill=fill);
  })
  class(pg) = c("bioshape", "list");
  invisible(pg);
}

#### Liposomes ####
#' @export
liposomes = function(n, r, center=c(0, 0), phi=c(0, 0), d=0, ...){
  C1 = circlesOnCircle(n=n[1], r=r, center=center, phi=phi[1])
  C2 = circlesOnCircle(n=n[2], r=r, center=center, phi=phi[2])
  R1 = attr(C1, "R")
  R2 = attr(C2, "R")
  R1 = R1 - r
  R2 = R2 + r
  d2 = (R1-R2-d)/2
  p1 = pointsCircle(n=n[1], r=R1, phi=phi[1], center=center)
  p2 = pointsCircle(n=n[2], r=R2, phi=phi[2], center=center)
  fn = function(id, p1, p2, d){
    p1 = c(p1$x[id], p1$y[id])
    slope = slope(x=c(p1[1], p2[1]), y=c(p1[2], p2[2]))
    if(p1[1] > center[1]){
      d = -d;
    }
    pp = shiftPoint(p1, slope=slope, d=d)
    data.frame(x=c(p1[1], pp[1]), y=c(p1[2], pp[2]), id=id)
  }
  # Lines: Side Chains
  l1 = lapply(seq(n[1]), fn, p1, center, d2);
  l2 = lapply(seq(n[2]), fn, p2, center, -d2);
  l1 = do.call(rbind, l1)
  l2 = do.call(rbind, l2)
  l2$id = l2$id + nrow(l1)
  l = rbind(l1, l2);
  # Liposome:
  obj = list(C1=C1, C2=C2, l=l);
  class(obj) = c("liposome", "list");
  return(obj);
}

### Ring ###
#' @export

ring = function(n, R = c(5, 7), center = c(0,0),
                col = 1, fill = NULL, phi = c(0,0)){
  return(duct(n = n, R = R, nc.r = NULL, center = center,
              col = col, fill = fill, phi = phi))
}

### Glandular Duct ###
#' @export
duct = function(n, R = c(5, 7), nc.r=1/2, center=c(0,0),
                col = 1, fill = NULL,  nc.fill = "#E03232", phi=c(0,0)) {

  if(length(phi) == 1) phi = c(phi, phi);
  c1 = pointsCircle(n, r=R[1], center=center, phi=phi[1]);
  c2 = pointsCircle(n, r=R[2], center=center, phi=phi[2]);

  ### X
  cbx = rbind(c2$x, c1$x)
  cbx = cbind(cbx[, -1], cbx[,1])

  cells.x = rbind(c1$x, c2$x, cbx, c1$x)
  cells.x = lapply(seq(ncol(cells.x)), function(nc){
    cells.x[, nc]
  })


  ### Y
  cby = rbind(c2$y, c1$y)
  cby = cbind(cby[, -1], cby[,1])

  cells.y = rbind(c1$y, c2$y, cby, c1$y)
  cells.y = lapply(seq(ncol(cells.y)), function(nc){
    cells.y[, nc]
  })

  cells = lapply(seq(n), function(id){
    l = list(x = cells.x[[id]], y = cells.y[[id]])
    l$col = col; l$fill = fill;
    class(l) = c("polygon", "list")
    return(l)
  })

  # plot.base(xlim=c(-10,10), ylim=c(-10,10))
  # polygon(cells.x, cells.y)

  ### Nuclei:
  if(is.null(nc.r)){
    class(cells) = c("bioshape", "list")
    return(cells);
  }

  shift = function(x) c(x[-1], x[1]);

  mid1.x = (c1$x + shift(c1$x))/2;
  mid2.x = (c2$x + shift(c2$x))/2;
  mid.x = (mid1.x + mid2.x)/2;

  mid1.y = (c1$y + shift(c1$y))/2;
  mid2.y = (c2$y + shift(c2$y))/2;
  mid.y = (mid1.y + mid2.y)/2;
  # Centers:
  nuclei = list(center = cbind(mid.x, mid.y), r = nc.r, fill = nc.fill);
  class(nuclei) = c("circle", "list")
  cells = c(cells, list(nuclei));
  # TODO
  # testFilledCircle(nuclei,r=nc.r, add=TRUE, line=FALSE)
  class(cells) = c("bioshape", "list")
  return(cells);
}

# r[1] = outer radius;
# r[c(2, 3)] = outer & inner radius of cell component;
#' @export
duct.complex = function(center = c(0, 0), r = c(7, 5, 2), n = 8,
                        h.scale = 1, R.scale = c(1, 1), dr = 0.25) {
  radius = r;
  # Basement Membrane
  lst = draw_blood_cell(radius = radius[1], center = center,
                        col = "#bf9d9b", fill = "#f0e19e", lwd = 7);
  tmp = list(r = radius[1], center = center, col = "white", lwd = 5);
  class(tmp) = c("circle", "list");
  lst = c(lst, list(tmp));

  # Duct
  tmp = duct(n, radius[c(2,3)], center = center, phi=pi/2/n,
             fill = "#f0b25b", nc.fill = "#6f618f");
  lst = c(lst, list(tmp));

  # Lens
  h = radius[2] * sin(pi/n) * h.scale;
  h = c(h, h);
  pos = c(0, 1);
  n2 = n/2; ninv = 1/(2*n);
  x = (radius[2] + dr) * cos((seq(n)/n2 - ninv) * pi) + center[1];
  y = (radius[2] + dr) * sin((seq(n)/n2 - ninv) * pi) + center[2];
  for(o in seq(n)) {
    if(o <= n2) {
      tmp = lens.group(x=c(x[o], x[o+n2]), y=c(y[o], y[o+n2]),
                       h=h, pos=pos, l.scale = R.scale, fill = "#FF3220");
      lst = c(lst, tmp); # list(tmp)
    }
    tmp = list(center = c(x[o],y[o]), r = 0.2, fill = "#8a4965")
    class(tmp) = c("circle", "list");
    lst = c(lst, list(tmp));
  }
  lst = as.bioshape(lst);
  invisible(lst);
}

# Generates a virus
#' @export
virus = function(R = 1, center = c(0,0), n.spike = 10, off.spike = c(0, 1),
                 r.spike=0.4, ngon.spike=4, phi.spike = 0, lwd = 6, lwd.spike = lwd/2,
                 col = "#23b5cc", col.spike = "#a60845", fill.spike = "#f26d07") {
  ### Spikes
  virus = spikes(R = R, center = center, n.spike = n.spike, off.spike = off.spike,
                  r.spike=r.spike, ngon.spike = ngon.spike[1], phi.spike = 0,
                  lwd = lwd.spike, lwd.spike = lwd.spike,
                  col = col.spike, fill = fill.spike);

  ### Envelope
  vc = list(r = R, center = center, lwd = lwd, col = col);
  class(vc) = c("circle", "list");
  virus$virus = vc;

  class(virus) = c("bioshape", "list");
  return(invisible(virus));
}

# off.spike = starting & ending of the spike;
#' @export
spikes = function(R = 1, center = c(0,0), n.spike = 10, off.spike = c(0, 1),
                  r.spike=0.25, ngon.spike=4, phi.spike = 0, lwd = 1, lwd.spike = 2*lwd,
                  col = "#D06432", fill = NULL) {
  # TODO: rename to r.ngon;
  r.ngon = r.spike;
  ### Spikes
  c1 = pointsCircle(n.spike, r = R + off.spike[1], center = center, phi = phi.spike);
  c2 = pointsCircle(n.spike, r = R + off.spike[2], center = center, phi = phi.spike);
  spike = lapply(seq(n.spike), function(id) {
    list(x = c(c1$x[id], c2$x[id]), y = c(c1$y[id], c2$y[id]));
  });
  spike$lwd = lwd.spike;
  spike$col = col;
  class(spike) = c("lines.list", "list");
  virus = list(spike);
  ### Head of Spikes
  if( ! is.null(r.ngon) && r.ngon > 0 ) {
    R3 = R + off.spike[2] + r.ngon;
    if(ngon.spike > 0) {
      h = ngon.circle(ngon.spike, N = n.spike, r = r.ngon, R = R3,
                      center = center, phi = phi.spike, fill = fill);
      virus$head = h;
    } else {
      cc = pointsCircle(n.spike, r = R3,
                        center = center, phi = phi.spike);
      cc = cbind(x = cc$x, y = cc$y);
      h  = list(r = r.ngon, center = cc,
                lwd = lwd, col = col, fill = fill);
      class(h) = c("circle", "list");
      virus$head = h;
    }
  }
  class(virus) = c("bioshape", "list");
  return(invisible(virus));
}

virus2 = function(R = 2, center = c(0,0), n.spike = 10, off.spike = c(0, 1, 1.5),
                  r.spike=0.25, ngon.spike=c(4, 0), phi.spike = 0,
                  lwd = 8, lwd.spike = lwd/2, col = "#D06432",
                  col.spike = c("#c90e40", "red"),
                  fill.spike = c("#5128c9", "yellow")) {

  v1 = virus(R = R, center = center, n.spike = n.spike, off.spike = off.spike,
             r.spike=r.spike, ngon.spike = ngon.spike[1], phi.spike = phi.spike,
             lwd = lwd, lwd.spike = lwd.spike[1],
             col = col, col.spike = col.spike[1], fill.spike = fill.spike[1])

  ### Spikes
  if(length(lwd.spike) <= 1) lwd.spike = rep(lwd.spike[1], 2);
  off.spike0 = off.spike[1];
  off.spike_i = c(off.spike0, off.spike[3]);
  spikes2 = spikes(R = R, center = center, n.spike = n.spike,
                   off.spike = off.spike_i, r.spike=r.spike,
                   ngon.spike = ngon.spike[2], phi.spike = phi.spike + pi/n.spike,
                   lwd = lwd.spike[2], lwd.spike = lwd.spike[2],
                   col = col.spike[2], fill = fill.spike[2]);

  virus = c(spikes2, v1);
  virus = as.bioshape(virus);
  return(virus)
}

### Convex lens
#' @export
lens = function(x, y, R = NULL, scale = c(1,1),
                lwd=1, col=1, fill = NULL) {
  d2 = (x[2] - x[1])^2 + (y[2] - y[1])^2;
  d  = sqrt(d2);
  ### Lens Radius
  if(is.null(R)) {
    R = d;
    R = c(R, R);
  } else if(length(R) == 1) R = c(R, R);
  R = R * scale;

  ### Center of Base
  mid.x = (x[1] + x[2]) / 2;
  mid.y = (y[1] + y[2]) / 2;

  ### Step 3
  dc = R^2 - d2/4;
  if(any(dc < 0)) stop("Invalid Lens Radius: too small!");
  dc = sqrt(dc) * sign(R);
  # TODO: Proper sign of d1 and d2
  d1 = dc[1];
  d2 = - dc[2];

  # Step 4: Circle Centers
  slope = slope(x, y);
  c1 = shiftLine(mid.x, mid.y, slope = slope, d = d1);
  c2 = shiftLine(mid.x, mid.y, slope = slope, d = d2);

  # Step 5: Circle Segments
  R = abs(R);
  alpha = asin(d/(2*R));
  phi1 = atan2((mid.y - c1$y), (mid.x - c1$x));
  phi2 = atan2((mid.y - c2$y), (mid.x - c2$x));
  clock = - c(-1,1);
  phi1 = phi1 + clock*alpha[1];
  phi2 = phi2 + clock*alpha[2];
  phi1 = rev(phi1);
  phi2 = rev(phi2);

  # TODO: handle data.frame?
  c1 = unlist(c1[1, c(1, 2)]);
  c2 = unlist(c2[1, c(1, 2)]);
  lst1 = list(r = R[1], center = c1, phi = phi1);
  lst2 = list(r = R[2], center = c2, phi = phi2);
  class(lst1) = c("circle.arc", "list");
  class(lst2) = c("circle.arc", "list");
  # TODO: col and fill

  lst = list(lst1, lst2);
  class(lst) = c("bioshape", "list");
  return(lst);
}

# x, y = endpoints of Axis (group of lens);
#' @export
lens.group = function(x, y, h, pos=NULL, R=NULL, l.scale=1,
                      lwd=1, col=1, fill=NULL, lty=1) {
  len = length(h);
  if(is.null(pos)) {
    if(len > 2) stop("Specify the position of all the lenses!");
    pos = c(0, 1);
  }
  slope = slope(x, y);
  ### Lens Positions
  xL = (1 - pos)*x[1] + pos*x[2];
  yL = (1 - pos)*y[1] + pos*y[2];
  #
  isR = ! is.null(R);
  lst = lapply(seq(len), function(id) {
    d  = h[id]; d = c(-d, d);
    xy = shiftLine(xL[id], yL[id], slope=slope, d=d);
    R  = if(isR) R[id] else NULL;
    scR = if(length(l.scale) == 1) l.scale else l.scale[id];
    lwd = if(length(lwd) == 1) lwd else lwd[id];
    col = if(length(col) == 1) col else col[id];
    lst = lens(xy$x, xy$y, R=R, scale=scR, lwd=lwd, col=col);
    lty = if(length(lty) == 1) lty else lty[id];
    lst[[1]]$lty = lty;
    lst[[2]]$lty = lty;
    if( ! is.null(fill)) {
      fill = if(length(fill) == 1) fill else fill[id];
      lst[[1]]$fill = fill;
      lst[[2]]$fill = fill;
    }
    return(lst);
  });
  lst = do.call(c, lst);
  return(as.bioshape(lst));
}

### Neuron body
#' @export
neuron = function(center = c(0, 0), n = 5, r = 2, phi = 0,
                  axon.length = 3 * r, dendrite.length = ~ r/2, r.nucl = ~ (R - r)/2,
                  r.synapse = r*2/3, type.syn = c("Solid", "Tree", "Detail", "Radial"),
                  col = "red", col.nucl = 1, fill.nucl = NULL, col.synapse = col) {
  body = neuron.body(center = center, n = n, r = r, phi = phi, col = col);
  ### Init
  axon.length = axon.length; # force = scale * r;
  R = r;
  r = attr(body, "r");
  phin = 2 * pi/n;
  phi0 = - phin/2;
  pi20 = pi/2; pi32 = 3*pi/2;
  phiD = seq(n) * phin + phi;
  phiR = phiD - (2*pi) * floor(phiD / (2*pi));
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
    axon = list(x = c(x0[n], xy[1]), y = c(y0[n], xy[2]));
    axon = list(axon);
    lenDendites = n - 1;
    if( ! is.null(type.syn)) {
      syn  = synapse(xy, phi = phiA, type = type.syn, l = r.synapse, fill = col.synapse);
      axon = c(axon, syn);
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
  ### Nucleus
  r.nucl = eval.formula(r.nucl);
  if(r.nucl > 0) {
    nucleus = list(r = r.nucl, center = center, col = col.nucl, fill = fill.nucl);
    class(nucleus) = c("circle", "list");
    nucleus = list(nucleus);
  } else nucleus = NULL;

  ### Neuron
  neuron = c(Body = body, Axon = axon,
             Dendrites = dend, Nucleus = nucleus);
  if( ! inherits(neuron, "bioshape")) {
    class(neuron) = c("bioshape", "list");
  }
  return(neuron);
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
neuron.body = function(center = c(0, 0), n = 5, r = 3, phi = 0, col = 1){
  phi0 = phi + pi/n;
  cc = circlesOnFixedCircle(n = n, r = r, center = center, phi = phi0);
  R = r;
  r = attr(cc, "r");
  phin = 2*pi/n;
  id = seq(n);
  a1 = pi/2 + phi + id * phin;
  a2 = pi/2 + phi + (id + (n - 2)/2) * phin;
  lst = lapply(id, function(id){
    lst = list(
      r = r, center = c(x = cc$x[id], y = cc$y[id]),
      phi = c(a1[id], a2[id]) );
    lst$col = col;
    class(lst) = c("circle.arc", "list");
    return(lst);
  })
  attr(lst, "r") = r;
  class(lst) = c("bioshape", "list");
  return(lst);
}

# phi = angle of axon (in radians);
# slope = slope of axon, as alternative to angle phi;
# alpha = angle of tree-cone (in degrees);
#' @export
synapse = function(p, phi, type = c("Solid", "Tree", "Detail", "Radial"),
                   slope=NULL, col=1, fill="#808080", l=1/2, n=4, alpha = 120,
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
    return(as.bioshape(lst));
  }
  # Tree Types
  if(phi >= 0) {
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
  } else if(type == "Solid") {
    px = l * cos(th) + p[1];
    py = l * sin(th) + p[2];
    xy = list(x = c(p[1], px), y = c(p[2], py), col=col, fill=fill);
    class(xy) = c("polygon", "list");
  } else stop("Not yet implemented!");
  xy = as.bioshape(list(xy));
  return(xy);
}

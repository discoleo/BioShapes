#######################################
#
# BioShapes
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


### Spiral
# p1, p2 = endpoints of spiral;
#  OR if slope specified = starting point + length;
# n = number of turns;
# slope = slope of the inner core;
# N = number of points used to draw the spiral;
#' @export
spirals = function(p1, p2, n=5.5, A=1, phi=0, slope=NULL, N=128) {
  if(is.null(slope)) {
    x = c(p1[1], p2[1]);
    y = c(p1[2], p2[2]);
    slope = slope(x, y);
    l = sqrt((p1[1]- p2[1])^2 + (p1[2]- p2[2])^2);
  } else {
    if(length(p2) > 1) stop("Provide either: length and slope, or the 2 endpoints!")
    l = p2;
  }
  #
  n = 2*n*pi;
  ninv = 1 / n;
  v = l * ninv;
  t = seq(0, n, length.out=N);
  # Rotation matrix: by column
  # rotm = matrix(sdiv * c(1, s, -s, 1), ncol=2, nrow=2);
  if(slope == -Inf || (slope != Inf && p1[2] > p2[1])){
    v = - v;
  } else { phi = phi + pi; }
  x  = v*t;
  y  = A*sin(t + phi);
  xc = x + cos(t + phi);
  #
  if(abs(slope) == Inf) {
    sgn = sign(slope);
    dx  = y;
    dy  = xc;
    lst = list(x = p1[1] + dx, y = p1[2] + dy);
    lst = list(lst);
    class(lst) = c("bioshape", class(lst));
    return(lst);
  }
  sdiv = 1 / sqrt(slope^2 + 1);
  #
  dx = (xc - slope*y) * sdiv; # + p1[1];
  dy = (slope*xc + y) * sdiv; # + p1[2];
  lst = list(x = p1[1] + dx, y = p1[2] + dy);
  lst = list(lst);
  class(lst) = c("bioshape", class(lst));
  return(lst);
}


# n = number of loops;
# N = number of points to draw curve;
# A = amplitude;
# phi = phase shift of sinusoid;
#' @export
helix = function(p1, p2, n = 3, A = 1, phi = 0, parts = 0, N = 128,
                 lwd = NULL, col = NULL, slope = NULL) {
  if(is.null(slope)) {
    x = c(p1[1], p2[1]);
    y = c(p1[2], p2[2]);
    slope = slope(x, y);
    l = sqrt((p1[1] - p2[1])^2 + (p1[2] - p2[2])^2);
  } else {
    if(length(p2) > 1)
      stop("Provide either: length and slope, or the 2 endpoints!");
    l = p2;
  }
  #
  npi = 2 * pi * n;
  ninv = 1 / npi;
  t = seq(0, npi, length.out=N);
  if(x[1] > x[2]) l = -l;
  sinusoid = function(tx, ty) {
    x  = l*tx;
    y  = A*sin(ty + phi);
    #
    if(abs(slope) == Inf) {
      sgn = sign(slope);
      dx = y; dy = x * sgn;
      lst = list(x = p1[1] + dx, y = p1[2] + dy);
      return(lst);
    }
    sdiv = 1 / sqrt(slope^2 + 1);
    # Rotation matrix: by column
    # rotm = matrix(sdiv * c(1, s, -s, 1), ncol=2, nrow=2);
    dx = (x - slope*y) * sdiv; # + p1[1];
    dy = (slope*x + y) * sdiv; # + p1[2];
    lst = list(x = p1[1] + dx, y = p1[2] + dy);
  }
  lst = sinusoid(t*ninv, t);
  if( ! is.null(lwd)) lst$lwd = lwd;
  if( ! is.null(col)) lst$col = col;
  lst = list(lst);
  if(parts > 0) {
    tl = helix.link(n, parts);
    tt = tl$id / tl$div;
    seg = sinusoid(tt, (2*pi)*tl$id);
    attr(lst, "segments") = seg;
  }
  class(lst) = c("bioshape", class(lst));
  return(lst);
}

#' @export
helix.rad = function(R=3, n=8, center=c(0,0), r=1, phi=0,
                     lwd = 1, col = NULL, fill = NULL, N=257) {
  id = seq(0, 1, length.out = N)
  rr = r * cos(2*pi*n*id + phi);
  RR = exp(2i*pi*id);
  RR = R*RR + rr*RR;
  x  = Re(RR) + center[1];
  y  = Im(RR) + center[2];
  xy = list(x=x, y=y);
  if( ! is.null(lwd)) xy$lwd = lwd;
  if( ! is.null(col)) xy$col = col;
  if( ! is.null(fill)) xy$fill = fill;
  class(xy) = c("polygon", "list");
  xy = list(xy);
  class(xy) = "bioshape";
  return(xy);
}

#' @export
helix.link = function(n, k=3, phi=pi/2) {
  if(n < 0) stop("Only positive n!");
  iN = floor(n);
  nTail = n - iN;
  len = (2*k + 2)*iN + 1;
  if(nTail >= 0.5) {
    len = len + k + 1;
    iN  = iN + 0.5;
    nTail = nTail - 0.5;
  }
  id = seq(0, iN, length.out = len);
  # TODO: sin(x + phi) = sin(x);
  # [incorrect] works only with integer k!
  # id = id[ - seq(0, len, by = k+1)];
  if(nTail > 0) {
    # TODO
  }
  #
  return(list(id=id, div = tail(id, 1)));
}


### DNA
# S1 = Options for Segment 1;
#  - L1: draw only 1 line if shorter than this;
#  - L2: draw 2 lines if length between L2 & L1;
#  - short: use enhanced algorithm for length < pi - short * d_theta;
#' @export
dna.new = function(x, y, n=3, phi=c(pi/2, pi) + pi/4, A=1, n.lines = 6,
			lwd=1, lwd.lines = lwd,
			col = c("red", "green"), col.lines = col,
			S1 = list(L1 = 1.5, L2 = 2, short = 0.5, LE2 = 2)) {
  phi = as.radians0(phi);
  p1 = c(x[1], y[1]); p2 = c(x[2], y[2]);
  h1 = helix(p1, p2, n=n, A=A, phi=phi[1], lwd=lwd, parts=0);
  h2 = helix(p1, p2, n=n, A=A, phi=phi[2], lwd=lwd, parts=0);
  # TODO: use helix();
  if( ! is.null(col)) {
    h1[[1]]$col = col[1]; h2[[1]]$col = col[2];
  }
  lst = c(Helix1 = h1, Helix2 = h2);
  if(n.lines == 0) return(as.bioshape(lst));
  #
  if(is.helix.rev(phi)) {
    col.lines = rev(col.lines);
  }
  # TODO: verify & correct;
  if(x[1] == x[2]) col.lines = rev(col.lines);
  pp  = which.intersect.sin(phi, n);
  len = length(pp$x0);
  pi2 = 2*pi*n;
  p0  = c(0, pp$x0 / pi2);
  Ln  = dist.xy(x, y);
  slope = slope(x, y);
  lenL1 = pi / (n.lines + 1);
  phix0 = pp$x0[[1]]; # first Intersection/Segment
  hasL1 = phix0 >= lenL1;
  lstLL = list();
  # print(cbind(phi, c(pp$x0[[1]], pp$x1[[1]])));
  as.lines.dna = function(pp) {
    ### pp = pp[- c(1, 6)];
    py = A * sin(pi2*pp + phi[1]);
    h1 = rotate(pp * Ln, py, slope=slope, p1);
    py = A * sin(pi2*pp + phi[2]);
    h2 = rotate(pp * Ln, py, slope=slope, p1);
    lenL = length(h1$x);
    tmp = lapply(seq(2, lenL - 1), function(iL) {
      data.frame(
        x = c(h1$x[iL], h2$x[iL]),
        y = c(h1$y[iL], h2$y[iL]), id = iL);
    });
    tmp = do.call(rbind, tmp);
  }
  for(i in seq(len)) {
    pS = p0[i];
    pE = p0[i + 1];
	# Starting-Segment:
	if(i == 1) {
		if( ! hasL1) next;
		lenS1 = round(phix0 / lenL1) + 2;
		if(phix0 < pi - S1$short * lenL1) {
			# very short segment;
			if(phix0 <= S1$L2*lenL1) lenS1 = ifelse(phix0 <= S1$L1*lenL1, 2, 3);
			# include 1st line close to start of segment;
			pp = seq(pS + lenL1/(4*pi2), pE, length.out = lenS1);
			pp = c(pS, pp);
		} else pp = seq(pS, pE, length.out = lenS1);
	} else {
		pp = seq(pS, pE, length.out = n.lines + 2);
	}
	tmp  = as.lines.dna(pp);
    colL = if(i %% 2 == 1) col.lines[1] else col.lines[2];
    tmp$lwd = lwd.lines;
    tmp$col = colL;
    lstLL = c(lstLL, list(tmp));
  }
  ### DNA-End
  # TODO: advanced-processing of DNA-end;
  pS  = p0[length(p0)]; pE = 1;
  dnE = pi2 * (pE - pS) / lenL1;
  # print(c(pS, dnE));
  if(dnE >= S1$LE2) {
    lenS1 = round(dnE) + 2;
    pp   = seq(pS, pE, length.out = lenS1);
    tmp  = as.lines.dna(pp);
    colL = if(len %% 2 == 1) col.lines[2] else col.lines[1];
    tmp$lwd = lwd.lines;
    tmp$col = colL;
    lstLL = c(lstLL, list(tmp));
  }
  #
  lst = c(Lines = lstLL, lst);
  return(as.bioshape(lst));
}

### Genome types
#' @export
genome = function(r, center = c(0,0), phi.dna = pi/4,
                  type = c("DNA", "helix", "arc", "circle", "cds", "css", "line", "none"),
                  phi.arc = pi/6,
                  lwd=2, col = "#64D000", ...) {
  type = match.arg(type);
  if(type == "none") return(NULL);
  as.xy = function() {
    p1x = r * cos(phi.dna) + center[1];
    p1y = r * sin(phi.dna) + center[2];
    # p1  = c(p1x, p1y);
    p2x = r * cos(phi.dna + pi) + center[1];
    p2y = r * sin(phi.dna + pi) + center[2];
    # p2  = c(p2x, p2y);
    return(list(x = c(p2x, p1x), y = c(p2y, p1y)));
  }
  as.pp = function(xy) list(p1 = c(xy$x[1], xy$y[1]), p2 = c(xy$x[2], xy$y[2]));
  #
  if(type == "DNA") {
    p = as.xy(); # TODO: col = col with 2 cols;
    lst = dna.new(p$x, p$y, lwd=lwd, ...);
    return(invisible(lst));
  }
  if(type == "helix") {
    p = as.xy();
    p = as.pp(p); # TODO: standardize;
    lst = helix(p$p1, p$p2, lwd=lwd, col=col, ...);
    return(invisible(lst));
  }
  if(type == "arc") {
    if(length(phi.arc) == 1) phi.arc = c(phi.arc, 2*pi - phi.arc);
    lst = list(r = r, center = center, phi = phi.arc, lwd=lwd, col=col);
    class(lst) = c("circle.arc", "list");
    return(as.bioshape(list(lst)));
  }
  if(type == "circle") {
    lst = list(r = r, center = center, lwd=lwd, col=col);
    class(lst) = c("circle", "list");
    return(as.bioshape(list(lst)));
  }
  if(type == "line") {
    p = as.xy();
    p$lwd = lwd; p$col = col;
    return(as.bioshape(list(p)));
  }
  if(type == "css") {
    lst = helix.rad(R = r, center = center, n=10, phi=0, r = r/8,
                    lwd=lwd, col=col);
    return(as.bioshape(list(lst)));
  }
  if(type == "cds") {
    l1 = helix.rad(R = r, center = center, n=10, phi=0, r = r/8,
                   lwd=lwd[1], col=col[1]);
    if(length(lwd) == 1) lwd = rep(lwd, 2);
    if(length(col) == 1) col = rep(col, 2);
    l2 = helix.rad(R = r, center = center, n=10, phi=pi/2, r = r/8,
                   lwd=lwd[2], col=col[2]);
    lst = list(H1 = l1[[1]], H2 = l2[[1]]);
    return(as.bioshape(lst));
  }
}

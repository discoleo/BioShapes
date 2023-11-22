#######################################
#
# Title: BioShapes
# Maintainer: L. Mada
#
# https://github.com/discoleo/BioShapes
#
# Continuation of:
# 1. BSc Thesis: Adrian Cotoc (2022-2023)
# Faculty of Mathematics and Informatics, UVT
#
# Coordinator:
#   Prof. Daniela Zaharie
#   Dr. med. Leonard Mada (Syonic SRL)
#   in collaboration with Syonic SRL
#   [old] GitHub: https://github.com/Adi131313/BioShapes
#
# 2. BSc Thesis: Darian Voda (2021-2022)


#### Cell-like Objects ####


### Smooth muscle cell or Fibroblast
#' @export
cell.SmoothMuscle = function(x, y, r = 1, r.nc = r/2.5, slope = NULL, t.nc = c(0.5, 0.5),
		lwd = 1, col = NULL, fill = NULL, col.nc = 1, fill.nc = NULL, N=128, phi = c(0, pi)) {
	if(is.null(slope)) slope = slope(x, y);
	d = sqrt((x[1] - x[2])^2 + (y[1] - y[2])^2);
	pL = seq(0, d, length.out = N);
	pp = shift.point(c(x[1], y[1]), d=pL, slope=slope);
	pL = seq(x[1], x[2], length.out = length(pL));
	px = pp[,1] - x[1]; px = px * pi / max(abs(px));
	# Margin 1:
	pS = r * sin(px + phi[1]) + pp[,2];
	lst = list(x = pL, y = pS);
	# Margin 2:
	pS = r * sin(px + phi[2]) + pp[,2];
	lst$x = c(lst$x, rev(pL));
	lst$y = c(lst$y, rev(pS));
	#
	lst$lwd = lwd;
	if( ! is.null(col)) lst$col = col;
	if( ! is.null(fill)) lst$fill = fill;
	class(lst) = c("polygon", "list");
	# Nucleus:
	if(r.nc != 0) {
		# TODO: ellipse;
		mid.x = (1-t.nc[1])*x[1] + t.nc[1]*x[2];
		mid.y = (1-t.nc[2])*y[1] + t.nc[2]*y[2];
		center = c(mid.x, mid.y);
		nc = list(r = r.nc, center = center, col = col.nc);
		if( ! is.null(fill.nc)) nc$fill = fill.nc;
		class(nc) = c("circle", "list");
		lst = list(Cell = lst, N = nc);
	} else lst = list(lst);
	lst = as.bioshape(lst);
	return(lst)
}


### Brush-Border Cells
# p1 = Base of Cell (Point 1);
# w, h = width, height;
# n = number of half-cycles of sine-wave;
#' @export
as.options.BrushBorder = function(...) {
	opt = list(...);
	len = length(opt);
	if(len == 1 && is.list(opt[[1]])) opt = opt[[1]];
	if(is.null(opt$n)) opt$n = 6.5;
	if(is.null(opt$A)) opt$A = 0.5;
	if(is.null(opt$N)) opt$N = 128;
	if(is.null(opt$phi)) opt$phi = 0;
	return(opt);
}
#' @export
cell.BrushBorder = function(p1, w, h, slope=0, r.nc = ~ 1/5, t.nc = c(1/2, 7/20),
		lwd = 1, lwd.nc = lwd, col = NULL, fill = NULL, col.nc = 1, fill.nc = NULL,
		brush.options = list(n=6.5, A=0.5, N=128, phi=0)) {
	# Cell:
	p11 = p1;
	if(is.null(slope)) {
		if(length(w) < 2) stop("Provide slope or the 2 base points!");
		x = c(p1[1], w[1]); y = c(p1[2], w[2]);
		p12 = w; w = dist.xy(x, y);
		slope = slope(x, y);
	} else {
		if(length(w) > 1) warning("Wrong width!");
		p12 = shiftPoint(p1, d=w, slope=slope);
	}
	p21 = shift.ortho(p1, d=h, slope=slope);
	p21 = unlist(p21);
	p22 = shiftPoint(p21, d=w, slope=slope);
	# Brush-Border:
	opt   = as.options.BrushBorder(brush.options);
	brush = helix(p21, p22, n = opt$n, A = opt$A, phi = opt$phi, N = opt$N);
	brush = brush[[1]];
	brush$x = c(brush$x, p22[1], p12[1], p11[1], p21[1]);
	brush$y = c(brush$y, p22[2], p12[2], p11[2], p21[2]);
	brush$lwd = lwd;
	if( ! is.null(col)) brush$col = col;
	if( ! is.null(fill)) brush$fill = fill;
	class(brush) = c("polygon", "list");
	lst = list(Cell = brush);
	# Nucleus:
	isFormula = inherits(r.nc, "formula");
	if(isFormula || r.nc != 0) {
		t.nc[2] = 1 - t.nc[2];
		center = center.p4(p11, p12, p21, p22, t = t.nc);
		if(isFormula) {
			d = sqrt((p11[1] - p12[1])^2 + (p11[2] - p12[2])^2);
			r = eval(r.nc[[2]]) * d;
		} else r = r.nc;
		nc = list(r = r, center = center, lwd = lwd.nc, col = col.nc);
		if( ! is.null(fill.nc)) nc$fill = fill.nc;
		class(nc) = c("circle", "list");
		lst$N = nc;
	}
  return(as.bioshape(lst));
}

###
#' @export
epithelium.brush = function(x, y, n=5, h = ~ 2, r.nc = ~ 1/5, t.nc = c(1/2, 7/20),
		lwd = 1, lwd.nc = lwd, col = NULL, fill = NULL, col.nc = 1, fill.nc = NULL,
		brush.options = list(n=6.5, A=0.5, N=128, phi=0)) {
	d0 = dist.xy(x, y);
	d1 = d0 / n;
	h  = if(inherits(h, "formula")) eval(h[[2]]) * d1 else h;
	xy = split.line(x, y, n=n);
	lst.args = list(h=h, r.nc=r.nc, t.nc=t.nc, slope = NULL, lwd=lwd, lwd.nc=lwd.nc,
		col=col, fill=fill, col.nc=col.nc, fill.nc=fill.nc, brush.options=brush.options);
	ep = lapply(seq(1, n), function(id) {
		p1 = c(xy$x[[id]], xy$y[[id]]);
		p2 = c(xy$x[[id + 1]], xy$y[[id + 1]]);
		lst.args$p1 = p1; lst.args$w = p2;
		do.call(cell.BrushBorder, lst.args);
	})
	invisible(as.bioshape(ep));
}


### Epithelium: Mono-Layer
# - cells are randomly deformed;
#' @export
epithelium.rMonoLayer = function(x, y, n = 10, h = 1,
		dx.sc = c(0.125, 0.08), dy.sc = c(0.125, 0.08)) {
	lst = polygons.rLinear(x=x, y=y, n=n, h=h, dx.sc=dx.sc, dy.sc=dy.sc);
	invisible(lst);
}

### Polygons:
# - linearly arranged & randomly deformed;
#' @export
polygons.rLinear = function(x, y, n, h = 1,
		dx.sc = c(0.125, 0.08), dy.sc = c(0.125, 0.08)) {
	L = dist.xy(x, y);
	x0 = seq(0, L, length.out = n+1);
	dx0 = L / n; dx = dx0 * dx.sc[1]; dxm = dx0 * dx.sc[2];
	# Bottom:
	xB = runif(n - 1, -dx, dx);
	xB = c(0, xB, 0);
	xB = x0 + xB;
	# Top:
	xM = runif(2, -dxm, dxm);
	xT = runif(n - 1, -dx, dx);
	xT = c(xM[1], xT, xM[2]);
	xT = x0 + xT;
	# yT:
	dy = abs(h) * dy.sc[1]; dym = abs(h) * dy.sc[2];
	y0 = rep(h, n + 1);
	yT = y0 + runif(n + 1, -dy, dy);
	idL = n + 1;
	xm = (xT[-1] + xT[- idL]) / 2;
	xm = xm + runif(n, -dxm, dxm);
	ym = (yT[-1] + yT[- idL]) / 2;
	ym = ym + runif(n, -dym, dym);
	# Rotation:
	phi = atan2(y[2] - y[1], x[2] - x[1]);
	cs = cos(phi); sn = sin(phi);
	mR = matrix(c(cs, sn, -sn, cs), nrow = 2);
	#
	yB  = rep(0, n + 1);
	xy0 = c(x[1], y[1]);
	xyB = mR %*% rbind(xB, yB) + xy0;
	xyT = mR %*% rbind(xT, yT) + xy0;
	xyM = mR %*% rbind(xm, ym) + xy0;
	xyB[1,1] = x[1]; xyB[1, n+1] = x[2];
	#
	lst = lapply(seq(n), function(id) {
		id1 = id + 1; id2 = c(id, id1);
		x = c(xyB[1, id2], xyT[1, id1], xyM[1, id], xyT[1, id]);
		y = c(xyB[2, id2], xyT[2, id1], xyM[2, id], xyT[2, id]);
		lst = list(x=x, y=y);
		class(lst) = c("polygon", "list");
		return(lst);
	})
	return(as.bioshape(lst));
}


#################

### Muscle tissue
#' @export
muscle = function(scale.x = 1.5, scale.R = c(1.5, 1.5),
                  x = c(-2, 2), y = c(1, 1), dy = 0.4, dx = 2, n = 6, fill = "red"){

  lst = list()
  for(iy in seq(0,n)){
    if(iy%%2 == 0){
      for(idx in seq(-2,2,2)){
        lens = lens(x = x + idx*dx, y = y + iy*dy, scale=scale.R, fill = fill)
        lst = c(lst, lens)
      }}
    else{
      for(idx in seq(-3,3,2)){
        lens = lens(x = x + idx*dx, y = y + iy*dy, scale=scale.R, fill = fill)
        lst = c(lst, lens)
      }
    }}
  lst = as.bioshape(lst)
  return (invisible(lst))
}


###################
### Blood Cells ###

# TODO

### Red blood cells
#' @export
draw_blood_cell = function(center = c(0, 0),
                           radius = 1, lwd = 10, col = "#901000", fill = "red"){
  #lim = c(-radius, radius) * 2;
  #plot.base(xlim = lim, ylim = lim)
  #circle <- seq(0, 2 * pi, length.out = 80)
  #x <- radius * sin(circle)
  #y <- radius * cos(circle)
  #polygon(x, y, col = "red", border = "#901000", lwd = lwd)
  lst = list(center = center, r = radius,
             fill = fill, col = col, lwd = lwd);
  class(lst) = c("circle", "list");
  lst = list(lst);
  class(lst) = c("bioshape", "list");
  return(lst)

}


########################
########################

### Tumors

# n = number of cells;
# r = radius of each cell;
# R = radius of tumor mass;
# r.nc = radius of nucleus;
#' @export
tumor.mass = function(n = 50, r = 0.25, R = 6.5*r, r.nc = r/4,
		center = c(0,0), phi = 0, fill = "#6480FF", fill.nc = "purple",
		top = c("2out", "in", "out"), adj.nc = TRUE, ...) {
	xy = uniform.circle(n, r = R, center=center, phi=phi, ...);
	top = match.arg(top);
	if(top == "2out") {
		id = c(seq(1, n, by=2), seq(2, n, by=2));
		xy = xy[id,];
	} else if(top == "in") xy = xy[seq(n, 1), ];
	# Cells:
	lst = list(r=r, center = xy, fill = fill);
	class(lst) = c("circle", "list;");
	lstAll = list(Cells = lst);
	# Nucleus:
	if( ! is.null(adj.nc)) {
		r.adj = if(is.logical(adj.nc)) r.nc / 4 else r.nc * adj.nc;
		xy[,1] = xy[,1] + runif(n, -r.adj, r.adj);
		xy[,2] = xy[,2] + runif(n, -r.adj, r.adj);
	}
	lst = list(r = r.nc, center = xy, fill = fill.nc);
	class(lst) = c("circle", "list;");
	lstAll$Nuclei = lst;
	lstAll = as.bioshape(lstAll);
	invisible(lstAll);
}

### Mixed Tumor
#' @export
tumor.mix = function(d = c(0.53, 0.49), R = 2, phi = c(0, pi),
		fill = c("#6480FF", "#F09648"), by = 3) {
	c1 = tumor.mass(60, R=R, d = d[1], phi = phi[1], fill = fill[[1]]);
	if(length(d) == 1) {
		c2  = c1;
		len = nrow(c1$Cells$center);
		c1$Cells$center  = c1$Cells$center[- seq(1, len, by=by), ];
		c1$Nuclei$center = c1$Nuclei$center[- seq(1, len, by=by), ];
		c2$Cells$center  = c2$Cells$center[seq(1, len, by=by), ];
		c2$Nuclei$center = c2$Nuclei$center[seq(1, len, by=by), ];
		c2$Cells$fill = fill[[2]];
		c.mix = as.bioshape(list(c1, c2));
	} else {
		c2 = tumor.mass(30, R=R, d = d[2], phi = phi[2], fill = fill[[2]]);
		c.mix = mix.cells(c1, c2, by=by);
	}
	invisible(c.mix);
}
# Mix 2 populations of cells
#' @export
mix.cells = function(x1, x2, by = 3) {
	len = nrow(x1$Cells$center);
	tmp = x1;
	x1$Cells$center = x1$Cells$center[seq(1, len, by = by), ];
	x1$Nuclei$center = x1$Nuclei$center[seq(1, len, by = by), ];
	tmp$Cells$center = tmp$Cells$center[- seq(1, len, by = by), ];
	tmp$Nuclei$center = tmp$Nuclei$center[- seq(1, len, by = by), ];
	lst = list(tmp, x2, x1);
	lst = as.bioshape(lst);
	invisible(lst);
}

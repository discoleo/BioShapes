

#' @export
phage.base = function(x, y = NULL, r = 1, phi = 0,
		col = NULL, fill = NULL, id.faces = FALSE) {
	center = xy.coords(x, y);
	center = c(center$x, center$y);
	phi = phi + pi/6 + pi;
	pp = points.circle(6, center=center, r=r, phi=phi);
	mid = function(id, pp) {
		list(
			x = (pp$x[id[1]] + pp$x[id[2]]) / 2,
			y = (pp$y[id[1]] + pp$y[id[2]]) / 2);
	}
	adj = function(m, px, py, t = 1/2) {
		ti = 1 - t;
		list(
			x = (t*m$x + ti*px),
			y = (t*m$y + ti*py));
	}
	mT = mid(c(1,3), pp);
	mL = mid(c(3,5), pp);
	mR = mid(c(5,1), pp);
	mT = adj(mT, pp$x[2], pp$y[2], t = 2/3)
	mL = adj(mL, pp$x[4], pp$y[4]);
	mL = adj(mL, pp$x[5], pp$y[5], t = 2/3);
	mR = adj(mR, pp$x[6], pp$y[6]);
	mR = adj(mR, pp$x[5], pp$y[5], t = 2/3);
	mm = list(x = c(mT$x, mL$x, mR$x), y = c(mT$y, mL$y, mR$y));
	as.t3 = function(idp, idm) {
		x = c(pp$x[idp], mm$x[idm]);
		y = c(pp$y[idp], mm$y[idm]);
		lst = list(x=x, y=y);
		class(lst) = c("polygon", "list");
		return(lst);
	}
	lst = list(
		as.t3(c(1,2), 1), as.t3(c(2,3), 1),
		as.t3(c(3,4), 2), as.t3(c(4,5), 2),
		as.t3(c(5,6), 3), as.t3(c(6,1), 3),
		as.t3(1, c(1,3)), as.t3(3, c(1,2)), as.t3(5, c(3,2)),
		as.t3(NULL, c(1,2,3))
		);
	if(! is.null(fill)) {
		len = length(fill);
		if(len > 1) {
			for(id in seq(len)) lst[[id]]$fill = fill[[id]];
		} else lst$fill = fill;
	}
	if(! is.null(col)) {
		len = length(col);
		if(len > 1) {
			for(id in seq(len)) lst[[id]]$col = col[[id]];
		} else lst$col = col;
	}
	# Number faces
	if(id.faces) {
		len = 10;
		cc = lapply(seq(len), function(id) {
			p = lst[[id]];
			center.xy(p$x, p$y);
			});
		cc = do.call(rbind, cc);
		lstID = list(x = cc[,1], y = cc[,2],
			labels = seq(len));
		class(lstID) = c("text", "list");
		lst$Count = lstID;
	}
	return(as.bioshape(lst));
}


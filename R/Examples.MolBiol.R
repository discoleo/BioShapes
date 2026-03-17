

### Western Blot

### Empty Gel
#' @export
example.gel = function(h = 8, h.well = 1.5, n = 8, x = c(0, 8), y = c(0,0),
		fill = c("#C2C2C2", "#3224E8")) {
	tmp = gel(x, y, n=n, h=h, h.well = h.well, fill = fill)
	plot.base(axt = NULL);
	lines(tmp);
	invisible(tmp);
}

# Slightly Skewed Gel:
# - Student was not careful enough when placing the gel
#   in the imaging device;
# - Also forgot to run the electrophoresis;
# example.gel(x = c(0, 8), y = c(0, 1), n = c(8, 6))


### Western Blot

### Empty Gel
#' @export
example.gel = function(h = 8, h.well = 1.5, x = c(0, 8), y = c(0,0),
		fill = "#C2C2C2") {
	tmp = gel(x, y, h=h, h.well = h.well, fill = fill)
	plot.base(axt = NULL);
	lines(tmp);
	invisible(tmp);
}

# Slightly Skewed Gel:
# - Student was not careful enough when placing the gel
#   in the imaging device;
# example.gel(x = c(0, 8), y = c(0, 1))
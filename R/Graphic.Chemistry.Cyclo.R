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

#' @export
circle.ngon = function(R=3, center=c(0,0), r=1/4, n=3, phi=0, clockwise=TRUE,
                       r.adj = -0.01, phi.adj=0) {
  if(n %% 2 == 0) { d = 2; }
  else d = (1 + cos(pi/n));
  # ng = number of ngons;
  th = acos(1 - (d*r/R)^2);
  ng = round(2*pi / th);
  shift = function(x) (x + c(x[-1], x[1])) / 2;
  th = 2*pi*seq(0, ng) / ng + phi;
  px = R * cos(th) + center[1];
  py = R * sin(th) + center[2];
  # Note: contact points are on the circle R;
  cx = shift(px);
  cy = shift(py);
  # corrected r:
  r  = sqrt((px[1] - cx[1])^2 + (py[1] - cy[1])^2) + r.adj;
  r  = 2*r / d;
  t0 = if(clockwise) 3*pi/2 else pi/2;
  th = t0 + th + pi/ng + phi.adj;
  pg = lapply(seq(ng), function(id) {
    ngon(n, r=r, center=c(cx[id], cy[id]), phi = th[id]);
  })
  attr(pg, "r") = r;
  class(pg) = "bioshape";
  return(pg);
}

#' @export
circle.spiro = function(n, R=3, center = c(0,0), ngon=6, phi=0,
                        type = c("auto", "in", "out", "clockwise", "anti-clockwise", "real-clock"),
                        r.adj=0, phi.adj=0) {
  # Note:
  # - Center of ngon is on circle,
  #   unlike "real-clock" where contact points are on circle;
  type = match.arg(type);
  shift = function(x) (x + c(x[-1], x[1])) / 2;
  th = 2*pi*seq(0, n) / n + phi;
  # Center is on the Circle R;
  cx = R * cos(th) + center[1];
  cy = R * sin(th) + center[2];
  # Compute r:
  isEven = (ngon %% 2 == 0);
  if(type == "auto") {
    # Default:
    #  Odd => In; Even => Clockwise;
    type = if(isEven) "clockwise" else "in";
  } else if(isEven && (type == "anti-clockwise" || type == "real-clock")) {
    type = "clockwise"; # the same;
  }
  if(isEven) {
    r = R * sin(pi/n);
  } else if(type == "in" || type == "out") {
    alfa = if(type == "in") 2*pi/ngon else 3*pi/ngon;
    alfa = pi - pi/n - alfa;
    # Note: overlaps possible for ngon > 5;
    r = R * sin(pi/n) / sin(alfa);
  } else {
    # Odd ngon: both clockwise & anti-clockwise;
    cs = cos(pi/ngon);
    r  = 2*R*sin(pi/n) / sqrt(1 + cs^2 - 2*cs*cos(pi*(1 - 2/n)));
    if(type == "real-clock") {
      return(circle.ngon(R=R, r=r, center=center, n = ngon,
                         phi=phi, r.adj=r.adj, phi.adj=phi.adj));
    }
    cs = r * (1 - cs) / (2*R); # shift
    sh = if(type == "clockwise") - acos(cs) else (pi - acos(cs));
    phi.adj = phi.adj + sh;
  }
  r = r + r.adj;
  if(type == "clockwise") {
    t0 = if(isEven) pi/2 else 3*pi/2 + pi/n;
  } else pi/2;
  if(type == "anti-clockwise") t0 = t0 + pi/n;
  if(type == "in") { t0 = pi;} else if(type == "out") t0 = 0;
  th = t0 + th + phi.adj;
  pg = lapply(seq(n), function(id) {
    ngon(ngon, r=r, center=c(cx[id], cy[id]), phi = th[id]);
  })
  attr(pg, "r") = r;
  class(pg) = "bioshape";
  return(pg);
}

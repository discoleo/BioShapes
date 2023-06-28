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
examples.arrows = function(new.plot = TRUE, lwd = 2) {
  ### Plot
  if(new.plot) plot.base();

  ### Row 1:
  y = c(7, 9); yt = 6.7;
  abline(h = y[2], col="green")
  text(4, y[2] + 1,
       "Arrows with different ArrowHeads")

  # Simple ArrowHead
  x = c(-1, -1);
  d = -0.5;
  d.head = c(-0.5,0.5)
  a1 = arrowSimple(x, y, d=d, d.head=d.head, lwd=lwd);
  text(-1, yt,
       "Simple", cex = 0.75);

  # Double ArrowHead
  x = c(1, 1);
  d = -0.3;
  d.head=-0.5
  a2 = arrowDouble(x, y, d=d, d.head=d.head, lwd=lwd);
  text(1, yt,
       "Double", cex = 0.75);

  # Inverted ArrowHead
  x = c(3, 3);
  d = -0.5;
  d.head=c(-0.5,0.5)
  a3 = arrowInverted(x, y, d=d, d.head=d.head, lwd=lwd);
  text(3, yt,
       "Inverted", cex = 0.75);

  # Diamond ArrowHead
  x = c(5, 5);
  d.head = c(-0.5, 0.5);
  d = -0.5;
  arrowDiamond(x, y, d=d, d.head=d.head, lwd=lwd, join=0);
  text(5, yt,
       "Diamond", cex = 0.75);

  # T Shape ArrowHead
  x = c(7, 7);
  arrowT(x, y, d=-0.75, lwd=lwd);
  text(7, yt,
       "T Shape", cex = 0.75);

  # Measurement ArrowHead
  x = c(9.5, 9.5);
  arrowMeasure(x, y, d=-0.5, lwd=lwd);
  text(9.5, yt,
       "Measurement", cex = 0.70);


  ### Row 2
  y = c(3, 5); yt = 2.7;
  abline(h = y[2], col="green")

  # X Shape ArrowHead
  x = c(-1, -1);
  arrowX(x, y, d=0.5, lwd=lwd);
  text(-1, yt,
       "X Shape", cex = 0.70);

  # Square Shape ArrowHead
  x = c(1, 1);
  arrowSquare(x, y, d=-0.5, lwd=lwd);
  text(1, yt,
       "Square Shape", cex = 0.70);

  # Square Flag
  x = c(3, 3);
  arrowSquare(x, y, d=-0.5, d.head=c(0, 2*d), lwd=lwd);
  text(3.5, yt,
       "Square Flag", cex = 0.70);

  # Multiple-Lined ArrowHead
  n = 3; d = 0.5;
  x = c(6, 6);
  arrowN(x, y, n=n, d=d, lwd=lwd);
  text(6, yt,
       "Multiple-Lined", cex = 0.70);

  # Double Lined Inverted ArrowHead
  x = c(9, 9);
  arrowDoubleInverted(x, y, d=-0.3, lwd=lwd);
  text(9, yt,
       "Double-Lined Inverted", cex = 0.70);


  ### Row 3
  y = c(-1, 1); yt = -1.3;
  abline(h = y[2], col="green")

  # Solid Rectangle ArrowHead
  x = c(-1, -1);
  arrowSolidSquare(x, y, d=-0.5, lwd=lwd, col="darkred", fill="red");
  text(-1, yt,
       "Solid Rectangle", cex = 0.70);

  # Triangle ArrowHead
  x = c(1, 1);
  d = -0.5;
  a1 = arrowTriangle(x, y, d=d, lwd=lwd);
  text(1, yt,
       "Triangle", cex = 0.70);

  # Solid Circle ArrowHead
  x = c(3, 3);
  arrowCircle(x, y, r=0.5, lwd=lwd, fill="#FFB0A0");
  text(3, yt,
       "Solid Circle", cex = 0.70);

  # Simple Circle ArrowHead
  x = c(5, 5);
  arrowCircle(x, y, r=0.5, lwd=lwd);
  text(5, yt,
       "Simple Circle", cex = 0.70);
}

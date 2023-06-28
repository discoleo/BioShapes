### Type = "in"
x = c(1,7); y = c(1, 5);
plot.base()
measure(x, y)
points(x, y, col="green")

### Type = "out"
x = c(1,7); y = c(1, 5);
plot.base()
measure(x, y, type="out")
points(x, y, col="green")

### Asymmetric
x = c(1,7); y = c(1, 5);
plot.base()
measure(x, y, type="out", d=c(-0.5, 1.5))
points(x, y, col="green")

### Asymmetric & ---
x = c(1,7); y = c(1, 5);
plot.base()
measure(x, y, type="out", d=c(-0.5, 1.5), lty=2)
points(x, y, col="green")


### Start point on head & ---
x = c(1,7); y = c(1, 5);
plot.base()
measure(x, y, type="out", d=1.5, lty=2)
points(x, y, col="green")

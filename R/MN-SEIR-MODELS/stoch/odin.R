library(odin)
library(dde)
model = odin("MN-SEIR-MODELS/try_stoch.R")
x = model()

seirds_col <- c("#8c8cd9", "#e67300", "#d279a6", "#ff4d4d", "#999966", "#660000")

set.seed(1)
x_res <- x$run(0:365)
par(mar = c(4.1, 5.1, 0.5, 0.5), las = 1)
matplot(x_res[, 1], x_res[, -1], xlab = "Time", ylab = "Number of individuals",
        type = "l", col = seirds_col, lty = 1)
legend("left", lwd = 1, col = seirds_col, legend = c("S", "E", "Ir", "Id", "R", "D"), bty = "n")

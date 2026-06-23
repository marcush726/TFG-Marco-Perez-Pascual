#MODELO SIS
library(ggplot2)

N     <- 20
beta  <- 0.5
gamma <- 0.2
R0    <- beta / gamma
cat(sprintf("N = %d,  beta = %.2f,  gamma = %.2f,  R0 = beta/gamma = %.2f\n",
            N, beta, gamma, R0))

lambda <- function(i) beta * i * (N - i) / N
mu     <- function(i) gamma * i
lam  <- sapply(1:N, lambda)
mu_v <- sapply(1:N, mu)


QT <- matrix(0, N, N)
for (i in 1:N) {
  if (i > 1) QT[i, i-1] <- mu_v[i]
  if (i < N) QT[i, i+1] <- lam[i]
  QT[i, i]  <- -(mu_v[i] + lam[i])
}


ev    <- eigen(t(QT))
idx   <- which.max(Re(ev$values))
alpha <- -Re(ev$values[idx])
v_raw <- Re(ev$vectors[, idx])
if (sum(v_raw) < 0) v_raw <- -v_raw
v_raw[v_raw < 0] <- 0
qsd <- v_raw / sum(v_raw)

cat(sprintf("\nTasa de decaimiento: alpha = %.6f\n", alpha))
cat(sprintf("Verificacion sum(QSD) = %.10f\n", sum(qsd)))
residuo <- as.vector(qsd %*% QT) + alpha * qsd
cat(sprintf("Residuo ||u.QT + alpha.u||_inf = %.2e\n", max(abs(residuo))))

ord    <- order(Re(ev$values), decreasing = TRUE)
alpha2 <- -Re(ev$values[ord[2]])
Ru     <- 2 * alpha / alpha2
cat(sprintf("alpha2 = %.6f | Ru = 2*alpha/alpha2 = %.6f  (%s)\n",
            alpha2, Ru, ifelse(Ru < 1, "QSD adecuada", "QSD inadecuada")))
cat(sprintf("E_qsd[T] = 1/alpha = %.2f unidades de tiempo\n", 1/alpha))


pi_coef    <- numeric(N); pi_coef[1] <- 1
for (i in 2:N) pi_coef[i] <- pi_coef[i-1] * lam[i-1] / mu_v[i]
p0 <- pi_coef / sum(pi_coef)

tau_coef    <- numeric(N); tau_coef[1] <- 1
for (i in 2:N) tau_coef[i] <- tau_coef[i-1] * lam[i-1] / mu_v[i-1]
p1 <- tau_coef / sum(tau_coef)

F_p0  <- cumsum(p0)
F_qsd <- cumsum(qsd)
F_p1  <- cumsum(p1)
cat(sprintf("\nOrden estocastico p^(0) <=_st QSD: %s\n", all(F_p0 >= F_qsd - 1e-10)))
cat(sprintf("Orden estocastico QSD  <=_st p^(1): %s\n", all(F_qsd >= F_p1 - 1e-10)))

pi_ini    <- rep(1/N, N)
invQT_neg <- solve(-QT)
num <- as.vector(matrix(pi_ini, nrow = 1) %*% invQT_neg)
RE  <- num / sum(num)

cat("\nTabla comparativa (estados 1 a 10):\n")
print(data.frame(
  estado = 1:10,
  QSD    = round(qsd[1:10], 6),
  RE     = round(RE[1:10], 6),
  p0     = round(p0[1:10], 6),
  p1     = round(p1[1:10], 6)
))
cat(sprintf("\nMedia E[I] bajo QSD:   %.4f\n", sum((1:N) * qsd)))
cat(sprintf("Media E[I] bajo RE:    %.4f\n", sum((1:N) * RE)))
cat(sprintf("Media E[I] bajo p^(0): %.4f\n", sum((1:N) * p0)))
cat(sprintf("Media E[I] bajo p^(1): %.4f\n", sum((1:N) * p1)))

#1
df <- data.frame(
  estado       = rep(1:N, 4),
  probabilidad = c(qsd, RE, p0, p1),
  distribucion = factor(rep(c("QSD (exacta)", "RE", "p0", "p1"), each = N),
                        levels = c("QSD (exacta)", "RE", "p0", "p1"))
)

ggplot(df, aes(x = estado, y = probabilidad,
               color = distribucion, shape = distribucion)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.2) +
  scale_color_manual(
    name   = NULL,
    values = c("#E41A1C", "#377EB8", "#4DAF4A", "#FF7F00"),
    labels = c("QSD (exacta)", "RE", expression(p^{(0)}), expression(p^{(1)}))
  ) +
  scale_shape_manual(
    name   = NULL,
    values = c(16, 17, 15, 18),
    labels = c("QSD (exacta)", "RE", expression(p^{(0)}), expression(p^{(1)}))
  ) +
  guides(
    color = guide_legend(override.aes = list(shape = c(16, 17, 15, 18))),
    shape = "none"
  ) +
  labs(
    title    = "Modelo SIS: Comparacion de distribuciones",
    x = "Numero de infectados (i)",
    y = "Probabilidad"
  ) +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(),
        legend.text = element_text(size = 14),
        legend.key.size = unit(1.2, "cm"))
ggsave("fig_sis_01_distribuciones.pdf", width = 7, height = 6)

#2
df_cdf <- data.frame(
  estado       = rep(1:N, 4),
  acumulada    = c(F_p0, F_qsd, F_p1, cumsum(RE)),
  distribucion = factor(rep(c("p0", "QSD (exacta)", "p1", "RE"), each = N),
                        levels = c("p0", "QSD (exacta)", "RE", "p1"))
)

ggplot(df_cdf, aes(x = estado, y = acumulada,
                   color = distribucion, linetype = distribucion)) +
  geom_step(linewidth = 0.9) +
  scale_color_manual(
    values = c("#4DAF4A", "#E41A1C", "#377EB8", "#FF7F00"),
    labels = c(expression(p^{(0)}), "QSD (exacta)", "RE", expression(p^{(1)})),
    name   = NULL
  ) +
  scale_linetype_manual(
    values = c("dashed", "solid", "dotted", "dashed"),
    labels = c(expression(p^{(0)}), "QSD (exacta)", "RE", expression(p^{(1)})),
    name   = NULL
  ) +
  guides(
    color    = guide_legend(override.aes = list(
      linetype = c("dashed", "solid", "dotted", "dashed"))),
    linetype = "none"
  ) +
  labs(
    title    = "Modelo SIS: Orden estocastico entre distribuciones",
    x        = "Numero de infectados (k)",
    y        = "F. de distribucion F(k)"
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position  = "bottom",
    legend.text      = element_text(size = 14),
    legend.key.size  = unit(1.2, "cm"),
    panel.grid.minor = element_blank()
  )
ggsave("fig_sis_02_fda.pdf", width = 7, height = 6)

#3
set.seed(42)
x0     <- sample(1:N, size = 1, prob = qsd)
t      <- 0; x <- x0
tiempo <- t; estado <- x

while (x > 0 && t < 150) {
  total <- mu_v[x] + lam[x]
  dt    <- rexp(1, rate = total)
  t     <- t + dt
  x     <- sample(c(x - 1, x + 1), 1, prob = c(mu_v[x], lam[x]) / total)
  tiempo <- c(tiempo, t)
  estado <- c(estado, x)
}

df_traj   <- data.frame(tiempo = tiempo, estado = estado)
media_qsd <- sum((1:N) * qsd)
q25_qsd   <- min(which(cumsum(qsd) >= 0.25))
q75_qsd   <- min(which(cumsum(qsd) >= 0.75))

ggplot(df_traj, aes(x = tiempo, y = estado)) +
  geom_step(color = "#2166AC", linewidth = 0.7) +
  geom_hline(yintercept = media_qsd, linetype = "dashed",
             color = "#E41A1C", linewidth = 0.7) +
  geom_hline(yintercept = c(q25_qsd, q75_qsd), linetype = "dotted",
             color = "#E41A1C", linewidth = 0.6) +
  labs(
    title    = "Modelo SIS - Trayectoria",
    x = "Tiempo",
    y = "Numero de infectados I(t)"
  ) +
  theme_bw(base_size = 13) +
  theme(panel.grid.minor = element_blank())
ggsave("fig_sis_03_trayectoria.pdf", width = 7, height = 6)

#4
set.seed(143)
N_traj <- 600
t_max  <- 150

trajs <- vector("list", N_traj)
for (k in 1:N_traj) {
  x0_k <- 3
  t <- 0; x <- x0_k
  tiempo <- t; estado <- x
  while (x > 0 && t < t_max) {
    total <- mu_v[x] + lam[x]
    dt    <- rexp(1, rate = total)
    t     <- t + dt
    x     <- sample(c(x - 1, x + 1), 1, prob = c(mu_v[x], lam[x]) / total)
    tiempo <- c(tiempo, t)
    estado <- c(estado, x)
  }
  trajs[[k]] <- data.frame(tiempo = tiempo, estado = estado, traj = k)
}

t_eval <- c(5, 20, 50, 100, 1000)
conv_list <- lapply(t_eval, function(tt) {
  estados_vivos <- c()
  for (tr in trajs) {
    idx <- max(which(tr$tiempo <= tt))
    if (tr$estado[idx] > 0) estados_vivos <- c(estados_vivos, tr$estado[idx])
  }
  if (length(estados_vivos) == 0) return(NULL)
  tbl <- table(factor(estados_vivos, levels = 1:N))
  data.frame(
    estado = 1:N,
    prob   = as.numeric(tbl) / sum(tbl),
    t      = tt,
    label  = paste0("t = ", tt, "  (", length(estados_vivos), " vivas)")
  )
})
conv_list <- Filter(Negate(is.null), conv_list)
df_conv   <- do.call(rbind, conv_list)

df_qsd_ref <- data.frame(estado = 1:N, prob = qsd, t = NA, label = "QSD")
etiquetas_conv <- sapply(conv_list, function(d) unique(d$label))

df_conv_plot <- rbind(df_conv, df_qsd_ref)
df_conv_plot$label <- factor(df_conv_plot$label,
                             levels = c(etiquetas_conv, "QSD"))

n_labels <- length(t_eval)
colores  <- c(colorRampPalette(c("#BDD7EE", "#1F497D"))(n_labels), "#E41A1C")
tipos    <- c(rep("solid", n_labels), "dashed")

ggplot(df_conv_plot, aes(x = estado, y = prob, color = label, linetype = label)) +
  geom_line(linewidth = 0.9, na.rm = TRUE) +
  geom_point(size = 1.8, na.rm = TRUE) +
  scale_color_manual(values = colores, name = NULL) +
  scale_linetype_manual(values = tipos, name = NULL) +
  labs(
    title    = "Modelo SIS - Convergencia a la QSD",
    subtitle = bquote(.(N_traj) ~ "trayectorias iniciadas aleatoriamente"),
    x        = "Numero de infectados (j)",
    y        = "Probabilidad condicional"
  ) +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())
ggsave("fig_sis_04_convergencia.pdf", width = 7, height = 6)


NT   <- 50
gamT <- 1.0
R0s  <- c(0.5, 0.9, 1.0, 1.3, 1.5, 2.0)

build_QT <- function(N, beta, gamma) {
  lam <- beta * (1:N) * (N - (1:N)) / N
  muv <- gamma * (1:N)
  Q <- matrix(0, N, N)
  for (i in 1:N) {
    if (i > 1) Q[i, i-1] <- muv[i]
    if (i < N) Q[i, i+1] <- lam[i]
    Q[i, i] <- -(muv[i] + lam[i])
  }
  Q
}
get_QSD <- function(Q) {
  e <- eigen(t(Q)); o <- order(Re(e$values), decreasing = TRUE)
  a  <- -Re(e$values[o[1]]); a2 <- -Re(e$values[o[2]])
  v <- Re(e$vectors[, o[1]]); if (sum(v) < 0) v <- -v; v[v < 0] <- 0
  list(qsd = v / sum(v), alpha = a, alpha2 = a2, Ru = 2 * a / a2)
}
get_RE <- function(Q) { M <- -solve(Q); M / rowSums(M) }
get_p0 <- function(N, beta, gamma) {
  lam <- beta * (1:N) * (N - (1:N)) / N; muv <- gamma * (1:N)
  p <- numeric(N); p[1] <- 1
  for (i in 2:N) p[i] <- p[i-1] * lam[i-1] / muv[i]
  p / sum(p)
}
get_p1 <- function(N, beta, gamma) {
  lam <- beta * (1:N) * (N - (1:N)) / N; muv <- gamma * (1:N)
  tau <- numeric(N); tau[1] <- 1
  for (i in 2:N) tau[i] <- tau[i-1] * lam[i-1] / muv[i-1]
  tau / sum(tau)
}

get_ET <- function(Q) as.vector(-solve(Q) %*% rep(1, nrow(Q)))
caract <- function(p) {
  k <- seq_along(p); m <- sum(k * p)
  list(media = m, sd = sqrt(sum((k - m)^2 * p)), p1 = p[1], moda = which.max(p))
}

R <- list()
for (r0 in R0s) {
  b  <- r0 * gamT; Q <- build_QT(NT, b, gamT)
  qs <- get_QSD(Q); bb <- get_RE(Q)
  p0 <- get_p0(NT, b, gamT); p1 <- get_p1(NT, b, gamT); ET <- get_ET(Q)
  R[[as.character(r0)]] <- list(
    Ru = qs$Ru, u = caract(qs$qsd), p1 = caract(p1),
    b1 = caract(bb[1, ]), bN = caract(bb[NT, ]),
    E1T = ET[1], ENT = ET[NT],
    d_p1 = max(abs(p1 - qs$qsd)),
    d_b1 = max(abs(bb[1, ] - qs$qsd)), d_bN = max(abs(bb[NT, ] - qs$qsd))
  )
}

#tabla1
fmtnum <- function(x) formatC(x, format = "g", digits = 5)
cat("\n--- TABLA 1 (formato Artalejo, N=50) ---\n")
hdr <- c("", sprintf("R0=%.1f", R0s))
cat(sprintf("%-10s", hdr), "\n")
cat(sprintf("%-10s", "Ru"));   for (r0 in R0s) cat(sprintf("%11s", fmtnum(R[[as.character(r0)]]$Ru))); cat("\n")
for (dist in c("u", "b1", "bN", "p1")) {
  lab <- c(u = "QSD", b1 = "b_1", bN = "b_N", p1 = "p(1)")[dist]
  cat("--", lab, "--\n")
  for (st in c("media", "sd", "p1", "moda")) {
    cat(sprintf("%-10s", paste0(lab, ".", st)))
    for (r0 in R0s) {
      v <- R[[as.character(r0)]][[dist]][[st]]
      cat(sprintf("%11s", if (st == "moda") sprintf("%d", v) else fmtnum(v)))
    }
    cat("\n")
  }
}
cat(sprintf("%-10s", "E1[T]")); for (r0 in R0s) cat(sprintf("%11s", fmtnum(R[[as.character(r0)]]$E1T))); cat("\n")
cat(sprintf("%-10s", "EN[T]")); for (r0 in R0s) cat(sprintf("%11s", fmtnum(R[[as.character(r0)]]$ENT))); cat("\n")


#tabla2
cat(sprintf("%-10s", "")); for (r0 in R0s) cat(sprintf("%11s", sprintf("R0=%.1f", r0))); cat("\n")
for (d in c("d_b1", "d_p1", "d_bN")) {
  lab <- c(d_b1 = "|b1-u|", d_p1 = "|p1-u|", d_bN = "|bN-u|")[d]
  cat(sprintf("%-10s", lab))
  for (r0 in R0s) cat(sprintf("%11s", formatC(R[[as.character(r0)]][[d]], format = "g", digits = 4)))
  cat("\n")
}

tex1 <- file("tabla1_sis.tex", "w")
writeLines(c(
  "\\begin{table}[H]", "  \\centering",
  paste0("  \\begin{tabular}{l", paste(rep("r", length(R0s)), collapse = ""), "}"),
  "    \\hline",
  paste0("    & ", paste(sprintf("$R_0=%.1f$", R0s), collapse = " & "), " \\\\"),
  "    \\hline",
  paste0("    $R_{\\mathbf u}$ & ",
         paste(sapply(R0s, function(r0) fmtnum(R[[as.character(r0)]]$Ru)), collapse = " & "), " \\\\"),
  "    \\hline"), tex1)
for (dist in c("u", "b1", "bN", "p1")) {
  lab <- c(u = "\\mathbf{u}", b1 = "\\mathbf{b}_1", bN = "\\mathbf{b}_N", p1 = "\\mathbf{p}^{(1)}")[dist]
  labg <- paste0("{", lab, "}")   # agrupado, para evitar doble subindice (b_N_1)
  rows <- c(
    sprintf("$E[%s]$", lab), sprintf("$\\sigma(%s)$", lab),
    sprintf("$%s_1$", labg),  sprintf("$%s_m$", labg))
  sts <- c("media", "sd", "p1", "moda")
  for (k in 1:4) {
    vals <- sapply(R0s, function(r0) {
      v <- R[[as.character(r0)]][[dist]][[sts[k]]]
      if (sts[k] == "moda") sprintf("%d", v) else fmtnum(v)
    })
    writeLines(paste0("    ", rows[k], " & ", paste(vals, collapse = " & "), " \\\\"), tex1)
  }
  writeLines("    \\hline", tex1)
}
e1 <- sapply(R0s, function(r0) fmtnum(R[[as.character(r0)]]$E1T))
eN <- sapply(R0s, function(r0) fmtnum(R[[as.character(r0)]]$ENT))
writeLines(c(
  paste0("    $E_1[T]$ & ", paste(e1, collapse = " & "), " \\\\"),
  paste0("    $E_N[T]$ & ", paste(eN, collapse = " & "), " \\\\"),
  "    \\hline", "  \\end{tabular}",
  "  \\caption{Modelo SIS ($N=50$, $\\gamma=1$): caracteristicas de la distribucion QS $\\mathbf{u}$, las distribuciones RE extremas $\\mathbf{b}_1$ y $\\mathbf{b}_N$, y la aproximacion $\\mathbf{p}^{(1)}$, junto con el cociente espectral $R_{\\mathbf u}$. Por el Teorema 4.6, $\\mathbf{b}_1=\\mathbf{p}^{(0)}$.}",
  "  \\label{tab:sis_caracteristicas}", "\\end{table}"), tex1)
close(tex1)

tex2 <- file("tabla2_sis.tex", "w")
writeLines(c(
  "\\begin{table}[H]", "  \\centering",
  paste0("  \\begin{tabular}{l", paste(rep("r", length(R0s)), collapse = ""), "}"),
  "    \\hline",
  paste0("    & ", paste(sprintf("$R_0=%.1f$", R0s), collapse = " & "), " \\\\"),
  "    \\hline"), tex2)
for (d in c("d_b1", "d_p1", "d_bN")) {
  lab <- c(d_b1 = "$|\\mathbf{b}_1-\\mathbf{u}|$",
           d_p1 = "$|\\mathbf{p}^{(1)}-\\mathbf{u}|$",
           d_bN = "$|\\mathbf{b}_N-\\mathbf{u}|$")[d]
  vals <- sapply(R0s, function(r0) formatC(R[[as.character(r0)]][[d]], format = "g", digits = 4))
  writeLines(paste0("    ", lab, " & ", paste(vals, collapse = " & "), " \\\\"), tex2)
}
writeLines(c(
  "    \\hline", "  \\end{tabular}",
  "  \\caption{Distancias maximas $\\max_j|\\cdot_j - u_j|$ respecto de la distribucion QS ($N=50$, $\\gamma=1$).}",
  "  \\label{tab:sis_distancias}", "\\end{table}"), tex2)
close(tex2)
cat("\nTablas LaTeX guardadas: tabla1_sis.tex, tabla2_sis.tex\n")

#5
R0_grid <- seq(0.3, 2.5, by = 0.05)
Ru_grid <- sapply(R0_grid, function(R) get_QSD(build_QT(NT, R * gamT, gamT))$Ru)
df_Ru <- data.frame(R0 = R0_grid, Ru = Ru_grid)

ggplot(df_Ru, aes(x = R0, y = Ru)) +
  geom_line(color = "#377EB8", linewidth = 0.9) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "#E41A1C") +
  annotate("text", x = 2.2, y = 1.06,
           label = "R[u] == 1", parse = TRUE, color = "#E41A1C", size = 4.5) +
  labs(
    title    = "Modelo SIS finito: cociente espectral",
    subtitle = expression(R[u] < 1 ~ "=> regimen QSD alcanzado antes que la extincion (QSD adecuada)"),
    x        = expression(R[0] == beta / gamma),
    y        = expression(R[u] == 2 * alpha / alpha[2])
  ) +
  theme_bw(base_size = 13) +
  theme(panel.grid.minor = element_blank())
ggsave("fig_sis_05_Ru.pdf", width = 7, height = 5)


#SIS INFINITO
mu_inf <- 1   # normalizamos mu = 1, de modo que R = lam

# -- QSD geometrica y RE explicita (Ejemplo 4.4); se definen antes de la tabla
qs_geom <- function(jmax, R) (1 - R) * R^((1:jmax) - 1)
RE_lineal <- function(i, jmax, R) {
  re <- numeric(jmax)
  S  <- if (i >= 2) sum((1 - R^((1:(i - 1)) - i)) / (1:(i - 1))) else 0
  den <- S - (1 - R^i) / R^i * log(1 - R)
  for (j in 1:jmax) {
    re[j] <- if (j < i) (1 - R^j) / j / den
    else       R^(j - i) * (1 - R^i) / j / den
  }
  re
}


#tabla3
i0_re <- 5; jmax_tab <- 200
cat("\n--- Tabla Cap. 3: decaimiento y existencia de la QSD (mu = 1) ---\n")
tabla_c3 <- data.frame()
for (R in c(0.3, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99, 1.0, 1.2, 1.5, 2.0)) {
  a_u   <- if (R < 1) mu_inf * (1 - R) else 0     # = alpha0 = alpha si R<1
  E_qsd <- if (R < 1) 1 / (1 - R) else NA
  E_re  <- if (R < 1) sum((1:jmax_tab) * RE_lineal(i0_re, jmax_tab, R)) else NA
  l_1   <- min(1, R^(-1)); l_5 <- min(1, R^(-5))
  existe <- if (R < 1) "Si"
  else if (R == 1) "No (a_u=0)"
  else "No (abs. no segura)"
  tabla_c3 <- rbind(tabla_c3, data.frame(
    R = R,
    a_u = round(a_u, 4), alpha0 = round(a_u, 4), alpha = round(a_u, 4),
    l_1 = round(l_1, 4), l_5 = round(l_5, 4),
    E_QSD = ifelse(is.na(E_qsd), "---", sprintf("%.3f", E_qsd)),
    E_RE  = ifelse(is.na(E_re),  "---", sprintf("%.3f", E_re)),
    existe_QSD = existe
  ))
}
print(tabla_c3, row.names = FALSE)
cat("Las tres tasas coinciden: a_u = alpha0 = alpha = mu(1-R). E[QSD] crece\n")
cat("hacia infinito mas rapido que E[RE] cuando R->1 (cf. Ejemplo 4.4).\n")
write.csv(tabla_c3, "tabla_sis_cap3.csv", row.names = FALSE)

texc3 <- file("tabla_cap3_sis.tex", "w")
writeLines(c(
  "\\begin{table}[H]", "  \\centering", "  \\begin{tabular}{rcccccc}",
  "    \\hline",
  "    $R$ & $a_{\\mathbf u}=\\alpha_0=\\alpha$ & $l_1$ & $l_5$ & $E[\\mathbf{c}]$ & $E[\\mathbf{b}_{i_0}]$ & Existe QS \\\\",
  "    \\hline"), texc3)
for (k in seq_len(nrow(tabla_c3))) {
  r  <- tabla_c3$R[k]
  au <- sprintf("%.2f", tabla_c3$a_u[k])
  l1 <- sprintf("%.3f", tabla_c3$l_1[k]); l5 <- sprintf("%.3f", tabla_c3$l_5[k])
  eq <- tabla_c3$E_QSD[k]; er <- tabla_c3$E_RE[k]
  eq <- ifelse(eq == "---", "---", eq)
  er <- ifelse(er == "---", "---", er)
  ex <- gsub("_", "\\\\_", tabla_c3$existe_QSD[k])
  ex <- ifelse(grepl("^Si$", ex), "S\\'i", ex)
  writeLines(sprintf("    %.2f & %s & %s & %s & %s & %s & %s \\\\",
                     r, au, l1, l5, eq, er, ex), texc3)
}
writeLines(c(
  "    \\hline", "  \\end{tabular}",
  "  \\caption{Modelo SIS de poblaci\\'on infinita (proceso lineal, $\\mu=1$): tasas de decaimiento, probabilidad de absorci\\'on $l_i$, medias de las distribuciones QS y RE ($i_0=5$), y existencia de la distribuci\\'on QS seg\\'un $R$. Las tres tasas coinciden, $\\alpha=\\alpha_0=a_{\\mathbf u}=\\mu(1-R)$.}",
  "  \\label{tab:sis_cap3}", "\\end{table}"), texc3)
close(texc3)


#6
R <- 0.75; jmax <- 50
df6 <- rbind(
  data.frame(estado = 1:jmax, prob = qs_geom(jmax, R),    dist = "QSD exacta"),
  data.frame(estado = 1:jmax, prob = RE_lineal(5,  jmax, R), dist = "RE (i0 = 5)"),
  data.frame(estado = 1:jmax, prob = RE_lineal(30, jmax, R), dist = "RE (i0 = 30)")
)
df6$dist <- factor(df6$dist, levels = c("QSD exacta", "RE (i0 = 5)", "RE (i0 = 30)"))

ggplot(df6, aes(x = estado, y = prob, color = dist)) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A"), name = NULL) +
  labs(
    title    = "Modelo SIS (N infinito) - QSD vs distribucion RE",
    subtitle = expression(beta == 0.75 ~ "," ~ gamma == 1 ~ "," ~ R[0] == 0.75),
    x = "Estado i (numero de infectados)", y = "Probabilidad"
  ) +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())
ggsave("fig_sis_06_infinita_qsd_RE.pdf", width = 7, height = 5.5)

cat(sprintf("\nR=0.75: media QSD geom = 1/(1-R) = %.4f\n", 1/(1-R)))
cat(sprintf("        media RE(i0=5)  = %.4f\n", sum((1:jmax) * RE_lineal(5,  jmax, R))))
cat(sprintf("        media RE(i0=30) = %.4f\n", sum((1:jmax) * RE_lineal(30, jmax, R))))


sim_lineal <- function(lam, mu, x0, tmax, xcap = 200) {
  t <- 0; x <- x0; ts <- 0; xs <- x0
  while (t < tmax && x > 0 && x < xcap) {
    li <- lam * x; mi <- mu * x; tot <- li + mi
    if (tot <= 0) break
    t <- t + rexp(1, tot)
    if (t > tmax) break
    if (runif(1) < li / tot) x <- x + 1 else x <- x - 1
    ts <- c(ts, t); xs <- c(xs, x)
  }
  data.frame(t = ts, I = xs)
}

#7
set.seed(2025)
mu <- 1; n_each <- 12; tmax <- 30; xcap <- 200
df7 <- data.frame()
for (R in c(0.6, 1.0, 1.5)) {
  lam <- R * mu
  for (k in 1:n_each) {
    tr <- sim_lineal(lam, mu, x0 = 10, tmax = tmax, xcap = xcap)
    fin <- tail(tr$I, 1)
    est <- if (fin == 0) "Absorbida"
    else if (fin >= xcap) "Escapa a infinito"
    else "En curso"
    etiqueta <- if (R < 1) " (absorcion segura)"
    else if (R == 1) " (caso critico)"
    else " (absorcion no segura)"
    df7 <- rbind(df7, data.frame(
      t = tr$t, I = tr$I, id = paste(R, k),
      R = sprintf("R0 = %.1f%s", R, etiqueta),
      estado = est))
  }
}
ggplot(df7, aes(x = t, y = I, group = id, color = estado)) +
  geom_step(alpha = 0.8, linewidth = 0.5) +
  facet_wrap(~R, ncol = 1, scales = "free_y") +
  scale_color_manual(values = c("Absorbida" = "#4DAF4A",
                                "En curso" = "#377EB8",
                                "Escapa a infinito" = "#E41A1C"), name = NULL) +
  labs(
    title = "Modelo SIS (N infinito): Trayectorias segun R0",
    x = "Tiempo", y = "N. de infectados I(t)"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())
ggsave("fig_sis_07_trayectorias_R.pdf", width = 7, height = 7)

cat("\n--- Resumen de los tres regimenes (mu = 1) ---\n")
for (R in c(0.6, 1.0, 1.5)) {
  l1 <- min(1, R^(-1)); l10 <- min(1, R^(-10))
  msg <- if (R < 1) "absorcion segura, a_u > 0 : EXISTE QSD"
  else if (R == 1) "absorcion segura, a_u = 0 : NO existe QSD"
  else "absorcion no segura : NO existe QSD"
  cat(sprintf("R = %.1f : l_1=%.3f l_10=%.3f -> %s\n", R, l1, l10, msg))
}

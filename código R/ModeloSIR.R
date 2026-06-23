# MODELO SIR

library(ggplot2)

N     <- 20
m     <- 3            # I(0)
n     <- N - m        # S(0)
beta  <- 0.03
gamma <- 0.2
R0    <- beta * n / gamma
t_max <- 200
cat(sprintf("SIR: N=%d, (m,n)=(%d,%d), beta=%.3f, gamma=%.2f -> R0=%.3f\n",
            N, m, n, beta, gamma, R0))


gillespie_sir <- function(I0, S0, beta, gamma, t_max) {
  I <- I0; S <- S0; t <- 0
  traj <- data.frame(tiempo = t, I = I, S = S)
  while (I > 0 && t < t_max) {
    tasa_inf <- beta * I * S
    tasa_rec <- gamma * I
    total    <- tasa_inf + tasa_rec
    if (total <= 0) break
    dt <- rexp(1, total); t <- t + dt
    if (runif(1) < tasa_inf / total) { I <- I + 1; S <- S - 1 } else { I <- I - 1 }
    traj <- rbind(traj, data.frame(tiempo = t, I = I, S = S))
  }
  traj
}

#1
set.seed(7)
traj_ej   <- gillespie_sir(m, n, beta, gamma, t_max)
df_puntos <- traj_ej[, c("I", "S")]
df_seg <- data.frame(
  I_ini = head(df_puntos$I, -1), S_ini = head(df_puntos$S, -1),
  I_fin = tail(df_puntos$I, -1), S_fin = tail(df_puntos$S, -1)
)
estados_fondo <- expand.grid(I = 0:N, S = 0:N)
estados_fondo <- estados_fondo[estados_fondo$I + estados_fondo$S <= N, ]

ggplot() +
  geom_point(data = estados_fondo, aes(x = I, y = S), color = "gray85", size = 1.5) +
  geom_segment(data = df_seg, aes(x = I_ini, y = S_ini, xend = I_fin, yend = S_fin),
               color = "#2166AC", linewidth = 0.6,
               arrow = arrow(length = unit(0.2, "cm"), type = "closed")) +
  geom_point(data = df_puntos, aes(x = I, y = S), color = "#2166AC", size = 2.5) +
  geom_point(data = df_puntos[1, ], aes(x = I, y = S), color = "#4DAF4A", size = 5, shape = 17) +
  geom_point(data = tail(df_puntos, 1), aes(x = I, y = S), color = "#E41A1C", size = 5, shape = 15) +
  annotate("text", x = df_puntos$I[1] + 0.5, y = df_puntos$S[1] + 0.5,
           label = "Inicio", color = "#4DAF4A", size = 4) +
  annotate("text", x = tail(df_puntos$I, 1) + 0.5, y = tail(df_puntos$S, 1) + 0.5,
           label = "Absorcion", color = "#E41A1C", size = 4) +
  scale_x_continuous(breaks = 0:N) + scale_y_continuous(breaks = 0:N) +
  labs(title = "Modelo SIR: Trayectoria en el espacio de estados (I,S)",
       x = "Infectados", y = "Susceptibles") +
  theme_bw(base_size = 13) + theme(panel.grid.minor = element_blank())
ggsave("fig_sir_01_trayectoria_2d.pdf", width = 7, height = 6)


#2
set.seed(15)
escenarios <- list(
  list(beta = 0.005, label = "R0 < 1  (extincion rapida)"),
  list(beta = 0.03,  label = "R0 > 1  (epidemia grande)")
)
df_esc <- do.call(rbind, lapply(escenarios, function(e) {
  tr <- gillespie_sir(m, n, e$beta, gamma, t_max)
  data.frame(I = tr$I, S = tr$S, escenario = e$label)
}))
ggplot(df_esc, aes(x = I, y = S, color = escenario)) +
  geom_path(linewidth = 0.9) + geom_point(size = 2) +
  geom_point(data = data.frame(I = m, S = n), aes(x = I, y = S),
             color = "#4DAF4A", size = 5, shape = 17, inherit.aes = FALSE) +
  annotate("text", x = m + 0.5, y = n + 0.3, label = "Inicio", color = "#4DAF4A", size = 4) +
  scale_color_manual(values = c("#2166AC", "#D6604D"), name = NULL) +
  scale_x_continuous(breaks = 0:N) + scale_y_continuous(breaks = 0:N) +
  labs(title = "Modelo SIR: Trayectorias en el espacio de estados (I, S)",
       x = "Infectados", y = "Susceptibles") +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom", legend.text = element_text(size = 14),
        legend.key.size = unit(1.2, "cm"), panel.grid.minor = element_blank())
ggsave("fig_sir_02_trayectorias_2d.pdf", width = 7, height = 6)




build_sir <- function(N, beta, gamma) {
  S <- list(); idx <- list(); k <- 0
  for (i in 1:N) for (j in 0:(N - i)) { k <- k + 1; S[[k]] <- c(i, j); idx[[paste(i, j)]] <- k }
  ns <- k; Q <- matrix(0, ns, ns)
  for (s in 1:ns) {
    i <- S[[s]][1]; j <- S[[s]][2]; li <- beta * i * j; mi <- gamma * i
    if (j >= 1) { d <- idx[[paste(i + 1, j - 1)]]; if (!is.null(d)) Q[s, d] <- Q[s, d] + li }
    if (i - 1 >= 1) { d <- idx[[paste(i - 1, j)]]; if (!is.null(d)) Q[s, d] <- Q[s, d] + mi }
    Q[s, s] <- -(li + mi)
  }
  list(Q = Q, S = S, idx = idx, ns = ns)
}
sir <- build_sir(N, beta, gamma)
cat(sprintf("Estados transitorios |S_T| = %d  (= N(N+1)/2 = %d)\n", sir$ns, N*(N+1)/2))


ev <- eigen(t(sir$Q)); o <- order(Re(ev$values), decreasing = TRUE)
alpha <- -Re(ev$values[o[1]]); alpha2 <- -Re(ev$values[o[2]]); Ru <- 2 * alpha / alpha2
qs <- Re(ev$vectors[, o[1]]); if (sum(qs) < 0) qs <- -qs; qs[qs < 0] <- 0; qs <- qs / sum(qs)
k10 <- sir$idx[["1 0"]]
cat(sprintf("alpha=%.5f alpha2=%.5f Ru=%.5f | QS masa en (1,0)=%.6f\n",
            alpha, alpha2, Ru, qs[k10]))


M  <- -solve(sir$Q); k0 <- sir$idx[[paste(m, n)]]
ET <- sum(M[k0, ]); re <- M[k0, ] / ET
cat(sprintf("E_(%d,%d)[T]=%.4f ; RE reparte masa en %d estados\n", m, n, ET, sum(re > 1e-6)))
dfprob <- data.frame(i = sapply(sir$S, `[`, 1), j = sapply(sir$S, `[`, 2), QS = qs, RE = re)


#3
ggplot(dfprob, aes(x = i, y = j, fill = RE)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient(low = "#F7FBFF", high = "#08306B", name = "Prob. RE") +
  scale_x_continuous(breaks = seq(0, N, 2)) + scale_y_continuous(breaks = seq(0, N, 2)) +
  labs(title = "Modelo SIR: distribucion RE sobre el espacio (I,S)",
       subtitle = bquote("Inicio (" ~ I[0] == .(m) ~ "," ~ S[0] == .(n) ~ ")," ~ R[0] == .(round(R0,2))),
       x = "Infectados (I)", y = "Susceptibles (S)") +
  theme_bw(base_size = 13) + theme(panel.grid.minor = element_blank())
ggsave("fig_sir_03_RE_mapa.pdf", width = 7, height = 5)

#4
margI_RE <- tapply(dfprob$RE, dfprob$i, sum)
margI_QS <- tapply(dfprob$QS, dfprob$i, sum)
dfmarg <- rbind(
  data.frame(I = as.integer(names(margI_RE)), prob = as.numeric(margI_RE), dist = "RE"),
  data.frame(I = as.integer(names(margI_QS)), prob = as.numeric(margI_QS), dist = "QS (degenerada)")
)
dfmarg$dist <- factor(dfmarg$dist, levels = c("QS (degenerada)", "RE"))
ggplot(dfmarg, aes(x = I, y = prob, color = dist, shape = dist)) +
  geom_line(linewidth = 0.9) + geom_point(size = 2.3) +
  scale_color_manual(values = c("#E41A1C", "#377EB8"), name = NULL) +
  scale_shape_manual(values = c(16, 17), name = NULL) +
  scale_x_continuous(breaks = seq(0, N, 2)) +
  labs(title = "Modelo SIR: marginal de I bajo la QS y la RE",
       subtitle = "La QS concentra toda la masa en (1,0); la RE la reparte",
       x = "Infectados (I)", y = "Probabilidad marginal") +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())
ggsave("fig_sir_04_RE_marginal.pdf", width = 7, height = 5)



#tabla
Ru_directo <- function(N, beta, gamma) {
  Q <- build_sir(N, beta, gamma)$Q; e <- sort(Re(eigen(Q)$values), decreasing = TRUE)
  2 * (-e[1]) / (-e[2])
}
Ru_formula <- function(beta, gamma) if (gamma <= beta) 1 else 2 * gamma / (beta + gamma)
R0_vals <- c(0.5, 1.0, 1.5, 2.0, 2.5, 4.0)
cat("\n--- Ru del SIR (N=20, gamma=0.2) ---\nR0   beta    Ru\n")
tab <- data.frame()
for (R0v in R0_vals) {
  bb <- R0v * gamma / n; ruf <- Ru_formula(bb, gamma)
  cat(sprintf("%.1f  %.4f  %.5f\n", R0v, bb, ruf))
  tab <- rbind(tab, data.frame(R0 = R0v, beta = bb, Ru = ruf))
}
texr <- file("tabla_sir_Ru.tex", "w")
writeLines(c("\\begin{table}[H]", "  \\centering", "  \\begin{tabular}{rrr}", "    \\hline",
             "    $R_0$ & $\\beta$ & $R_{\\mathbf u}=2\\alpha/\\alpha_2$ \\\\", "    \\hline"), texr)
for (k in seq_len(nrow(tab)))
  writeLines(sprintf("    %.1f & %.4f & %.4f \\\\", tab$R0[k], tab$beta[k], tab$Ru[k]), texr)
writeLines(c("    \\hline", "  \\end{tabular}",
             "  \\caption{Modelo SIR ($N=20$, $\\gamma=0.2$, $n=17$): cociente espectral $R_{\\mathbf u}$ para distintos valores de $R_0$. En todos los casos $R_{\\mathbf u}\\ge 1$, de modo que la distribuci\\'on QS nunca describe adecuadamente el sistema antes de la absorci\\'on.}",
             "  \\label{tab:sir_Ru}", "\\end{table}"), texr)
close(texr)

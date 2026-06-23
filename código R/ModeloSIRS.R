#MODELO SIRS


library(ggplot2)

N<-20; beta<-0.03; gamma<-0.2; delta<-0.1; m<-3; n<-N-m
R0<-beta*n/gamma; t_max<-300; N_traj<-1000
cat(sprintf("SIRS N=%d beta=%.3f gamma=%.2f delta=%.2f R0=%.2f\n",N,beta,gamma,delta,R0))

build <- function(N,beta,gamma,delta){
  S<-list();idx<-list();k<-0
  for(i in 1:N)for(j in 0:(N-i)){k<-k+1;S[[k]]<-c(i,j);idx[[paste(i,j)]]<-k}
  ns<-k;Q<-matrix(0,ns,ns)
  for(s in 1:ns){i<-S[[s]][1];j<-S[[s]][2];Rc<-N-i-j
  li<-beta*i*j;mi<-gamma*i;di<-delta*Rc
  if(j>=1){d<-idx[[paste(i+1,j-1)]];if(!is.null(d))Q[s,d]<-Q[s,d]+li}
  if(i>=2){d<-idx[[paste(i-1,j)]];if(!is.null(d))Q[s,d]<-Q[s,d]+mi}
  if(Rc>=1){d<-idx[[paste(i,j+1)]];if(!is.null(d))Q[s,d]<-Q[s,d]+di}
  Q[s,s]<--(li+mi+di)}
  list(Q=Q,S=S,idx=idx,ns=ns)
}
sirs<-build(N,beta,gamma,delta)
Iof<-sapply(sirs$S,`[`,1); Sof<-sapply(sirs$S,`[`,2)

ev<-eigen(t(sirs$Q)); o<-order(Re(ev$values),decreasing=TRUE)
alpha<--Re(ev$values[o[1]]); alpha2<--Re(ev$values[o[2]]); Ru<-2*alpha/alpha2
qs<-Re(ev$vectors[,o[1]]);if(sum(qs)<0)qs<--qs;qs[qs<0]<-0;qs<-qs/sum(qs)
km<-which.max(qs)
cat(sprintf("alpha=%.6f alpha2=%.6f Ru=%.5f | QS moda (%d,%d)=%.4f%%\n",
            alpha,alpha2,Ru,sirs$S[[km]][1],sirs$S[[km]][2],qs[km]*100))

M<--solve(sirs$Q)
re_from<-function(i0,j0){k<-sirs$idx[[paste(i0,j0)]];M[k,]/sum(M[k,])}
re0<-re_from(m,n)   # RE desde (3,17)

set.seed(42)
gillespie<-function(I0,S0,t_max){
  I<-I0;S<-S0;t<-0;tr<-data.frame(tiempo=t,I=I,S=S)
  while(I>0 && t<t_max){
    Rc<-N-I-S; ti<-beta*I*S; trc<-gamma*I; tw<-delta*Rc; tot<-ti+trc+tw
    if(tot<=0)break; t<-t+rexp(1,tot)
    u<-runif(1)*tot
    if(u<ti){I<-I+1;S<-S-1} else if(u<ti+trc){I<-I-1} else {S<-S+1}
    tr<-rbind(tr,data.frame(tiempo=t,I=I,S=S))
  }
  tr
}
trajs<-vector("list",N_traj); for(k in 1:N_traj) trajs[[k]]<-gillespie(m,n,t_max)


#1
set.seed(7); tej<-gillespie(m,n,t_max); dp<-tej[,c("I","S")]; dp$tiempo<-tej$tiempo
ds<-data.frame(I_ini=head(dp$I,-1),S_ini=head(dp$S,-1),
               I_fin=tail(dp$I,-1),S_fin=tail(dp$S,-1),tiempo=head(dp$tiempo,-1))
ef<-expand.grid(I=0:N,S=0:N); ef<-ef[ef$I+ef$S<=N,]
ggplot()+
  geom_point(data=ef,aes(I,S),color="gray85",size=1.5)+
  geom_segment(data=ds,aes(I_ini,S_ini,xend=I_fin,yend=S_fin,color=tiempo),
               linewidth=0.6,arrow=arrow(length=unit(0.15,"cm"),type="closed"))+
  geom_point(data=dp,aes(I,S,color=tiempo),size=1.8)+
  scale_color_gradient(low="#92C5DE",high="#08306B",name="Tiempo")+
  geom_point(data=dp[1,],aes(I,S),color="#4DAF4A",size=5,shape=17)+
  geom_point(data=tail(dp,1),aes(I,S),color="#E41A1C",size=5,shape=15)+
  annotate("text",x=dp$I[1]+0.6,y=dp$S[1]+0.5,label="Inicio",color="#4DAF4A",size=4)+
  annotate("text",x=tail(dp$I,1)+0.6,y=tail(dp$S,1)+0.5,label="Absorcion",color="#E41A1C",size=4)+
  scale_x_continuous(breaks=seq(0,N,2))+scale_y_continuous(breaks=seq(0,N,2))+
  labs(title="Modelo SIRS: Trayectoria en el espacio de estados (I,S)",x="Infectados",y="Susceptibles")+
  theme_bw(base_size=13)+theme(panel.grid.minor=element_blank())
ggsave("fig_sirs_01_trayectoria_2d.pdf",width=7,height=5)


#2
dfqs<-data.frame(I=Iof,S=Sof,qsd=qs)
ggplot(dfqs,aes(I,S,fill=qsd))+
  geom_tile(color="white",linewidth=0.3)+
  scale_fill_gradient(low="#DEEBF7",high="#08306B",name="Probabilidad",
                      labels=scales::percent_format(accuracy=0.1))+
  scale_x_continuous(breaks=seq(0,N,2))+scale_y_continuous(breaks=seq(0,N,2))+
  labs(title="Modelo SIRS: QS sobre el espacio de estados (I,S)",x="Infectados I",y="Susceptibles S")+
  theme_bw(base_size=13)+theme(panel.grid.minor=element_blank())
ggsave("fig_sirs_02_qsd_mapa.pdf",width=7,height=6)


#3
dfre<-data.frame(I=Iof,S=Sof,re=re0)
ggplot(dfre,aes(I,S,fill=re))+
  geom_tile(color="white",linewidth=0.3)+
  scale_fill_gradient(low="#DEEBF7",high="#08306B",name="Probabilidad",
                      labels=scales::percent_format(accuracy=0.1))+
  scale_x_continuous(breaks=seq(0,N,2))+scale_y_continuous(breaks=seq(0,N,2))+
  labs(title="Modelo SIRS: RE sobre el espacio de estados (I,S)",
       subtitle="Estado inicial (3,17)",x="Infectados I",y="Susceptibles S")+
  theme_bw(base_size=13)+theme(panel.grid.minor=element_blank())
ggsave("fig_sirs_02b_re_mapa.pdf",width=7,height=6)

#4
qsm<-as.numeric(tapply(qs,Iof,sum))
t_eval<-c(5,15,30,90)
cl<-lapply(t_eval,function(tt){
  Iv<-c(); for(tr in trajs){it<-max(which(tr$tiempo<=tt)); if(tr$I[it]>0)Iv<-c(Iv,tr$I[it])}
  if(!length(Iv))return(NULL)
  tb<-table(factor(Iv,levels=1:N))
  data.frame(I=1:N,prob=as.numeric(tb)/sum(tb),label=paste0("t = ",tt))
})
cl<-Filter(Negate(is.null),cl); dc<-do.call(rbind,cl)
dref<-data.frame(I=1:N,prob=qsm,label="QS")
dcp<-rbind(dc,dref)
dcp$label<-factor(dcp$label,levels=c(paste0("t = ",t_eval),"QS"))
cols<-c(colorRampPalette(c("#C7E9C0","#006D2C"))(length(t_eval)),"#E41A1C")
ggplot(dcp,aes(I,prob,color=label))+
  geom_line(linewidth=0.9,na.rm=TRUE)+geom_point(size=1.8,na.rm=TRUE)+
  scale_color_manual(values=cols,name=NULL)+
  labs(title="Modelo SIRS: Convergencia de P(I(t)=k | T>t) a la QS",
       subtitle=paste0(N_traj," trayectorias, inicio uniforme"),
       x="Numero de infectados (k)",y="Probabilidad condicional")+
  theme_bw(base_size=13)+theme(panel.grid.minor=element_blank())
ggsave("fig_sirs_03_convergencia.pdf",width=7,height=6)


#5
dvs<-c(0.05,0.1,0.2,0.5)
dd<-do.call(rbind,lapply(dvs,function(dv){
  Qd<-build(N,beta,gamma,dv); evd<-eigen(t(Qd$Q));od<-which.max(Re(evd$values))
  v<-Re(evd$vectors[,od]);if(sum(v)<0)v<--v;v[v<0]<-0;v<-v/sum(v)
  mI<-as.numeric(tapply(v,sapply(Qd$S,`[`,1),sum))
  data.frame(I=1:N,prob=mI,delta=factor(dv,levels=dvs))
}))
# Etiquetas con la delta griega real (expression) para la leyenda
lab_delta <- lapply(dvs,function(dv) bquote(delta == .(dv)))
ggplot(dd,aes(I,prob,color=delta))+
  geom_line(linewidth=0.9)+geom_point(size=2)+
  scale_color_manual(values=c("#92C5DE","#2166AC","#F4A582","#D6604D"),
                     name=NULL,labels=lab_delta)+
  labs(title=expression("Modelo SIRS: efecto de "*delta*" sobre la QS"),
       x="Numero de infectados (k)",y="P(I=k) bajo la QS")+
  theme_bw(base_size=13)+theme(legend.position="bottom",panel.grid.minor=element_blank())
ggsave("fig_sirs_04_efecto_delta.pdf",width=7,height=5)


#6
rem<-as.numeric(tapply(re0,Iof,sum))
dfm<-rbind(data.frame(I=1:N,prob=qsm,dist="QS"),
           data.frame(I=1:N,prob=rem,dist="RE (desde (3,17))"))
dfm$dist<-factor(dfm$dist,levels=c("QS","RE (desde (3,17))"))
ggplot(dfm,aes(I,prob,color=dist,shape=dist))+
  geom_line(linewidth=0.9)+geom_point(size=2.3)+
  scale_color_manual(values=c("#E41A1C","#377EB8"),name=NULL)+
  scale_shape_manual(values=c(16,17),name=NULL)+
  scale_x_continuous(breaks=seq(0,N,2))+
  labs(title="Modelo SIRS: marginal de I bajo la QS y la RE",
       subtitle="Ambas distribuciones son muy proximas (Ru<1)",
       x="Numero de infectados (k)",y="Probabilidad marginal")+
  theme_bw(base_size=13)+theme(legend.position="bottom",panel.grid.minor=element_blank())
ggsave("fig_sirs_05_QS_RE_marginal.pdf",width=7,height=5)


margI<-function(p)as.numeric(tapply(p,Iof,sum))
caract<-function(p){mI<-margI(p);ii<-1:N;md<-sum(ii*mI)
list(media=md,sd=sqrt(sum(ii^2*mI)-md^2),p1=mI[1],moda=which.max(mI))}
cq<-caract(qs); cr<-caract(re0)

#tabla 1
t1<-file("tabla_sirs_caract.tex","w")
writeLines(c("\\begin{table}[H]","  \\centering","  \\begin{tabular}{lrr}","    \\hline",
             "    & QS $\\mathbf{u}$ & RE $\\mathbf{b}_{(3,17)}$ \\\\","    \\hline",
             sprintf("    $E[I]$ & %.4f & %.4f \\\\",cq$media,cr$media),
             sprintf("    $\\sigma(I)$ & %.4f & %.4f \\\\",cq$sd,cr$sd),
             sprintf("    $P(I=1)$ & %.4f & %.4f \\\\",cq$p1,cr$p1),
             sprintf("    moda & %d & %d \\\\",cq$moda,cr$moda),
             "    \\hline","  \\end{tabular}",
             "  \\caption{Modelo SIRS ($N=20$, $\\beta=0.03$, $\\gamma=0.2$, $\\delta=0.1$): caracteristicas de la marginal de $I$ bajo la distribucion QS y la distribucion RE (estado inicial $(3,17)$). Ambas distribuciones son muy proximas.}",
             "  \\label{tab:sirs_caract}","\\end{table}"),t1)
close(t1)

#tabla 2
sir<-build(N,beta,gamma,0); Msir<--solve(sir$Q); Msirs<-M
inis<-list(c(1,1),c(1,5),c(1,10),c(3,17),c(5,10),c(5,5))
t2<-file("tabla_sirs_dist.tex","w")
writeLines(c("\\begin{table}[H]","  \\centering","  \\begin{tabular}{ccrrrr}","    \\hline",
             "    $I_0$ & $S_0$ & $E[T]_{\\text{SIR}}$ & $E[T]_{\\text{SIRS}}$ & Ratio & $\\max_j|b_j-u_j|$ \\\\",
             "    \\hline"),t2)
for(e in inis){
  ks<-sirs$idx[[paste(e[1],e[2])]]; kr<-sir$idx[[paste(e[1],e[2])]]
  re<-Msirs[ks,]/sum(Msirs[ks,]); ets<-sum(Msirs[ks,]); etr<-sum(Msir[kr,])
  writeLines(sprintf("    %d & %d & %.3f & %.3f & %.3f & %.5f \\\\",
                     e[1],e[2],etr,ets,ets/etr,max(abs(re-qs))),t2)
}
writeLines(c("    \\hline","  \\end{tabular}",
             "  \\caption{Modelo SIRS ($N=20$, $\\beta=0.03$, $\\gamma=0.2$, $\\delta=0.1$): tiempos medios de absorcion en los modelos SIR y SIRS, su cociente, y distancia maxima entre la distribucion RE (desde cada estado inicial) y la distribucion QS. Las distancias, pequenas, confirman la proximidad entre ambas distribuciones.}",
             "  \\label{tab:sirs_dist}","\\end{table}"),t2)
close(t2)

#tabla3
delta_grid <- c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5)
t3 <- file("tabla_sirs_Ru.tex", "w")
writeLines(c("\\begin{table}[H]","  \\centering","  \\begin{tabular}{rrrr}","    \\hline",
             "    $\\delta$ & $\\alpha$ & $\\alpha_2$ & $R_{\\mathbf u}=2\\alpha/\\alpha_2$ \\\\","    \\hline"),t3)
for(dv in delta_grid){
  Qd <- build(N,beta,gamma,dv)$Q
  e <- sort(Re(eigen(Qd)$values),decreasing=TRUE)
  a<--e[1];a2<--e[2];Ru<-2*a/a2
  cat(sprintf("delta=%.2f  alpha=%.6f  alpha2=%.6f  Ru=%.5f\n",dv,a,a2,Ru))
  writeLines(sprintf("    %.2f & %.6f & %.6f & %.4f \\\\",dv,a,a2,Ru),t3)
}
writeLines(c("    \\hline","  \\end{tabular}",
             "  \\caption{Modelo SIRS ($N=20$, $\\beta=0.03$, $\\gamma=0.2$): cociente espectral $R_{\\mathbf u}$ en funcion de la tasa de perdida de inmunidad $\\delta$. Para $\\delta$ pequeno el sistema se comporta como el SIR ($R_{\\mathbf u}\\ge 1$); a partir de $\\delta\\approx 0.02$ se tiene $R_{\\mathbf u}<1$ y la distribucion QS pasa a ser adecuada.}",
             "  \\label{tab:sirs_Ru}","\\end{table}"),t3)
close(t3)

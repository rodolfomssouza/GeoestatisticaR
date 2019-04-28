# Rotina para análise de geoestatística utilizando o pacote geoR ----------
# Rodolfo Marcondes Silva Souza - rodolfomssouza@gmail.com
# Versão 3



# Carregando pacotes ------------------------------------------------------
require('geoR')
require('Hmisc')
require('fields')
require('fBasics')
require('plot3D')
require('viridis')

# Funções -----------------------------------------------------------------
mesc = function(x){
    if(is.null(x)){return(c('Esf', 'Exp', 'Gau', 'Sent*'))}
    else if(x=='spherical'){return(c('Esf*', 'Exp', 'Gau', 'Sent'))}
    else if(x=='exponential'){return(c('Esf', 'Exp*', 'Gau', 'Sent'))}
    else if(x=='gaussian'){return(c('Esf', 'Exp', 'Gau*', 'Sent'))}
    else{return(c('Esf', 'Exp', 'Gau', 'Sent'))}
}


por_area = function(m, xl){
    # m = matrix ou vetor de dados
    # xl = vetor com os classes/ níveis de corte
    nv = length(m)
    nc = length(xl) - 1
    tnc = rep(NA, nc)
    cls = rep(NA, nc)
    for(i in 1:nc){
        cr = c(xl[i], xl[i+1])
        if(i<nc){
            cls[i] = paste(round(cr[1],3),' <- C < ',round(cr[2],3), sep='')
            tnc[i] = round(100*(length(which(cr[1]<=m&m<cr[2]))/nv),2)
        }
        else{
            cls[i] = paste(round(cr[1],3),' < C ', sep='')
            # tnc[i] = round(100*(length(which(cr[1]>=m))/nv),2)
            tnc[i] = (100 - sum(tnc, na.rm=T))
        }
        
    }
    tnc = data.frame(cls, tnc)
    colnames(tnc) = c('Classes', 'Porcentagem')
    return(tnc)
}



# Alterando diretório de trabalho -----------------------------------------
setwd(dir = '~/')
setwd('Programming/R/GeoestatísticaYT/')  # Mude para o seu diretório de trabalho


# Carregando os dados -----------------------------------------------------
dgeo=read.table('Dados_geo.txt', h=T)
attach(dgeo); names(dgeo)

# Nome inicial para salver os arquivos
nome_analise = 'Erosividade'


# Carregando os dados no módulo de geoestatística (geoR) ------------------
dGeo=as.geodata(dgeo, coords.col = 1:2, data.col = 3)


#  Análise estatística ----------------------------------------------------
estb=basicStats(Z, ci=0.95); estb
tks=ks.test(Z, 'pnorm', mean=mean(Z), sd=sd(Z)); tks
shapiro.test(Z)
jarqueberaTest(Z)
boxplot(Z)


# Iniciando geoestatística ------------------------------------------------
# Máxima distância
mx=max(X)
my=max(Y)
mdta=sqrt(mx^2+my^2);mdta; mdta/3


# Gerando e plotando o semivariograma -------------------------------------
plot(variog(dGeo))
mdt=140; nlg=6 # máxima distância e número de logs
svgteorico=variog(dGeo, max.dist=mdt, uvec=nlg)
svgteorico1=variog(dGeo, max.dist=mdt, uvec=nlg, direction=pi/8)
svgteorico2=variog(dGeo, max.dist=mdt, uvec=nlg, direction=pi/4)
svgteorico3=variog(dGeo, max.dist=mdt, uvec=nlg, direction=pi/2)
par(mfrow=c(2,2))
plot(svgteorico); plot(svgteorico1, main=expression(pi/8)); plot(svgteorico2, main=expression(pi/4)); plot(svgteorico3, main=expression(pi/2)); layout(1)
h=svgteorico$u
v=svgteorico$v
npar=svgteorico$n
ltmx=(max(h)+0.4*max(h))
ltmy=(max(v)+0.4*max(v))
plot(svgteorico, las=1, xaxs='i', yaxs='i',pch=16, col='red', ylim=c(0,ltmy), xlim=c(0,ltmx))
textxy(h,v,npar, cex=0.7)
sm=data.frame(h,v,npar, row.names=NULL);sm


# Fazendo ajuste do semivariograma teórico --------------------------------
esferico=variofit(svgteorico, cov.model='sph', max.dist=max(h), messages=F)
exponencial=variofit(svgteorico, cov.model='exponential', max.dist=max(h), messages=F)
gaussiano=variofit(svgteorico, cov.model='gaussian', max.dist=max(h),messages=F)
X11(); sentimento=eyefit(svgteorico, silent=F)
stp = unlist(sentimento)
sigmasq=as.numeric(stp[2]); phi=as.numeric(stp[3]); tausq=as.numeric(stp[4]);


# Plotando todos o semivariograma e seus ajustes --------------------------
par(mfrow=c(2,2),mar=c(5,5,2,2))
plot(svgteorico, las=1, type='p',pch=19, cex=1.4, col='black' , xlab='', ylim=c(0,ltmy), xlim=c(0,ltmx),ylab='', xaxs='i', yaxs='i', cex.lab=1.4, cex.axis=1.2, main='Esférico')
lines.variomodel(esferico, col='red' , lwd=2, lty=1)
plot(svgteorico, las=1, type='p',pch=19, cex=1.4, col='black' , xlab='', ylim=c(0,ltmy), xlim=c(0,ltmx),ylab='', xaxs='i', yaxs='i', cex.lab=1.4, cex.axis=1.2, main='Exponencial')
lines.variomodel(exponencial, col='red' , lwd=2, lty=1)
plot(svgteorico, las=1, type='p',pch=19, cex=1.4, col='black' , xlab='', ylim=c(0,ltmy), xlim=c(0,ltmx),ylab='', xaxs='i', yaxs='i', cex.lab=1.4, cex.axis=1.2, main='Gaussiano')
lines.variomodel(gaussiano, col='red' , lwd=2, lty=1)
plot(svgteorico, las=1, type='p',pch=19, cex=1.4, col='black' , xlab='', ylim=c(0,ltmy), xlim=c(0,ltmx),ylab='', xaxs='i', yaxs='i', cex.lab=1.4, cex.axis=1.2, main='Sentimento')
lines(sentimento, col='red' , lwd=2, lty=1)


dev.copy(pdf, paste(nome_analise, '_Semivariogramas.pdf', sep = ''),  width=8, height=6)
dev.off()
layout(1)


# Parâmetros ajustados para os modelos ------------------------------------ 
c0.esf=esferico$nugget; c.esf=esferico$cov.pars[1]; a.esf=esferico$cov.pars[2]
c0.exp=exponencial$nugget; c.exp=exponencial$cov.pars[1]; a.exp=exponencial$cov.pars[2]
c0.gau=gaussiano$nugget; c.gau=gaussiano$cov.pars[1]; a.gau=gaussiano$cov.pars[2]
C0.snt = sigmasq; C.snt = tausq; a.snt = phi
ide.esf = 100*(c0.esf/(c0.esf+c.esf)); ide.exp = 100*(c0.exp/(c0.exp+c.exp)); ide.gau = 100*(c0.gau/(c0.gau+c.gau)); ide.sent = 100*(C.snt/(C0.snt+C.snt))

# Validação cruzada dos ajustes -------------------------------------------
cv.esf=xvalid(dGeo, model=esferico); zsco.esf=cv.esf$std.error
cv.exp=xvalid(dGeo, model=exponencial); zsco.exp=cv.exp$std.error
cv.gau=xvalid(dGeo, model=gaussiano); zsco.gau=cv.gau$std.error
cv.sent=xvalid(dGeo, model=sentimento); zsco.sent=cv.sent$std.error


# Média e variância do erro reduzido --------------------------------------
jkmed.esf=round(mean(zsco.esf),5); jkvar.esf=round(var(zsco.esf),5)
jkmed.exp=round(mean(zsco.exp),5); jkvar.exp=round(var(zsco.exp),5)
jkmed.gau=round(mean(zsco.gau),5); jkvar.gau=round(var(zsco.gau),5)
jkmed.sent=round(mean(zsco.sent),5); jkvar.sent=round(var(zsco.sent),5)


# R2 krigagem -------------------------------------------------------------
r2k.esf = cor(cv.esf$predicted,cv.esf$data)^2
r2k.exp = cor(cv.exp$predicted,cv.exp$data)^2
r2k.gau = cor(cv.gau$predicted,cv.gau$data)^2
r2k.sent = cor(cv.sent$predicted,cv.sent$data)^2


# Resumo das análises para escolha do melhor modelo -----------------------
modelos=c('Esf', 'Exp','Gau','Sent')
m.jk=rbind(jkmed.esf,jkmed.exp, jkmed.gau, jkmed.sent)
v.jk=rbind(jkvar.esf, jkvar.exp, jkvar.gau, jkvar.sent)
c0.smfit=rbind(c0.esf, c0.exp, c0.gau, tausq)
c.smfit=rbind(c.esf, c.exp, c.gau, sigmasq)
a.smfit=rbind(a.esf, a.exp, a.gau, phi)
ide = rbind(ide.esf, ide.exp, ide.gau, ide.sent)
r2k = rbind(r2k.esf, r2k.exp, r2k.gau, r2k.sent)

resumo=data.frame(row.names=modelos, c0.smfit, c.smfit, a.smfit, ide, m.jk, v.jk, r2k)
resumo


# Escolha do melhor modelo ------------------------------------------------
# esferico; exponencial; gaussiano; sentimento
smfit = esferico

row.names(resumo) = mesc(smfit$cov.model)
write.table(x = resumo, file = paste(nome_analise, '_Resumo_Geo.txt', sep = ''))


# Gerando o grid para interpolações ---------------------------------------
ndiv = 200 # Tamanho do intervalo para interpolação
x.range <- as.integer(range(X))
y.range <- as.integer(range(Y))
grid.map=expand.grid(x=seq(from=x.range[1], to=x.range[2], by= (x.range[2] - x.range[1])/ndiv),
                     y=seq(from=y.range[1], to=y.range[2], by=(y.range[2] - y.range[1])/ndiv))



# Carregando limites da borda ---------------------------------------------
lmt=read.table('Dados_Contorno.txt', h=T)
dlmt=read.geodata('Dados_Contorno.txt', h=T, coords.col=1:2, data.col=NULL)


# Fazendo Krigagem --------------------------------------------------------
# krg=krige.conv(dGeo, locations=grid.map, krige=krige.control(obj.model=smfit))
krg=krige.conv(dGeo, locations=grid.map, krige=krige.control(obj.model=smfit), borders=dlmt)


# Gráfico de contornos ----------------------------------------------------
contour(krg, f=T, col=viridis(10), nlevels=10)


# Configurações para salvar o semivariograma ------------------------------
par(family='serif', las=1, xaxs='i', yaxs='i', mar=c(5,5,2,2), cex.axis = 1.0, cex.lab = 1.2)
plot(svgteorico, type='p',pch=19,col='black' , ylim=c(0,(max(sm$v)+max(sm$v)*0.15)), xlim=c(0,(max(sm$h)+max(sm$h)*0.15)), ylab='', xlab='h (m)')#, yaxt='n')
vr=var(dgeo$Z)
abline(v=NULL, h=vr, lty=2, lwd=1, untf=3)
lines(smfit, col='black' , lwd=2, lty=1)
mtext(text=expression(gamma(h)), side=2, line=3, las=3, cex=1.2)
dev.copy(pdf, paste(nome_analise, '_smg.pdf', sep = ''), width=8, height=6)
dev.off()


# Configurações para o mapa -----------------------------------------------
numlevel=9
# Sequências para a legenda
sl=seq(min(krg$predict), max(krg$predict),by=(max(krg$predict)-min(krg$predict))/numlevel)
sll=formatC(sl, digits=2, format='f', decimal.mark = ".")
zlia = range(sl)

xl = c(500, 750)      # Intervalo do eixo X
yl = c(9000, 9250)    # Intervalo do eixo Y

# Mapa
par(family='serif', las=1, xaxs='i', yaxs='i', mar=c(5,5,2,6.5), cex.axis = 1.0, cex.lab = 1.2)
plot(dgeo$X,dgeo$Y, type='n', xlab='X (UTM)', ylab='Y (UTM)', yaxt = 'n',
     xlim = xl, ylim = yl)
axis(side = 2, at = seq(min(Y), max(Y), (max(Y)-min(Y))/5), las = 3)
image(krg, add=T, col=rev(viridis(numlevel)), zlim=zlia)
contour(krg, add=T, levels=sl, labcex=0.5, lwd=1, labels=sll, drawlabels=F)
lines(lmt, lwd=2)
box()
points(dgeo$X,dgeo$Y, pch='+', col='red')
colkey(side = 4, length = 0.7, padj = 0.5, shift = 0, dist = -0.01, add = T, cex.axis = 0.9, cex.clab = 1,
       at = sl, mgp = c(1, 0.3, 0), tcl = -0.2, clim=zlia, labels = sll, col = rev(viridis(numlevel)),
       side.clab = 2, line.clab = -4.2, clab = 'Variavel')

dev.copy(pdf, paste(nome_analise, '_mapa.pdf', sep = ''), width=8, height=8)
dev.off()


# Porcentagens das classes do mapa
pcm = por_area(m = krg$predict, xl = sl); pcm
write.csv(x = pcm, file = paste(nome_analise, '_classes_mapa.csv', sep = ''))



# Mapa do desvio padrão
esc.var=matrix(krg$krige.var, ncol=1)
esc.desv=matrix(sqrt(krg$krige.var), ncol=1)
#numlevel=4

# Sequências para a legenda
slsd=seq(min(esc.desv), max(esc.desv),by=(max(esc.desv)-min(esc.desv))/numlevel)
sllsd=formatC(slsd, digits=2, format='f', decimal.mark = ".")
zlsd = range(slsd)

# Mapa
par(family='serif', las=1, xaxs='i', yaxs='i', mar=c(5,5,2,6.5), cex.axis = 1.0, cex.lab = 1.2)
plot(dgeo$X,dgeo$Y, type='n', xlab='X (UTM)', ylab='Y (UTM)', yaxt = 'n',
     xlim = xl, ylim = yl)
axis(side = 2, at = seq(min(Y), max(Y), (max(Y)-min(Y))/5), las = 3)
image(krg, val=(krg$krige.var)^0.5, add=T, col=rev(viridis(numlevel)), zlim=c(min(esc.desv),max(esc.desv)))
lines(lmt, lwd=2)
box()
points(dgeo$X,dgeo$Y, pch='+', col='red')

colkey(side = 4, length = 0.7, padj = 0.5, shift = 0, dist = 0.05, add = T, cex.axis = 0.9, cex.clab = 1,
       at = slsd, mgp = c(1, 0.3, 0), tcl = -0.2, clim=zlsd, labels = sllsd, col = rev(viridis(numlevel)),
       side.clab = 2, line.clab = -4.2, clab = 'Desvio padrão')
dev.copy(pdf, paste(nome_analise, '_mapa_desvio_padrao.pdf', sep = ''), width=8, height=8)
dev.off()

# Porcentagem das classes do desvio padrão
pcmsd = por_area(m = (krg$krige.var)^0.5, xl = slsd); pcmsd
write.csv(x = pcmsd, file = paste(nome_analise, '_classes_mapa_desvio.csv', sep = ''))


# Salvando sessão do R ----------------------------------------------------
save.image(paste(nome_analise, '.RData', sep = ''))

# Fim

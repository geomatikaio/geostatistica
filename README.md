# geostatistica
Essa pasta é feita para guardar codigos interessantes de geoestatítica.

# Código em R

#GeoestatÝstica_1.0-Modelagem Espacial#

#Mudar para o diret¾rio de trabalho#

#Carregar os pacotes necessßrios(Programa R.3.1.1)
library(fBasics)
library(geoR)
library(gstat)
library(lattice)
library(sp)

#carregar arquivo de dados
Geo <- read.table("jura.txt", head=T)
names(Geo) #mostrar o nome das varißveis do banco Geo

#Geo	#mostrar os dados (varißveis) na area de trabalho
#Resumo do arquivo Geo
summary(Geo)

#Selecionar as variaveis quantitativas
Cd <- as.geodata(Geo,coords.col = 1:2, data.col = 5)
Co <- as.geodata(Geo,coords.col = 1:2, data.col = 6)
Ni <- as.geodata(Geo,coords.col = 1:2, data.col = 9)

#Mapa de pontos
#par(mfrow = c(1,2))

#definfindo o tamanho dos caracteres
x11()
par(cex = .75)
points(Cd, xlab = "W-E", main = "Cadmio", ylab = "N-S", cex.min = 1.2,cex.max =1.5, pt.divide = "quart", x.leg=c(0.2,1.8),y.leg=c(4.3,5.7))

x11()
par(cex = .75)
points(Co, xlab = "W-E", main = "Cobalto", ylab = "N-S", cex.min = 1.2,cex.max =1.5, pt.divide = "quart", x.leg=c(0.2,1.8),y.leg=c(4.3,5.7))

x11()
par(cex = .75)
points(Ni, xlab = "W-E", main = "NÝquel", ylab = "N-S", cex.min = 1.2, cex.max=1.5, pt.divide = "quart", x.leg=c(0.2,2.3), y.leg=c(4.3,5.7))


###VerificaþÒo da estacionaridade (ausencia de tendencia)
#Cadmio
#par(mfrow = c(1,2))
x11()
#par(mfrow = c(1,2)) # set up the graphics
plot(Cd$coords[, 1], Cd$data, xlab = "W-E", ylab = "Cadmio", pch = 20, cex = 1.5, main = "Cadmio")
lines(lowess(Cd$data ~ Cd$coords[, 1]))

#Cobalto
#par(mfrow = c(1,2))
x11()
#par(mfrow = c(1,2)) # set up the graphics
plot(Co$coords[, 1], Co$data, xlab = "W-E", ylab = "Cobalto", pch = 20, cex = 1.5, main = "Cobalto")
lines(lowess(Co$data ~ Co$coords[, 1]))

#NÝquel
#par(mfrow = c(1,2))
x11()
#par(mfrow = c(1,2)) # set up the graphics
plot(Ni$coords[, 1], Ni$data, xlab = "W-E", ylab = "NÝquel", pch = 20, cex = 1.5, main = "NÝquel")
lines(lowess(Ni$data ~ Ni$coords[, 1]))


###Semivariogramas experimentais direcionais
#par(mfrow = c(1,2))

#Variavel cadmio

v4Cd <- variog4(Cd, uvec = seq(0, 2, length = 20), direction = c(0, 45, 90, 135), tolerance = 45,
unit.angle = "degrees")

x11()
plot(v4Cd)
title("Cadmio")

#Variavel Cobalto

v4Co <- variog4(Co, uvec = seq(0, 2, length = 20), direction = c(0, 45, 90, 135), tolerance = 45,
unit.angle = "degrees")

x11()
plot(v4Co)
title("Cobalto")


#Variavel NÝquel

v4Ni <- variog4(Ni, uvec = seq(0, 2, length = 20), direction = c(0, 45, 90, 135), tolerance = 45,
unit.angle = "degrees")

x11()
plot(v4Ni)
title("NÝquel")

###Semivariograma de superficie
#Variavel Cadmio
trellis.par.set(sp.theme()) 
vgcd = gstat(NULL, "Cadmio", I(Cd)~1, Geo)
x11()
plot(variogram(vgcd, cutoff=2, width=.2, map=TRUE),cex=1.5,
    main = "(cross) semivariance maps: Cadmio")

#Variavel Cobalto
trellis.par.set(sp.theme()) 
vgco = gstat(NULL, "Cobalto", I(Co)~1, Geo)
x11()
plot(variogram(vgco, cutoff=2, width=.2, map=TRUE),cex=1.5,
    main = "(cross) semivariance maps: Cobalto")


#Variavel NÝquel
coordinates(Geo) = ~ x + y;
trellis.par.set(sp.theme())
vgni = gstat(NULL, "NÝquel", I(Ni)~1, Geo)
x11()
plot(variogram(vgni, cutoff=2, width=.2, map=TRUE),cex=1.5,
    main = "(cross) semivariance maps: NÝquel")

###Modelamento semivariografico
###Analise semivariografica: Niquel
#par(mfrow = c(1,2)) # set up the graphics
variogNi <- variog(Ni, uvec = seq(0, 2, length = 20), option = "bin")
#Semivariograma unidirecional experimental
plot(variogNi, xlab = "DistÔncia", ylab = "SemivariÔncia",type="l")
title("Semivariograma Experimental: Niquel")

#Semivariograma unidirecional e modelo ajustado
plot(variogNi, xlab = "DistÔncia", ylab = "SemivariÔncia")
lines.variomodel(cov.model = "spherical" , cov.pars = rbind(c(9,0.12),c(64,1.4)), nugget = 11, max.dist = 2, lwd = 2, col = "red")
text(1,0.1, "MODELO: 11 + 9,0 Sph(h/0,12) + 64 Sph (h/1,4)")
title("Semivariograma Experimental e modelo ajustado: Niquel")


#Analise semivariografica: Cobalto
#Semivariogramas experimentais nas duas direþ§es
par(mfrow = c(1,2))
coordinates(Geo) = ~ x + y
vgmCo <- variogram(Co~1, Geo,coords.col = 1:2, alpha=c(67.5,157.5))
plot(vgmCo,type="b", lty=2, col="black",main="Var. Exper. Cobalto")

#Modelo anisotropico ajustado
par(mfrow = c(1,2))
model.Co <- vgm(11.5,"Sph",1.5,2.0,anis=c(67.5,0.67))
plot(vgmCo,model=model.Co,type="p", col=c("black"),main="Modelo ajustado Cobalto")

##dados originais
par(mfrow = c(1,3))
cd.vario <- variog(Cd, lam=0, max.dist=6)
plot(cd.vario)
title("Cadmio")
co.vario <- variog(Co, lam=0, max.dist=6)
plot(cd.vario)
title("Cobalto")
ni.vario <- variog(Ni, lam=0, max.dist=6)
plot(ni.vario)
title("NÝquel")


##ou entÒo use a funþÒo interativa eyefit()
co.ef <- eyefit(co.vario)

cd.ef <- eyefit(cd.vario)

ni.ef <- eyefit(ni.vario)


##ajustando (estimando parÔmetros) usando mÝnimos quadrados
args(variofit)

##minimos quadrados ordinßrios
cd.ols <- variofit(cd.vario, ini=c(1.75, 2.9), wei="equal", min="optim")
cd.ols

##Estimando parametros por maxima verossimilhanca

cd.ml <- likfit(Cd, ini=c(1.75,2.9), lambda=1)
cd.ml
summary(cd.ml)
plot(cd.vario)
lines(cd.ml, lty=2, col="blue", lwd=2, type="l")




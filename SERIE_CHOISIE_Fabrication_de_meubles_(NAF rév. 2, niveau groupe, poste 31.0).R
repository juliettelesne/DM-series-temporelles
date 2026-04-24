#A FAIRE : RAPPELER LES HYPOTHESES POUR KPSS ADF PP
#DEMANDER A REDA SI CA VAUT LE COUP DE METTRE BOX-PIERCE (NORMALEMENT, IL FAUT UN ECHANTILLON DE TAILLE CONSIDERABLE POUR QUE CE TEST SOIT VALABLE.)
#CONFIRMER AVEC REDA LE CHOIX DE L'AR2 (sur la base des acfs et pacfs)
#REGARDER CE QU'IL SE PASSE SUR LA MINIMISATION DE LA MIN SQUARE ERROR TRAINING VS TEST
#FINIR LA PARTIE MODELISATION : GROS DEBAT A PARTIR DE LA 310.
#Données accessibles via ce lien : https://www.insee.fr/fr/statistiques/serie/010768180#Telechargement

require(zoo) #format de serie temporelle pratique et facile d'utilisation (mais plus volumineux)
require(tseries) #diverses fonctions sur les series temporelles
library(fUnitRoots)
#PREMIERE ETAPE : RENDRE LA SERIE EXPLOITABLE (13)
#DEUXIEME ETAPE : RENDRE LA SERIE TEMPORELLE STATIONNAIRE (43)
#TROISIEME ETAPE : SOUS L'HYPOTHESE D'UNE SERIE STATIONNAIRE AUTOUR D'UNE TENDANCE LINEAIRE (72)
#QUATRIEME ETAPE : SOUS L'HYPOTHESE D'UNE RACINE UNITAIRE AVEC TENDANCE (342)

#PREMIERE ETAPE : RENDRE LA SERIE EXPLOITABLE

getwd() #affiche le wd
list.files() #liste les éléments du wd

datafile <- "valeurs_mensuelles.csv" #définit le fichier de données

data <- read.csv(datafile,sep=";") #importe un fichier .csv dans un objet de classe data.frame

fm.source <- zoo(data[[2]]) #On convertit le deuxième élément de data en serie temporelle de type "zoo".
print(fm.source) #On observe que les premières valeurs sont les plus récentes.

fm.source_rev <- rev(fm.source) #On inverse l'ordre des valeurs pour obtenir un ordre chronologique.
print(fm.source_rev) #On observe que les 3 dernières valeurs du vecteur ne sont pas des valeurs d'indices.

T <- length(fm.source_rev) #On prépare la suppression de ces trois dernières valeurs.
fm.raw <- fm.source_rev[1:(T-3-72)] #On supprime les 3 dernières valeurs qui ne correspondent pas à des indices, puis les 72 dernières afin d'étudier les années précédents la crise du Covid et les différentes guerres lesquelles annihilent l'hypothèse de stationnarité.

str(fm.raw) #On vérifie si les valeurs rentrées sont des caractères ou des nombres : ce sont des caractères.
fm.num <- as.numeric(fm.raw) #On convertit les caractères en nombres.
str(fm.num) #On vérifie qu'il s'agit bien du format numérique, nécessaire au traitement.
class(fm.num) #On observe que la série n'est plus au format "zoo".

fm <- zoo(fm.num) #On reconvertit en série temporelle de type "zoo".
class(fm) #On vérifie que l'on fait à nouveau face à une série temporelle de type "zoo".
print(fm) #On n'observe pas d'anomalie.







#DEUXIEME ETAPE : RENDRE LA SERIE TEMPORELLE STATIONNAIRE

#Premier tracé de la série
par(mfrow=c(1,1))
plot(fm)
#On présume l'existence d'une tendance décroissante. Il est cependant nécessaire de vérifier l'existence d'une racine unitaire.
#On présume l'existence d'une racine unitaire à cause de la forte inertie autour de ce qui semble être une "tendance" linéaire décroissante.

#On affiche les ACF et PACF.
par(mfrow=c(1,2))
pacf(fm, lag = 50)# On observe une PACF d'ordre 1 très proche de 1, ce qui nous fait suspecter l'existence d'une racine unitaire.
acf(fm, lag = 50) # On observe une très forte persistance des ACF, ce qui nous maintient dans notre suspicion.

#Test de stationnarité de la série

#KPSS : valable dans un cas très général
#H_0 : stationnarité
kpss.test(fm, null = "Level") # 
kpss.test(fm, null = "Trend") #On se fie à ce test, sachant que l'on observe une tendance linéaire.
#La p-valeur du test de KPSS est inférieure à 0.01 donc ce test rejette l'hypothèse de stationnarité.

#Test d'existence d'une racine unitaire.

#Perron-Phillips : valable dans un cas très général 
#H_0 : non-stationnarité
#Le test de Perron-Phillips est robuste aux ruptures, pas KPSS.
pp.test(fm)
#La p-valeur du test de Perron-Phillips est supérieure à 0.05, donc ce test ne rejette pas l'existence d'une racine unitaire.

#ADF
# H0 : existence d'une racine unitaire
Qtests <- function(series, k, fitdf=0) {
  pvals <- apply(matrix(1:k), 1, FUN=function(l) {
    pval <- if (l<=fitdf) NA else Box.test(series, lag=l, type="Ljung-Box", fitdf=fitdf)$p.value
    return(c("lag"=l,"pval"=pval))
  })
  return(t(pvals))
}
series <- fm; kmax <- 24; adftype="ct"
adfTest_valid <- function(series, kmax, adftype){
  k <- 0
  noautocorr <- 0
  while (noautocorr==0){
    cat(paste0("ADF with ",k," lags: residuals OK? "))
    adf <- adfTest(series, lags=k, type=adftype)
    pvals <- Qtests(adf@test$lm$residuals, 24, fitdf = length(adf@test$lm$coefficients))[,2]
    if (sum(pvals<0.05,na.rm=T)==0) {
      noautocorr <- 1; cat("OK \n")
    } else cat("nope \n")
    k <- k+1
  }
  return(adf)
}

#Etant donné le tracé de la série, on présume l'existence d'une tendance linéaire, ce qui nous fait choisir la spécification 'ct'.
adf <- adfTest_valid(fm,30,adftype="ct")

#On vérifie l'hypothèse de résidus blancs avec une série de tests Ljung-Box.
Qtests(adf@test$lm$residuals, 30, fitdf = length(adf@test$lm$coefficients))
#Les résidus sont bien décorrélés.

adf
##La p-valeur du test d'Augmented Dickey-Fuller est supérieure à 0.10 donc ce test ne rejette pas l'existence d'une racine unitaire.

#Bilan : on part du principe que la série a une racine unitaire, les tests étant unanimes.


















#On peut tout de même essayer de retrancher à la série ce qui semble être une tendance linéaire et vérifier que cette série transformée a une racine unitaire.

#On extrait la tendance linéaire.
tf <- time(fm)
model <- lm(fm ~ tf) 

#Série temporelle des résidus du modèle (obtenue par retranchement du modèle à la série initiale.)
resf <- model$residuals

#Tracé de cette série
par(mfrow=c(1,2))
plot(fm, col = "blue", lwd = 2)
abline(model, col = "red", lwd = 2) #On a a priori une tendance décroissante.
plot(resf) #Les résidus sont approximativement centrés de variance approximativement stable avec quelques outliers, mais on observe une très forte inertie, ce qui laisse supposer l'existence d'une racine unitaire.

#On affiche les ACFs et PACFs
par(mfrow=c(1,2))
pacf(resf, lag = 50)#
acf(resf, lag = 50) #On observe une très forte persistance des ACFs, ce qui laisse toujours présumer de l'existence d'une racine unitaire.

#Test de stationnarité de la série

#KPSS : valable dans un cas très général
#H_0 : stationnarité
kpss.test(resf, null = "Level") # On se fie à ce test, sachant que les résidus n'ont pas tendance par définition.
kpss.test(resf, null = "Trend") #
#La p-valeur du test de KPSS est inférieure à 0.01 donc ce test rejette l'hypothèse de stationnarité.

#Test d'existence d'une racine unitaire.

#Perron-Phillips : valable dans un cas très général 
#H_0 : non-stationnarité
#Le test de Perron-Phillips est robuste aux ruptures, pas KPSS.
pp.test(resf)
#La p-valeur du test de Perron-Phillips est supérieure à 0.05, donc ce test ne rejette pas l'existence d'une racine unitaire.

#ADF
# H0 : existence d'une racine unitaire

#Etant donné le tracé de la série des résidus, on choisit la spécification 'c'.
adf <- adfTest_valid(resf, 30,adftype="c")

#On vérifie l'hypothèse de résidus blancs avec une série de tests Ljung-Box.
Qtests(adf@test$lm$residuals, 30, fitdf = length(adf@test$lm$coefficients))
#Les résidus sont bien décorrélés.

adf
##La p-valeur du test d'Augmented Dickey-Fuller est supérieure à 0.10 donc ce test ne rejette pas l'existence d'une racine uniatire.

#Bilan : on observe que la majorité des tests prônent l'existence d'une racine unitaire, et on observe une très forte inertie de la série des résidus du modèle linéaire.








#Bilan : il nous faut différencier la série pour espérer la stationnariser.

#Ce faisant, la série est différenciée.
dfm <- diff(fm)

#Tracé de la série différenciée
par(mfrow=c(1,1))
plot(dfm) # La série différenciée semble osciller autour de 0 en évoluant dans un couloir, donc semble stationnaire.

#Tracé des ACF et PACF
par(mfrow=c(1,2))
pacf(dfm)#Les PACFs s'annulent rapidement, ce qui nous faire présumer que la série différenciée est bien stationnaire.
acf(dfm) #Il en est de même pour les ACF.

#Test de stationnarité de la série

#KPSS : valable dans un cas très général
#H_0 : stationnarité
kpss.test(dfm, null = "Level") #On ne prend en compte que ce test, sachant que l'on observe pas de tendance linéaire, mais des oscillations autour d'une constante.
kpss.test(dfm, null = "Trend") #
#La p-valeur du test de KPSS est supérieure à 0.10 donc ce test ne rejette pas l'hypothèse de stationnarité.

#Test de l'existence d'une racine unitaire

#Perron-Phillips : valable dans un cas très général 
#H_0 : non-stationnarité
#Le test de Perron-Phillips est robuste aux ruptures, pas KPSS.
pp.test(dfm)
#La p-valeur du test de Perron-Phillips est inférieure à 0.01 donc ce test rejette l'existence d'une racine unitaire.

#On n'observe pas de tendance nette, donc on choisit une spécification avec une simple constante pour le test d'Augmented Dickey-Fuller.
adf <- adfTest_valid(dfm,24,"c")

#On vérifie l'hypothèse de résidus blancs avec une série de tests Ljung-Box.
Qtests(adf@test$lm$residuals, 24, fitdf = length(adf@test$lm$coefficients))

# Les résidus sont biens décorrélés.
adf
##La p-valeur du test d'Augmented Dickey-Fuller est inférieure à 0.01 donc ce test rejette l'existence d'une racine unitaire.

#Bilan : la série différenciée semble bien stationnaire et les tests sont unanimes.

#Pour stationnariser la série initiale, on la différencie donc une fois.

#Tracé de la série fm et de sa transformation stationnaire (série différenciée) dfm
par(mfrow=c(2,1))
plot(fm)
plot(dfm)







#Modélisation de la série différenciée

signif <- function(estim){ #fonction de test des significations individuelles des coefficients d'un modèle
  coef <- estim$coef
  se <- sqrt(diag(estim$var.coef))
  t <- coef/se
  pval <- (1-pnorm(abs(t)))*2
  return(rbind(coef,se,pval))
}


#Tracé de la série différenciée
plot(dfm)
#On affiche les ACF et PACF
par(mfrow=c(1,2))
pacf(dfm)# Au-delà du second ordre, les PACFs s'annulent.
acf(dfm) # Au-delà du premier ordre, les ACFs s'annulent.

#Ce faisant, on essaye de modéliser la série différenciée par un ARMA(2,1)
arima201 <- arima(dfm,c(2,0,1)) #

# On vérifie que les résidus du grand modèle sont bien blancs, i.e. que le modèle est valide
Qtests(arima201$residuals, 24,3)#tests de LB pour les ordres 1 à 24 avec 2+1 = 3 degrés de libertés.
round(Qtests(arima201$residuals,24,fitdf=3),3)
#Toutes les p-valeurs des tests de Ljung-Box réalisés pour des lags allant de 4 à 24 sont supérieures à 0.05, donc l'hypothèse de décorrélation des résidus n'est pas rejetée.

#Test de Student (applicable car résidus blancs)
signif(arima201)
#On observe que les coefficients du modèle ne sont pas significatifs, indice d'une surparamétrisation. On va simplifier le modèle en appliquant la méthodologie du rasoir d'Occam.

#On essaye l'ordre 1 pour une moyenne mobile.
arima001 <- arima(dfm,c(0,0,1))
#On teste la blancheur des résidus du modèle avec des tests de Ljung-Box.
round(Qtests(arima001$residuals,50,fitdf=1),3)
#Toutes les p-valeurs des tests de Ljung-Box réalisés pour des lags allant de 2 à 24 sont supérieures à 0.05, donc l'hypothèse de décorrélation des résidus, ie de blancheur des résidus, n'est pas rejetée. Donc le modèle est valide.
#Test de Student (applicable car les résidus du modèle sont bien blancs.)
signif(arima001)
#La p-valeur pour le coefficient d'ordre 1 de moyenne mobile est inférieure à 0.01, donc le coefficient d'ordre 1 de moyenne mobile est significatif et le modèle est bien ajusté.
ma1 <- arima001 #On stocke ce modèle.

#On essaye l'ordre 1 pour un processus autoregressif.
arima100 <- arima(dfm,c(1,0,0))
#On teste la blancheur des résidus du modèle avec des tests de Ljung-Box.
round(Qtests(arima100$residuals,50,fitdf=1),3)
#La p-valeur du test au lag 2 est nulle ou presque nulle, donc les résidus du modèle ne sont pas blancs. Il reste des corrélations à expliquer et le modèle est invalide.

#On essaye l'ordre 2 pour une moyenne mobile.
arima200 <- arima(dfm,c(2,0,0))
#On teste la blancheur des résidus du modèle avec des tests de Ljung-Box.
round(Qtests(arima200$residuals,50,fitdf=2),3)
#Toutes les p-valeurs des tests de Ljung-Box réalisés pour des lags allant de 3 à 24 sont supérieures à 0.05, donc l'hypothèse de décorrélation des résidus, ie de blancheur des résidus, n'est pas rejetée. Donc le modèle est valide.
#Test de Student (applicable car les résidus du modèle sont bien blancs.)
signif(arima200)
#L'ensemble des p-valeurs des coefficients sont inférieure à 0.01, donc les coefficients sont significatifs au niveau 1% et le modèle est bien ajusté.
ar2 <- arima200 #On stocke ce modèle.

#Partie pouvant être passée, totalement inutile et visant seulement à confirmer une intuition.
#Le fait que le modèle AR2 et surtout le modèle MA1 soit valide nous fait présumer que tester d'autres moddèles avec davantage de paramètres pourra nous donner des modèles certes valides mais probablement mal-ajustés.
#On vérifie cette intuition.
arima101 <- arima(dfm,c(1,0,1)) 
round(Qtests(arima101$residuals,50,fitdf=2),3) #Comme prévu, le modèle est valide
signif(arima101) # Comme présumé, le coefficient AR1 est non-significatif, le modèle MA1 étant déjà valide et ajusté.

arima201 <- arima(dfm,c(2,0,1))
round(Qtests(arima201$residuals,50,fitdf=3),3) #Comme prévu, le modèle est valide.
signif(arima201) #Comme présumé, le coefficient MA1 est non-significatif, le modèle AR2 étant déjà valide et ajustée.
#Fin de la partie pouvant être passée.

#Bilan : les modèles valides et ajustés sont : arima001, arima200

#Calcul des AICs et BICs
models <- c("ar2","ma1"); names(models) <- models
apply(as.matrix(models),1, function(m) c("AIC"=AIC(get(m)), "BIC"=BIC(get(m))))




#OPTIONNEL 
      

#Affichage des ACF et PACF
par(mfrow=c(1,2))
pacf(arima001$residuals)
acf(arima001$residuals)
#Affichage de la série stationnarisée et de la série des résidus du modèle
par(mfrow=c(1,2))
plot(arima001$residuals)
plot(dfm)
#Comparaison de l'amplitude des variations
v = dfm**2
s = (arima001$residuals)**2
sqrt(mean(v))
sqrt(mean(s))


#Affichage des ACF et PACF
par(mfrow=c(1,2))
pacf(arima200$residuals)
acf(arima200$residuals)
#Affichage de la série stationnarisée et de la série des résidus du modèle
par(mfrow=c(1,2))
plot(arima200$residuals)
plot(dfm)
#Comparaison de l'amplitude des variations
v = dfm**2
s = (arima200$residuals)**2
sqrt(mean(v))
sqrt(mean(s))
#FIN DE L'OPTIONNEL

      
#On observe que le modèle MA(1) minimise à la fois les critères AIC et BIC. 

      
#Donc on choisit un modèle MA(1) = ARMA(0,1) pour modéliser la série différenciée.










      
#Bilan : on modélise la série de départ par un ARIMA(0,1,1)
arima011 <- arima(fm,c(0,1,1))
final_model <- arima011 #Renommage pour des raisons de lisibilité
      












arimafit <- function(estim){
  adjust <- round(signif(estim),3)
  pvals <- Qtests(estim$residuals,24,length(estim$coef)-1)
  pvals <- matrix(apply(matrix(1:24,nrow=6),2,function(c) round(pvals[c,],3)),nrow=6)
  colnames(pvals) <- rep(c("lag", "pval"),4)
  cat("tests de nullite des coefficients :\n")
  print(adjust)
  cat("\n tests d'absence d'autocorrelation des residus : \n")
  print(pvals)
}

estim <- arima(dfm,c(3,0,1)); arimafit(estim)
# valide mal ajusté


#On affiche les ACF et PACF
par(mfrow=c(2,1))
plot(dfm)
plot(ar2$residuals) #les résidus sont trop mims

acf(ar2$residuals, lag = 50)


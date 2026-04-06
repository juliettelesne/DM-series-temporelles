require(zoo) #format de serie temporelle pratique et facile d'utilisation (mais plus volumineux)
require(tseries) #diverses fonctions sur les series temporelles

install.packages("fUnitRoots")
library(fUnitRoots)
#PREMIERE ETAPE : RENDRE LA SERIE EXPLOITABLE (13)
#DEUXIEME ETAPE : RENDRE LA SERIE TEMPORELLE STATIONNAIRE (43)
#TROISIEME ETAPE : SOUS L'HYPOTHESE D'UNE SERIE STATIONNAIRE AUTOUR D'UNE TENDANCE LINEAIRE (72)
#QUATRIEME ETAPE : SOUS L'HYPOTHESE D'UNE RACINE UNITAIRE AVEC TENDANCE (342)
#
#Dans tous les cas, les résidus sont énormes : j'ai l'impression que c'est pas terrible terrible :(

#PREMIERE ETAPE : RENDRE LA SERIE EXPLOITABLE

getwd() #affiche le wd
list.files() #liste les elements du wd

datafile <- "valeurs_mensuelles.csv" #définit le fichier de données

data <- read.csv(datafile,sep=";") #importe un fichier .csv dans un objet de classe data.frame

fm.source <- zoo(data[[2]]) #convertit le deuxième élément de data en serie temporelle de type "zoo"
print(fm.source) #observe que les premières valeurs sont les plus récentes

fm.source_rev <- rev(fm.source) #inverse l'ordre des valeurs pour obtenir un ordre chronologique
print(fm.source_rev) #observe que les 3 dernières valeurs du vecteur ne sont pas des valeurs d'indices

T <- length(fm.source_rev) #prépare la suppression de ces trois dernières valeurs
fm.raw <- fm.source_rev[1:(T-3)] #supprime les 3 dernières valeurs  

str(fm.raw) #vérifie si les valeurs rentrées sont des caractères ou des nombres : ce sont des caractères
fm.num <- as.numeric(fm.raw) #convertit les caractères en nombres
str(fm.num) #vérifie qu'il s'agit bien du format numérique, nécessaire au traitement
class(fm.num) #est sorti du format "zoo"

fm <- zoo(fm.num) #reconvertit en série temporelle de type "zoo"
class(fm) #vérifie que l'on fait à nouveau face à une série temporelle de type "zoo"
print(fm) #pas d'anomalie observée




#DEUXIEME ETAPE : RENDRE LA SERIE TEMPORELLE STATIONNAIRE
#Indice CVS-CJO de la production industrielle (base 100 en 2021)
# - Métallurgie des autres métaux non ferreux (NAF rév. 2, niveau classe, poste 24.45) 

#Premier tracé de la série
plot(fm)

abline(v = 181, col = "red", lwd = 2, lty = 2) #rupture toruvée par Zivot Andrews
#On observe deux pics, mais a priori, les oscillations autour de la tendance sont d'amplitude relativement stable
#Le fait que des grands chocs n'altèrent pas durablement la valeur de la série
#laisse supposer que l'on est plus probablement dans le cas d'une série stationnaire
#avec peut-être une légère tendance
#que dans le cas d'une série à racine unitaire.

#On affiche les ACF et PACF
par(mfrow=c(1,2))
pacf(fm) #faible valeur inférieure à 0.4, donc il n'y a a priori pas de racine unitaire
acf(fm) # décroît lentement : signe peut-être d'une légère tendance ou d'une composante moyenne mobile




#Test rejetant l'hypothèse d'une racine unitaire :
#Début de crime : suppose de travailler sur un AR(1)
adfTest(fm, lags=0, type = "ct") # p-value <= 0.01 : on rejette au niveau 1% l'hypothèse de non-stationarité
#Fin de crime

#Général : le test de Perron-Phillips est valable dans un cas très général
pp.test(fm)

library(urca)
test_pp <- ur.pp(fm, type = "Z-tau", model = "constant", lags = "short")
summary(test_pp)
#Pour ce type de test : p-value <= 0.01 : on rejette au niveau 1% l'hypothèse de non-stationarité
#PP est robuste aux ruptures, pas KPSS

#KPSS
kpss.test(fm, null = "Level") # p-value = 0.02162 : on rejette au niveau 5% l'hypothèse de stationnarité
kpss.test(fm, null = "Trend") #p-value <= 0.01 : on rejette niveau 1% l'hypothèse de stationnarité avec tendance linéaire

#Test de l'existence d'une rupture
library(urca)
summary(ur.za(fm)) #potential break point at position : 181

#ON EST BLOQUES 


#TROISIEME ETAPE : SOUS L'HYPOTHESE D'UNE SERIE STATIONNAIRE AUTOUR D'UNE TENDANCE LINEAIRE

#On détermine la tendance :
#régression linéaire de fm sur le temps
tf <- time(fm)
model <- lm(fm ~ tf) 

#Série temporelle des résidus du modèle
resf <- model$residuals

#Tracé
par(mfrow=c(1,2))
plot(fm, col = "blue", lwd = 2)
abline(model, col = "red", lwd = 2) #On a peut-être une légère tendance croissante
plot(resf) #Les résidus sont approximativement centrés de variance approximativement stable avec quelques outliers

#Etude des ACF et PACF des résidus
par(mfrow=c(2,2))
acf(resf)
pacf(resf)
acf(fm)
pacf(fm)
#enlever la tendance ne change pas la forme des ACF et PACF

#Test de stationnarité de la série détrendée

#kpss toujours valable
kpss.test(resf) #p-valeur > 0.1 : on ne rejette plus l'hypothèse de stationnarité

#Début de crime
#adf suppose que l'on fait face à un AR(1) => pas fiable dans notre cas (on utilise 'nc' car on a corrigé la tendance)
adfTest(resf, lags=0, type = "nc") #p-value <= 0.01 : on rejette au niveau 1% l'hypothèse de non-stationarité
#Fin de crime

#Perron-Philips : très général
pp.test(resf) #p-value <= 0.01 : on rejette au niveau 1% l'hypothèse de non-stationarité

#Bilan : une fois la série détrendée, il n'existe pas de test allant à l'encontre de la stationnarité 



#Modélisons resf
#Ses PACF sont très faibles dès l'ordre 2 et a priori nulles après l'ordre 4 : on peut essayer un ordre d'AR(4) : celui-ci correspondrait à des ACF décroissante petit à petit
#Ses ACF diminuent lentement : on peut essayer un ordre MA de 2, au-delà duquel l'ACF réalise un premier saut vers le bas

arima402 <- arima(resf,c(4,0,2))

plot(arima402$residuals) # les résidus sont bien centrés de variance approximativement constante avec quelques outliers

#Les PACF et ACF correspondent bien à un bruit : PACF nulle, ACF nulle dès l'ordre 2
par(mfrow=c(1,2))
acf(arima402$residuals)
pacf(arima402$residuals)

#Test de blancheur des résidus du modèle

#on teste la nullité jointe des covariables du lag 1 au lag L entre les résidus
#grande p-valeur : on ne rejette pas H_0 : décorrélation du bruit    
#on fait croître la taille de la fenêtre pour éviter de ne pas prendre en compte un pic, ou au contraire de le gommer (pour se couvrir de la sensibilité).
# - petite fenètre : les pics à court terme ont un gros impact
# - grande fenêtre : perte de puissance de test (un gros pic est gommé)
#on choisit le lag max: on supose qu'au bout de deux ans, il ne devrait plus y avoir de corrélation dans les résidus
Qtests <- function(series, k, fitdf=0) {
  pvals <- apply(matrix(1:k), 1, FUN=function(l) {
    pval <- if (l<=fitdf) NA else Box.test(series, lag=l, type="Ljung-Box", fitdf=fitdf)$p.value
    return(c("lag"=l,"pval"=pval))
  })
  
  #fitdf : nb deg liberté
  return(t(pvals))
}
Qtests(arima402$residuals, 30, 6) #fitdf=6 : on a 6 paramètres libres dans l'estimation du modèle ARMA(4,2)
round(Qtests(arima402$residuals,30,fitdf=6),3)

#On ne rejette pas l'hypothèse de blancheur des résidus 

#: le modèle est ajusté

#On vérifie la significativité des coefficients :
signif <- function(estim){ #fonction de test des significations individuelles des coefficients
  coef <- estim$coef
  se <- sqrt(diag(estim$var.coef))
  t <- coef/se
  pval <- (1-pnorm(abs(t)))*2
  return(rbind(coef,se,pval))
}

signif(arima402) # Les coef ar3 et ar4 semblent peu significatifs : on essaye un ARMA(3,2)






arima302 <- arima(resf,c(3,0,2))

plot(arima302$residuals) # les résidus sont bien centrés de variance approximativement constante avec quelques outliers

#Les PACF et ACF correspondent bien à un bruit : PACF nulle, ACF nulle dès l'ordre 2
par(mfrow=c(1,2))
acf(arima302$residuals)
pacf(arima302$residuals)

#Test de blancheur des résidus du modèle

#on teste la nullité jointe des covariables du lag 1 au lag L entre les résidus
#grande p-valeur : on ne rejette pas H_0 : décorrélation du bruit    
#on fait croître la taille de la fenêtre pour éviter de ne pas prendre en compte un pic, ou au contraire de le gommer (pour se couvrir de la sensibilité).
# - petite fenètre : les pics à court terme ont un gros impact
# - grande fenêtre : perte de puissance de test (un gros pic est gommé)
#on choisit le lag max: on supose qu'au bout de deux ans, il ne devrait plus y avoir de corrélation dans les résidus
# probablement visible au bout de 2 cycle si pb
Qtests(arima302$residuals, 30, 5) #fitdf=5 : on a 5 paramètres libres dans l'estimation du modèle ARMA(4,2)
round(Qtests(arima302$residuals,30,fitdf=5),3)

#On ne rejette pas l'hypothèse de blancheur des résidus 

#: le modèle est ajusté

#On vérifie la significativité des coefficients :
signif(arima302)
# Le coef ar3 semble peu significatif : on peut essayer un ARMA(2,2)
# Le coef ma2 semble peu significatif : on peut essayer un ARMA(3,1)





arima202 <- arima(resf,c(2,0,2))

plot(arima202$residuals) # les résidus sont bien centrés de variance approximativement constante avec quelques outliers

#Les PACF et ACF correspondent bien à un bruit : PACF nulle, ACF nulle dès l'ordre 2
par(mfrow=c(1,2))
acf(arima202$residuals)
pacf(arima202$residuals)

#Test de blancheur des résidus du modèle

#on teste la nullité jointe des covariables du lag 1 au lag L entre les résidus
#grande p-valeur : on ne rejette pas H_0 : décorrélation du bruit    
#on fait croître la taille de la fenêtre pour éviter de ne pas prendre en compte un pic, ou au contraire de le gommer (pour se couvrir de la sensibilité).
# - petite fenètre : les pics à court terme ont un gros impact
# - grande fenêtre : perte de puissance de test (un gros pic est gommé)
#on choisit le lag max: on supose qu'au bout de deux ans, il ne devrait plus y avoir de corrélation dans les résidus
# probablement visible au bout de 2 cycle si pb
Qtests(arima202$residuals, 30, 4) #fitdf=5 : on a 4 paramètres libres dans l'estimation du modèle ARMA(4,2)
round(Qtests(arima202$residuals,30,fitdf=4),3)

#On ne rejette pas l'hypothèse de blancheur des résidus 

#: le modèle est ajusté

#On vérifie la significativité des coefficients :
signif(arima202) # Les coefficients au dernier ordre sont significatifs : on ne peut pas simplifier le modèle.




arima301 <- arima(resf,c(3,0,1))

plot(arima301$residuals) # les résidus sont bien centrés de variance approximativement constante avec quelques outliers

#Les PACF et ACF correspondent bien à un bruit : PACF nulle, ACF nulle dès l'ordre 2
par(mfrow=c(1,2))
acf(arima301$residuals)
pacf(arima301$residuals)

#Test de blancheur des résidus du modèle

#on teste la nullité jointe des covariables du lag 1 au lag L entre les résidus
#grande p-valeur : on ne rejette pas H_0 : décorrélation du bruit    
#on fait croître la taille de la fenêtre pour éviter de ne pas prendre en compte un pic, ou au contraire de le gommer (pour se couvrir de la sensibilité).
# - petite fenètre : les pics à court terme ont un gros impact
# - grande fenêtre : perte de puissance de test (un gros pic est gommé)
#on choisit le lag max: on supose qu'au bout de deux ans, il ne devrait plus y avoir de corrélation dans les résidus
Qtests(arima301$residuals, 30, 4) #fitdf=5 : on a 4 paramètres libres dans l'estimation du modèle ARMA(4,2)
round(Qtests(arima301$residuals,30,fitdf=4),3)

#On ne rejette pas l'hypothèse de blancheur des résidus 

#: le modèle est ajusté

#On vérifie la significativité des coefficients :
signif(arima301) #Le coefficient MA est significatif, mais pas le coefficient ar3. On peut essayer un ARMA(2,1)


arima201 <- arima(resf,c(2,0,1))

plot(arima201$residuals) # les résidus sont bien centrés de variance approximativement constante avec quelques outliers

#Les PACF et ACF correspondent bien à un bruit : PACF nulle, ACF nulle dès l'ordre 2
par(mfrow=c(1,2))
acf(arima201$residuals)
pacf(arima201$residuals)

#Test de blancheur des résidus du modèle

#on teste la nullité jointe des covariables du lag 1 au lag L entre les résidus
#grande p-valeur : on ne rejette pas H_0 : décorrélation du bruit    
#on fait croître la taille de la fenêtre pour éviter de ne pas prendre en compte un pic, ou au contraire de le gommer (pour se couvrir de la sensibilité).
# - petite fenètre : les pics à court terme ont un gros impact
# - grande fenêtre : perte de puissance de test (un gros pic est gommé)
#on choisit le lag max: on supose qu'au bout de deux ans, il ne devrait plus y avoir de corrélation dans les résidus
Qtests(arima201$residuals, 30, 3) #fitdf=5 : on a 3 paramètres libres dans l'estimation du modèle ARMA(4,2)
round(Qtests(arima201$residuals,30,fitdf=3),3)

#On ne rejette pas l'hypothèse de blancheur des résidus 

#: le modèle est ajusté

#On vérifie la significativité des coefficients :
signif(arima201) #Le coefficient MA est significatif, mais pas le coefficient ar2. On peut essayer un ARMA(1,1)



arima101 <- arima(resf,c(1,0,1))

plot(arima101$residuals) # les résidus sont bien centrés de variance approximativement constante avec quelques outliers

#Les PACF et ACF correspondent bien à un bruit : PACF nulle, ACF nulle dès l'ordre 2
par(mfrow=c(1,2))
acf(arima101$residuals)
pacf(arima101$residuals)

#Test de blancheur des résidus du modèle

#on teste la nullité jointe des covariables du lag 1 au lag L entre les résidus
#grande p-valeur : on ne rejette pas H_0 : décorrélation du bruit    
#on fait croître la taille de la fenêtre pour éviter de ne pas prendre en compte un pic, ou au contraire de le gommer (pour se couvrir de la sensibilité).
# - petite fenètre : les pics à court terme ont un gros impact
# - grande fenêtre : perte de puissance de test (un gros pic est gommé)
#on choisit le lag max: on supose qu'au bout de deux ans, il ne devrait plus y avoir de corrélation dans les résidus
Qtests(arima101$residuals, 30, 2) #fitdf=5 : on a 2 paramètres libres dans l'estimation du modèle ARMA(4,2)
round(Qtests(arima101$residuals,30,fitdf=2),3)

#On ne rejette pas l'hypothèse de blancheur des résidus 

#: le modèle est ajusté

#On vérifie la significativité des coefficients :
signif(arima101) #Le modèle ne peut pas être simplifié.






#L'ajout de paramètres ar3 et ar4 ne fait quasiment pas changer les résidus.
par(mfrow=c(3,3))
plot(arima202$residuals)
plot(arima302$residuals)
plot(arima402$residuals)
plot(arima301$residuals)
plot(arima201$residuals)
plot(arima101$residuals)
plot(resf)
#Les résidus semblent centrés, sans tendance, sans motifs périodiques, de variance semblant relativement constante dans le temps, bien que quelques outliers importants apparaissent.
#Ces signaux signal peuvent ressembler à des bruits blancs avec quelques valeurs extrêmes.
signif(arima202)
signif(arima302)
signif(arima402)
signif(arima301)
signif(arima201)
signif(arima101)
#On a tout de meme des coefficients AR proches de 1 et des résidus importants




#QUATRIEME ETAPE : SOUS L'HYPOTHESE D'UNE RACINE UNITAIRE AVEC TENDANCE

#Pour stationnariser, on différencie

variation <- fm-lag(fm,-1) #différencie
print(variation)
par(mfrow=c(1,1))
plot(variation) # a part outliers, semble stationnaire, évolue dans un couloir, pas de tendance aléatoire

#Tests de stationnarité de la série différenciée
kpss.test(variation) #vraiment robuste : p-valeur > 0.1 : ne rejette pas l'hypothèse de stationarité
#Début de crime
adf.test(variation) #hypothèse peut crédible à la vue des ACF et PACF : pas un AR(1) (p-valeur <= 0.01 : on rejete l'hypothèse de racine unitaire)
#Fin de crime
pp.test(variation) #p-valeur <= 0.01 : on rejete l'hypothèse de racine unitaire


#rien ne fait supposer une racine unitaire : PACF faible, pas de persistance de l'ACF
par(mfrow=c(1,2))
acf(variation) # ordre MA inférieur ou égal à 1
pacf(variation)# décroissance progressive de la PACF (exponentielle peut-être) on présume un ordre AR inférieru à 3



rw301 <- arima(variation,c(3,0,1))

plot(rw301$residuals) # les résidus sont bien centrés de variance approximativement constante avec quelques outliers

#Les PACF et ACF correspondent bien à un bruit : quasi nulle, ACF quasi nulle dès l'ordre 2
par(mfrow=c(1,2))
acf(rw301$residuals, lag = 100)
pacf(rw301$residuals, lag = 100)

#Test de blancheur des résidus du modèle

#on teste la nullité jointe des covariables du lag 1 au lag L entre les résidus
#grande p-valeur : on ne rejette pas H_0 : décorrélation du bruit    
#on fait croître la taille de la fenêtre pour éviter de ne pas prendre en compte un pic, ou au contraire de le gommer (pour se couvrir de la sensibilité).
# - petite fenètre : les pics à court terme ont un gros impact
# - grande fenêtre : perte de puissance de test (un gros pic est gommé)
#on choisit le lag max: on supose qu'au bout de deux ans, il ne devrait plus y avoir de corrélation dans les résidus
Qtests(rw301$residuals, 30, 4) #fitdf=5 : on a 4 paramètres libres dans l'estimation du modèle ARMA(4,2)
round(Qtests(rw301$residuals,30,fitdf=4),3)

#On ne rejette pas l'hypothèse de blancheur des résidus 

#: le modèle est ajusté

#On vérifie la significativité des coefficients :
signif(rw301) #Le coefficient ar3 n'est pas significatif : on peut essayer un ARMA(2,1)





rw201 <- arima(variation,c(2,0,1))

plot(rw201$residuals) # les résidus sont bien centrés de variance approximativement constante avec quelques outliers

#Les PACF et ACF correspondent bien à un bruit : quasi nulle, ACF quasi nulle dès l'ordre 2
par(mfrow=c(1,2))
acf(rw201$residuals, lag = 100)
pacf(rw201$residuals, lag = 100)

#Test de blancheur des résidus du modèle

#on teste la nullité jointe des covariables du lag 1 au lag L entre les résidus
#grande p-valeur : on ne rejette pas H_0 : décorrélation du bruit    
#on fait croître la taille de la fenêtre pour éviter de ne pas prendre en compte un pic, ou au contraire de le gommer (pour se couvrir de la sensibilité).
# - petite fenètre : les pics à court terme ont un gros impact
# - grande fenêtre : perte de puissance de test (un gros pic est gommé)
#on choisit le lag max: on supose qu'au bout de deux ans, il ne devrait plus y avoir de corrélation dans les résidus
Qtests(rw201$residuals, 30, 3) #fitdf=5 : on a 4 paramètres libres dans l'estimation du modèle ARMA(4,2)
round(Qtests(rw201$residuals,30,fitdf=3),3)

#On ne rejette pas l'hypothèse de blancheur des résidus 

#: le modèle est ajusté

#On vérifie la significativité des coefficients :
signif(rw201) #Le modèle ne peut pas être simplifié
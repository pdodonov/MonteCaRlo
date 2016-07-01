###Code written by Pavel Dodonov and used in the manuscript "Landscape-scale deforestation decreases gene flow distance of a keystone tropical palm, Euterpe Edulis Mart (Arecaceae)", by Alesandro S. Santos, Eliana Cazetta, Deborah Faria, Pavel Dodonov and Fernanda A. Gaiotto (in review).


###Functions to calculate significance of a logistic regression by unrestricted permutations
glm.signif<-function(x, y, Nperm=5000) {
	regr.orig<-glm(y~x, family=binomial)$deviance
	regr.perm<-numeric(Nperm)
	regr.perm[1]<-regr.orig
	for (i in 2:Nperm) {
		y.perm<-sample(y)
		regr.perm[i]<-glm(y.perm~x, family=binomial)$deviance
		}
	signif<-sum(regr.perm <= regr.orig) / Nperm
	return(signif)
	}

glm.signif.weights<-function(x, y, Nperm=5000, w) { #includes weights, i.e. number of observations corresponding to each proportion
	regr.orig<-glm(y~x, family=binomial, weights=w)$deviance
	regr.perm<-numeric(Nperm)
	regr.perm[1]<-regr.orig
	for (i in 2:Nperm) {
		x.perm<-sample(x)
		regr.perm[i]<-glm(y~x.perm, family=binomial, weights=w)$deviance
		}
	signif<-sum(regr.perm <= regr.orig) / Nperm
	return(signif)
	}

# Script for assessing significance of a linear model with restricted randomizations

x<-dados2$Cover #Forest cover in each landscape
y<-dados2$logGeneticFlow #Genetic flow observed for each plant (varies within landscapes)
group<-as.factor(dados2$Cover)
levels.use<-unique(group)
Nlevels<-length(levels.use)
regr.orig<-summary(lm(y~x))$r.squared
regr.perm<-numeric(5000)
regr.perm[1]<-regr.orig
for(i in 2:5000) {
	covers.rand<-sample(covers)
	x.rand<-x
	for(j in 1:Nlevels) {
		x.rand[x==levels.use[j]]=covers.rand[j]
		}
	regr.perm[i]<-summary(lm(y~x.rand))$r.squared
	print(i)
	}
hist(regr.perm)	
abline(v<-regr.orig)
signif.flow<-sum(regr.perm >= regr.orig)/5000
pred.flow<-predict(mod1)
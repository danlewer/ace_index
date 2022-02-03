library(data.table)
z <- function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)

d <- read.csv(url('https://raw.githubusercontent.com/danlewer/ace_index/main/ace_index_with_la_characteristics_3feb2022.csv'))
setDT(d)

# region reorder by ACE index

d$RGN17NM[d$RGN17NM == 'Yorkshire and The Humber'] <- 'Yorkshire & TH'
rgn <- d[, mean(ACE_rank), RGN17NM]
rgn <- rgn$RGN17NM[order(rgn$V1)]
d$RGN17NM <- factor(d$RGN17NM, rgn, rgn)

#-------------
# Scatterplots
#-------------

f_scatter <- function(x, dat = d, XX = 1, YX = T, TIT = '', logx = c(0.5, 1, 2, 5, 10, 20, 50, 100)) {
  dat$ind <- dat[,x, with = F]
  m <- lm(ACE_rank ~ ind, dat)
  R2 = format(round(summary(m)$r.squared, 3), nsmall = 3)
  TITLE <- paste0(TIT, ' (R2=', R2, ')')
  plot(dat$ind, dat$ACE_rank, axes = F, xlab = NA, ylab = NA, pch = 4, panel.first={points(0, 0, pch=16, cex=1e6, col="grey94")})
  title(main = TITLE, line = 0.5)
  if (XX == 1) axis(1, seq(0, 300, 50), las = 2)
  if (XX == 2) axis(1, seq(-324, 24, 50), seq(0, 300, 50), las = 2)
  if (XX == 3) axis(1, at = log(logx), logx, las = 2)
  if (YX) axis(2, las = 2)
  box(which = 'plot')
  abline(m, col = 'red', lwd = 2)
}

d$logDens <- log(d$density_keep)
d$rank_s7525 <- rank(d$s7525)

ds <- c("imd_rank_av_score", "income_rank_av_score", "employment_rank_av_score", "education_rank_av_score", "health_rank_av_score", "crime_rank_av_score", "barriers_rank_av_score", "living_rank_av_score", "idaci_rank_av_score", "rank_s7525", "logDens")
ds2 <- ds[1:9]
d[, (ds2) := lapply(.SD, function(x) 325 - rank(x)), .SDcols = ds2]
dst <- data.table(ds = ds, titles = c('IMD: index', 'IMD: income', 'IMD: employment', 'IMD: education', 'IMD: health', 'IMD: crime', 'IMD: barriers', 'IMD: living', 'IMD: child poverty','Inequality', 'Population/hectare'))
fr <- function(x) summary(lm(as.formula(paste0('ACE_rank~', x)), d))$r.squared
ds_r2 <- sapply(ds, fr)
ds <- ds[order(ds_r2, decreasing = T)]
titles <- dst$titles[match(ds, dst$ds)]

par(mfrow = c(4, 3), mar = c(3, 1, 2, 1), cex = 0.5, oma = c(5, 5, 0, 0), xpd = F)
mapply(f_scatter, x = ds, YX = rep(c(T, F, F), 4)[1:11], XX = c(2, 2, 2, 2, 2, 2, 2, 3, 2, 1, 1), TIT = titles)

m <- lm(ACE_rank ~ RGN17NM, d)
R2 = format(round(summary(m)$r.squared, 2), nsmall = 2)
y <- boxplot(ACE_rank ~ RGN17NM, d, axes = F, col = 'grey85', main = paste0('Region', ' (R2=', R2, ')'), whisklty = 1, medlwd = 2)
box(which = 'plot')
axis(1, 1:9, y$names, las = 2)

mtext('Rank of indicator (324 = worst)', side = 1, outer = T, line = 2, cex = 0.7)
mtext('ACE index (324 = highest occurence of ACEs)', side = 2, outer = T, line = 3, cex = 0.7)

#-----------
# regression
#-----------

# density is the population density per hectare in 2011.
# s7525 is a measure of within-local-authority inequality. It was calculated as "the ratio between the mean IDACI scores for the bottom and top quartile of neighbourhoods (Lower Super Output Areas, which are small areas with an average population of 1500 are) within a local authority. This measure ranged from 1.3 to 4.8. Population density was defined as the number of residents per hectare in 2011 (the most recent national census year)."

d[, idaci2 := -z(idaci_rank_av_score)]
d[, density2 := z(log(density_keep))]
d[, s75252 := z(rank(s7525))]

# example regression formula

m <- lm(ACE_rank ~ idaci2 + density2 + s75252 + RGN17NM, data = d)
summary(m) # main results
plot(m) # model diagnostics

# create table 1

vars <- c('idaci2', 'density2', 's75252', 'RGN17NM')

ind <- function(x) {
  form <- paste0('ACE_rank~', x)
  m <- lm(form, d)
  cbind(m$coef, confint(m), round(summary(m)$coef[,4], 5))
}
unadj <- do.call('rbind', lapply(vars, ind))

af <- as.formula(paste0('ACE_rank~', paste0(vars, collapse = '+')))
adj <- lm(af, d)
#ols_vif_tol(adj) # not real collinearity
adj <- cbind(adj$coef, confint(adj), round(summary(adj)$coef[,4], 5))

format_mod <- function(mr) {
  mr <- mr[row.names(mr) != '(Intercept)',]
  mr <- data.frame(mr)
  mr[, 1:3] <- format(round(mr[,1:3], 2), nsmall = 2)
  mr$res <- paste0(mr$V1, ' (', mr$X2.5.., ', ', mr$X97.5.., ') ')
  mr$p <- ifelse(mr$V4 < 0.001, 'p<0.001', paste0('p=', format(round(mr$V4, 3), nsmall = 3)))
  mr$res <- paste0(mr$res, mr$p)
  mr$res <- stri_replace_all_fixed(mr$res, '( ', '(')
  mr$res <- stri_replace_all_fixed(mr$res, '  ', ' ')
  mr$res <- stri_replace_all_fixed(mr$res, '  ', ' ')
  names(mr) <- LETTERS[1:6]
  mr[, c('A', 'B', 'C', 'D', 'F')] <- NULL
  mr
}

coefs <- cbind(format_mod(unadj), format_mod(adj))
#write.csv(coefs, 'Table1.csv')

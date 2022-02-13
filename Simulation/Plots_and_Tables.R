library(ggplot2)
library(plyr)
path = paste0(dirname(rstudioapi::getActiveDocumentContext()$path), '/simulation results/')

errors = c('N(0,1)', 't_2', 'Par(4,2)','Burr(1,2,1)')
theme_set(theme_bw())
### 4.1 Distributed AHR
# Line plot
load(paste0(path,'outputLowHetqua.Rdata'))
# for each of the 4 errors, for each of the 4 estimators, for each of the 6 
# different m, calculate mean l2 error and plot
for (i in 1:dim(output)[2]){
  df = data.frame(matrix(rep(0,3*length(m_list)*dim(output)[3]), ncol = 3))
  # df = data.frame(matrix(rep(0, length(m_list)*(dim(output)[3] + 1)), length(m_list)))
  # df$X1 = m_list
  colnames(df) = c('n_machine','estimation_error', 'estimator')
  for (j in 1:dim(output)[3]){
    for (k in 1:dim(output)[1]){
      betas = output[k,i,j,,1,]
      df[k + (j-1)*length(m_list),] = c(m_list[k], mean(apply(betas,1, function(x) norm(x - c(beta0.true, beta.true), type = '2'))), j)
    }
  }
  temp.plot = ggplot(data = df, aes(x = n_machine, y = estimation_error, group  = estimator)) + 
    geom_line(linetype = df$estimator, color = mapvalues(df$estimator, from = 1:5, to = c('black', 'blue1','steelblue','red', 'green4')), size = 0.8) +
    geom_point(shape = df$estimator + 14,color = mapvalues(df$estimator, from = 1:5, to = c('black', 'blue1','steelblue','red','green4')),size = 3.5) + 
    theme(panel.grid = element_blank(), legend.position = 'none', text = element_text(size=20), axis.text.x= element_text(size=20), axis.text.y= element_text(size=20),axis.ticks.length=unit(.25, "cm"),plot.title = element_text(hjust = 0.5)) + 
    labs(title = ifelse(i==2, expression(t[2]),errors[i]), x = "Number of Machines", y = "Estimation Error")
  ggsave(paste0(path, 'LowHetLine' , i,'.pdf'), plot = temp.plot, units="in", width=5, height=5, dpi=1200)
}



# Box plot

load(paste0(path, 'outputLowHetQua.Rdata'))
# for each of the 4 errors, for each of the 2 estimators (3: distributed-OLS and 4:Dist), for each of the 6 
# different m, calculate coverage and width, then boxplot them
for (i in 1:dim(output)[2]){
  df = data.frame(matrix(rep(0,3*length(m_list)*2*total_runs), ncol = 3))
  colnames(df) = c('n_machine','estimator', 'error')
  for (j in 3:4){
    for (k in 1:dim(output)[1]){
      betas = output[k,i,j,,1,]
      error = apply(betas,1, function(x) norm(x - c(beta0.true, beta.true), type = '2'))
      df[(1 + (k-1)*2*total_runs + (j-3)*total_runs):((k-1)*2*total_runs + (j-2)*total_runs),] = cbind(rep(m_list[k],total_runs), rep(j, total_runs), error)
    }
  }
  temp.plot = ggplot(data = df, aes(x = factor(n_machine), y = error, fill = factor(estimator))) +
    geom_boxplot() + 
    theme(panel.grid = element_blank(), legend.position = 'none', text = element_text(size=20), axis.text.x= element_text(size=20), axis.text.y= element_text(size=20),axis.ticks.length=unit(.25, "cm"),plot.title = element_text(hjust = 0.5)) + 
    labs(title = ifelse(i==2, expression(t[2]),errors[i]), x = "Number of Machines", y = "Estimation Error")
  
  ggsave(paste0(path, 'LowHetBox' , i,'.pdf'), plot = temp.plot, units="in", width=5, height=5, dpi=1200)
}






### 4.2 Distributed Confidence Construction



# Coverage plot   and    width plot
load(paste0(path, 'outputLowHetQua.Rdata'))
# for each of the 4 errors, for each of the 2 estimators (3: distributed-OLS and 4:Dist), for each of the 6 
# different m, calculate coverage and width, then boxplot them
for (i in 1:dim(output)[2]){
  nCoeff = dim(output)[6]-1
  df = data.frame(matrix(rep(0,4*length(m_list)*2*nCoeff), ncol = 4))
  colnames(df) = c('n_machine','estimator', 'cover_prob', 'width')
  for (j in 3:4){
    for (k in 1:dim(output)[1]){
      CIs = output[k,i,j,,2:3,-1]
      coverage = colMeans(sapply(1:nCoeff, function(x) CIs[,1,x] < beta.true[x] & CIs[,2,x] > beta.true[x]))
      width = colMeans(CIs[,2,] - CIs[,1,])
      df[(1 + (k-1)*2*nCoeff + (j-3)*nCoeff):((k-1)*2*nCoeff + (j-2)*nCoeff),] = cbind(rep(m_list[k],nCoeff), rep(j, nCoeff), coverage, width)
    }
  }
  
  temp.plot = ggplot(data = df, aes(x = factor(n_machine), y = cover_prob, fill = factor(estimator))) +
    geom_boxplot() + 
    theme(panel.grid = element_blank(), legend.position = 'none', text = element_text(size=20), axis.text.x= element_text(size=20), axis.text.y= element_text(size=20),axis.ticks.length=unit(.25, "cm"),plot.title = element_text(hjust = 0.5)) + 
    labs(title = ifelse(i==2, expression(t[2]),errors[i]), x = "Number of Machines", y = "Coverage")
  
  
  ggsave(paste0(path, 'CoverageHet' , i,'.pdf'), plot = temp.plot, units="in", width=5, height=5, dpi=1200)
  
  temp.plot = ggplot(data = df, aes(x = factor(n_machine), y = width, fill = factor(estimator))) +
    geom_boxplot() +     
    theme(panel.grid = element_blank(), legend.position = 'none', text = element_text(size=20), axis.text.x= element_text(size=20), axis.text.y= element_text(size=20),axis.ticks.length=unit(.25, "cm"),plot.title = element_text(hjust = 0.5)) + 
    labs(title = ifelse(i==2, expression(t[2]),errors[i]), x = "Number of Machines", y = "Width")
  
  
  ggsave(paste0(path, 'WidthHet' , i,'.pdf'), plot = temp.plot, units="in", width=5, height=5, dpi=1200)
}




# Table
# rows = 5 different m * 2 estimators
# columns = 4 errors * [coverage mean (sd), width mean (sd)] (first avg over simul runs, then over coeff)
load(paste0(path, 'outputLowHetQua.Rdata'))
cover.table = matrix(0,10,8)
names(cover.table) = list(c(''), c())
for (i in 1:dim(output)[2]){
  nCoeff = dim(output)[6]-1
  for (j in 3:4){
    for (k in 1:dim(output)[1]){
      CIs = output[k,i,j,,2:3,-1]
      coverage = colMeans(sapply(1:nCoeff, function(x) CIs[,1,x] < beta.true[x] & CIs[,2,x] > beta.true[x]))
      width = colMeans(CIs[,2,] - CIs[,1,])
      width.sd = apply(CIs[,2,] - CIs[,1,], 2, sd)
      entry1 = paste0(signif(mean(coverage),2), '(',signif(sd(coverage),2), ')')
      entry2 = paste0(signif(mean(width),2), '(',signif(mean(width.sd),2), ')') 
      cover.table[j-2 + (k-1)*2, ((i-1)*2 + 1):((i-1)*2 + 2)] = c(entry1, entry2)
    }
  }
}

# latex outputing
preText = c('\\multirow{2}{*}{m=50}  & Dist-OLS &',' & Dist-AHR &',
            '\\multirow{2}{*}{m=100}  & Dist-OLS &',' & Dist-AHR &',
            '\\multirow{2}{*}{m=200}  & Dist-OLS &',' & Dist-AHR &',
            '\\multirow{2}{*}{m=300}  & Dist-OLS &',' & Dist-AHR &',
            '\\multirow{2}{*}{m=400}  & Dist-OLS &',' & Dist-AHR &')
for (i in 1:nrow(cover.table)){
  string = paste(cover.table[i,], '&', collapse = ' ')
  cat(paste0(preText[i], substr(string,1,nchar(string) -1),'\\\\\n'))
}



### 4.3 Distributed regularized AHR


# line plot
try(dev.off(), silent = T)
load(paste0(path, 'outputHighHetQua.Rdata'))
# for each of the 4 errors, for each of the 4 estimators, for each of the 6 
# different m, calculate mean l2 error and plot
for (i in 1:dim(output)[2]){
  df = data.frame(matrix(rep(0,3*length(m_list)*dim(output)[3]), ncol = 3))
  # df = data.frame(matrix(rep(0, length(m_list)*(dim(output)[3] + 1)), length(m_list)))
  # df$X1 = m_list
  colnames(df) = c('n_machine','estimation_error', 'estimator')
  for (j in 1:dim(output)[3]){
    for (k in 1:dim(output)[1]){
      betas = output[k,i,j,,]
      df[k + (j-1)*length(m_list),] = c(m_list[k], mean(apply(betas,1, function(x) norm(x - c(beta0.true, beta.true), type = '2'))), j)
    }
  }
  temp.plot = ggplot(data = df, aes(x = n_machine, y = estimation_error, group  = estimator)) + 
    geom_line(linetype = df$estimator, color = mapvalues(df$estimator, from = 1:4, to = c('black', 'blue1','steelblue','red')), size = 0.8) +
    geom_point(shape = df$estimator + 14,color = mapvalues(df$estimator, from = 1:4, to = c('black', 'blue1','steelblue','red')),size = 3.5) + 
    theme(panel.grid = element_blank(), legend.position = 'none', text = element_text(size=20), axis.text.x= element_text(size=20), axis.text.y= element_text(size=20),axis.ticks.length=unit(.25, "cm"),plot.title = element_text(hjust = 0.5)) + 
    labs(title = ifelse(i==2, expression(t[2]),errors[i]), x = "Number of Machines", y = "Estimation Error")
  ggsave(paste0(path, 'HighHetLine' , i,'.pdf'), plot = temp.plot, units="in", width=5, height=5, dpi=1200)
}





# Box plot

load(paste0(path, 'outputHighHetQua.Rdata'))
# for each of the 4 errors, for each of the 2 estimators (3: distributed-OLS and 4:Dist), for each of the 6 
# different m, calculate coverage and width, then boxplot them
for (i in 1:dim(output)[2]){
  df = data.frame(matrix(rep(0,3*length(m_list)*2*total_runs), ncol = 3))
  colnames(df) = c('n_machine','estimator', 'error')
  for (j in 3:4){
    for (k in 1:dim(output)[1]){
      betas = output[k,i,j,,]
      error = apply(betas,1, function(x) norm(x - c(beta0.true, beta.true), type = '2'))
      df[(1 + (k-1)*2*total_runs + (j-3)*total_runs):((k-1)*2*total_runs + (j-2)*total_runs),] = cbind(rep(m_list[k],total_runs), rep(j, total_runs), error)
    }
  }
  temp.plot = ggplot(data = df, aes(x = factor(n_machine), y = error, fill = factor(estimator))) +
    geom_boxplot() + 
    theme(panel.grid = element_blank(), legend.position = 'none', text = element_text(size=20), axis.text.x= element_text(size=20), axis.text.y= element_text(size=20),axis.ticks.length=unit(.25, "cm"),plot.title = element_text(hjust = 0.5)) + 
    labs(title = ifelse(i==2, expression(t[2]),errors[i]), x = "Number of Machines", y = "Estimation Error")
  
  ggsave(paste0(path, 'HighHetBox' , i,'.pdf'), plot = temp.plot, units="in", width=5, height=5, dpi=1200)
}




# create line legend
index = 2
colors = c('black', 'blue1','steelblue','red','green4')
df0 = data.frame(x = 1:12, y = rep(5,12))
ggplot(data =df0, aes(x= x, y=y)) +
  geom_line(linetype = index, color = colors[index], size = 0.8) +
  geom_point(shape = index + 14, color = colors[index], size = 3.5) +
  theme(panel.grid = element_blank())
ggsave(paste0(path, 'legend2.png'), units="in", width=5, height=5, dpi=1200)
# 
# 
# # create box legend
# theme_set(theme_bw())
# df0 = data.frame(x = rep(1:2, each = 200), y = c(rep(seq(0.1,0.9,length.out =  100),2), rep(seq(0.1,2.9,length.out =  100),2)), z = factor(rep(c(rep(3,100),rep(4,100)),2) ))
# ggplot(data = df0, aes(x = factor(x), y = y, fill = z)) +
#   geom_boxplot() + xlab("Number of Machines") + ylab("Coverage") +
#   theme(panel.grid = element_blank(), legend.position = 'none')
# 
# ggsave(paste0(path, 'legend6.png'), units="in", width=5, height=5, dpi=1200)

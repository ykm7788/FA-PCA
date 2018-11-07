 
library(ggbiplot)
library(psych)


####### 输入数据 #######
otu = read.table("/home/umaru/Desktop/PCA-FA/input.txt",sep="\t",head=T,row.names=1)
data = t(otu)
group = read.table("/home/umaru/Desktop/PCA-FA/group.txt",sep="\t",row.names= 1)
sample = rownames(group)
class = as.factor(group$V2)


####### 判断变量相关性并选择成分/因子数量 #######
kmo = KMO(data)
MSA = kmo$MSA

pdf("test.parallel.pdf",width = 5,height = 5)
fa.parallel(data)
dev.off()

parallel = fa.parallel(data)
nfact = parallel$nfact
ncomp = parallel$ncomp
judgement = cbind(MSA,nfact,ncomp)
write.table (judgement, file ="test.judgement.xls", sep ="\t", row.names =FALSE,col.names =TRUE, quote =FALSE)

######### PCA #########
data.pca <- prcomp(data, scale. = TRUE)  ##R mode看因素之间关系只能用princomp##
summary.pca = summary(data.pca)
importance.pca = summary.pca$importance
rotation.pca = summary.pca$rotation
xy.pca = summary.pca$x
result.pca = rbind(importance.pca,rotation.pca,xy.pca)
write.table (result.pca, file ="test.pcadata.xls", sep ="\t", row.names =TRUE, col.names =TRUE, quote =FALSE)
 
theme_zg <- function(..., bg='white'){
  require(grid)
  theme_classic(...) +
    theme(rect=element_rect(fill=bg),
          plot.margin=unit(rep(0.5,4), 'lines'),
          panel.background=element_rect(fill='transparent', color='black'),
          panel.border=element_rect(fill='transparent', color='transparent'),
          panel.grid=element_blank(),
          axis.title.x =element_text(size=14), 
          axis.title.y=element_text(size=14),
          axis.title = element_text(color='black',vjust=0.1),
          legend.title=element_blank(),
          legend.key=element_rect(fill='transparent', color='transparent'))}

ggbiplot(data.pca, choices = 1:2,obs.scale = 1, var.scale = 1,groups = class, 
         ellipse = F, circle = F,var.axes = F,ellipse.prob = 0.68,varname.size = 1) +
  scale_color_discrete(l = 80, c = 200,h = c(0, 360)) +
  theme_zg() 
ggsave("test.PCA.pdf",width = 5,height = 5)
  
######### FA #########
data.fa<-fa(data,nfactors = 2,rotate = 'varimax',fm='ml')
loadings.fa = data.fa$loadings[1:nrow(data.fa$loadings),1:2]
weights.fa = data.fa$weights
xy.fa = data.fa$scores
result.fa = rbind(loadings.fa,weights.fa,xy.fa)
write.table (result.fa, file ="test.fadata.xls", sep ="\t", row.names =TRUE, col.names =TRUE, quote =FALSE)

pdf("test.FAvarimax.pdf",width = 5,height = 5)
fa.diagram(data.fa,simple = FALSE,main = 'fa.varimax',r = T,sort = T)
dev.off()


ggplot(as.data.frame(xy.fa),aes(x=ML1,y=ML2,colour = class)) +
  geom_point() +
  scale_color_discrete(l = 80, c = 200,h = c(0, 360)) +
  theme_zg()
ggsave("test.FA..pdf",width = 5,height = 5)


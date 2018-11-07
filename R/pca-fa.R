##安装包##
#install.packages("psych")
#install.packages("devtools", repo="http://cran.us.r-project.org")
#library(devtools)
#install_github("vqv/ggbiplot")


##加载包##
library("ggbiplot")
library("psych")


####### 输入数据 #######
otu = read.table("input.txt",sep="\t",head=T,row.names=1)
data = t(otu)
group = read.table("group.txt",sep="\t",row.names= 1)
sample = rownames(group)
class = as.factor(group$V2)


####### 判断变量相关性并选择成分/因子数量 #######
KMO(data)
fa.parallel(data)


######### PCA #########
data.pca <- prcomp(data, scale. = TRUE)  ##R mode看因素之间关系只能用princomp##
summary.pca = summary(data.pca)
importance.pca = summary.pca$importance
rotation.pca = summary$rotation
xy.pca = summary$x
result.pca = rbind(importance.pca,summary.pca$rotation,xy.pca)
write.table (result.pca, file ="pcadata.xls", sep ="\t", row.names =TRUE, col.names =TRUE, quote =FALSE)
 
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
         ellipse = F, circle = F,var.axes = T,ellipse.prob = 0.68,varname.size = 1) +
  scale_color_discrete(l = 80, c = 200,h = c(0, 360)) +
  theme_zg() 
ggsave("PCA.pdf",width = 5,height = 5)
  
######### FA #########
data.fa<-fa(data,nfactors = 2,rotate = 'varimax',fm='ml')#正交旋转的因子分析，使因子之间相互独立
loadings.fa = data.fa$loadings[1:nrow(data.fa$loadings),1:2]
weights.fa = data.fa$weights
xy.fa = data.fa$scores
result.fa = rbind(loadings.fa,weights.fa,xy.fa)
write.table (result.fa, file ="fadata.xls", sep ="\t", row.names =TRUE, col.names =TRUE, quote =FALSE)

pdf("FA.varimax.pdf",width = 5,height = 5)
fa.diagram(data.fa,simple = FALSE,main = 'fa.varimax',r = T,sort = T)
dev.off()


ggplot(as.data.frame(xy.fa),aes(x=ML1,y=ML2,colour = class)) +
  geom_point() +
  scale_color_discrete(l = 80, c = 200,h = c(0, 360)) +
  theme_zg()
ggsave("FA.pdf",width = 5,height = 5)

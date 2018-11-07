use List::Util qw/sum/;
use warnings;
use Getopt::Long;
my $VERSION = "1.0"; 
my $DATE = "2018-11-04";
my $author = 'ykm7788@qq.com';
my %opts;
GetOptions( \%opts,"i=s","g=s","o=s","pc=s","ellipse=s","circle=s","varaxes=s","varnamesize=f","rotate=s","fm=s","w=f","h=f");
my $usage = <<"USAGE";
        Program : $0
        Discription : PCA + FA Analysis 
        input data file (tab separated):
                SampleID  SampleA-1  SampleA-2  SampleA-3  SampleB-1  SampleB-2
                FactorA     xxx         xxx        xxx        xxx        xxx
                FactorB     xxx         xxx        xxx        xxx        xxx
                ...

        input group file (tab separated):
                #SampleID	Group (First row must be this)
                SampleA-1   A
                SampleA-2   A
                SampleA-3   A
                SampleB-1   B

        Version : $DATE $VERSION
                Author  : $author
        Usage : perl $0 [options]
	        -i      input data file for analysis
	        -g      input group file
	        -o      output file prefix
	        -w      the width of the graph default : 5
	        -h      the height of the graph default : 5

        PCA parameters :
                -pc             which pc to print 1:2/1:3/2:3/... default : 1:2
                -ellipse        T/F default : F
                -circle         T/F default : F
                -varaxes        T/F default : F
                -varnamesize    size of varaxes default : 1

        FA parameters :
                -rotate         rotate method: varimax quartimax bentlerT equamax varimin geominT bifactor promax oblimi simplimax bentlerQ geominQ biquartimin cluster
                -fm             fm method: minres uls ols wls gls ml minchi minrank alpha

        **  parameters without default are absolutely required  **
USAGE

die $usage if(!($opts{i}&&$opts{g}&&$opts{o}&&$opts{rotate}&&$opts{fm}));
$opts{w}=$opts{w}?$opts{w}:5;
$opts{h}=$opts{h}?$opts{h}:5;
$opts{pc}=$opts{pc}?$opts{pc}:'1:2';
$opts{ellipse}=$opts{ellipse}?$opts{ellipse}:'F';
$opts{circle}=$opts{circle}?$opts{circle}:'F';
$opts{varnamesize}=$opts{varnamesize}?$opts{varnamesize}:1;
$opts{varaxes}=$opts{varaxes}?$opts{varaxes}:'F';

my $data_file = $opts{i};
my $group_file = $opts{g};
my $prefix = $opts{o};
my $width = $opts{w};
my $height = $opts{h};


open CMD,">cmd.r";
select CMD ;
print CMD " 
library(ggbiplot)
library(psych)


####### 输入数据 #######
otu = read.table(\"$data_file\",sep=\"\\t\",head=T,row.names=1)
data = t(otu)
group = read.table(\"$group_file\",sep=\"\\t\",row.names= 1)
sample = rownames(group)
class = as.factor(group\$V2)


####### 判断变量相关性并选择成分/因子数量 #######
kmo = KMO(data)
MSA = kmo\$MSA

pdf(\"$prefix\.parallel.pdf\",width = $width,height = $height)
fa.parallel(data)
dev.off()

parallel = fa.parallel(data)
nfact = parallel\$nfact
ncomp = parallel\$ncomp
judgement = cbind(MSA,nfact,ncomp)
write.table (judgement, file =\"$prefix\.judgement.xls\", sep =\"\\t\", row.names =FALSE,col.names =TRUE, quote =FALSE)

######### PCA #########
data.pca <- prcomp(data, scale. = TRUE)  ##R mode看因素之间关系只能用princomp##
summary.pca = summary(data.pca)
importance.pca = summary.pca\$importance
rotation.pca = summary.pca\$rotation
xy.pca = summary.pca\$x
result.pca = rbind(importance.pca,rotation.pca,xy.pca)
write.table (result.pca, file =\"$prefix\.pcadata\.xls\", sep =\"\\t\", row.names =TRUE, col.names =TRUE, quote =FALSE)
 
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

ggbiplot(data.pca, choices = $opts{pc},obs.scale = 1, var.scale = 1,groups = class, 
         ellipse = $opts{ellipse}, circle = $opts{circle},var.axes = $opts{varaxes},ellipse.prob = 0.68,varname.size = $opts{varnamesize}) +
  scale_color_discrete(l = 80, c = 200,h = c(0, 360)) +
  theme_zg() 
ggsave(\"$prefix\.PCA\.pdf\",width = $width,height = $height)
  
######### FA #########
data.fa<-fa(data,nfactors = 2,rotate = \'$opts{rotate}\',fm=\'$opts{fm}\')
loadings\.fa = data.fa\$loadings[1:nrow(data\.fa\$loadings),1:2]
weights.fa = data.fa\$weights
xy.fa = data.fa\$scores
result.fa = rbind(loadings.fa,weights.fa,xy.fa)
write.table (result.fa, file =\"$prefix\.fadata\.xls\", sep =\"\\t\", row.names =TRUE, col.names =TRUE, quote =FALSE)

pdf(\"$prefix\.FAvarimax.pdf\",width = $width,height = $height)
fa.diagram(data.fa,simple = FALSE,main = \'fa.varimax\',r = T,sort = T)
dev.off()


ggplot(as.data.frame(xy.fa),aes(x=ML1,y=ML2,colour = class)) +
  geom_point() +
  scale_color_discrete(l = 80, c = 200,h = c(0, 360)) +
  theme_zg()
ggsave(\"$prefix\.FA.\.pdf\",width = $width,height = $height)

";
`R --restore --no-save < cmd.r`;

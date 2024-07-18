options("repos" = c(CRAN="https://mirrors.ustc.edu.cn/CRAN/")) #中科大镜像
options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/") #中科大镜像
# 网络优化
options(download.file.method = 'libcurl')
options(url.method='libcurl')
library(tidyverse)
library(data.table)

# count→tpm
func_count2tpm = function(count_data, annotation_gtf_dir, count2what = "TPM"){
  #annotation_gtf_dir = "./BZP02/00_raw_data/Xena/gencode.v22.annotation.gtf.gz"
  #count_data = exp_KIRC.count
  # 安装必备包
  pkg = c("txdbmaker", "GenomicFeatures", "DGEobj.utils", "AnnoProbe")
  sapply(pkg, simplify = T, function(x){
    if (! require(x,character.only=T) ) {
      BiocManager::install(x,ask = F,update = F) 
      require(x,character.only=T) 
    }})
  library(GenomicFeatures)
  library(DGEobj.utils)
  library(AnnoProbe)
  # 读取基因组注释文件
  txdb = makeTxDbFromGFF(annotation_gtf_dir, format = "gtf")
  # 获取每个基因id的外显子数据
  exons.list.per.gene = exonsBy(txdb, by="gene")
  # 对于每个基因，将所有外显子减少成一组非重叠外显子，计算它们的长度(宽度)并求和
  exonic.gene.sizes = sum(width(GenomicRanges::reduce(exons.list.per.gene)))
  # 得到geneid和长度数据
  gfe = data.frame(gene_id=names(exonic.gene.sizes),
                   length=exonic.gene.sizes)
  # 去除版本号
  gfe$gene_id = str_split(gfe$gene_id, "\\.", simplify = T)[,1]
  
  # 同步矩阵和gfe内容和顺序
  filter = intersect(rownames(count_data), gfe$gene_id)
  count_data = count_data[rownames(count_data) %in% filter, ]
  gfe = gfe[gfe$gene_id %in% filter, ]
  if(!identical(rownames(count_data), gfe$gene_id)){
    gfe = gfe[match(rownames(count_data), gfe$gene_id), ]
  }
  if(!identical(rownames(count_data), gfe$gene_id)){
    stop("错误！基因长度列表并未和矩阵行名完成同步！")
  }
  # 进行转换
  res = convertCounts(
    countsMatrix = as.matrix(count_data),
    unit = count2what,
    geneLength = gfe$length,
    log = FALSE,
    normalize = "none",
    prior.count = NULL
  )
  return(res)
}

# 获取连续及离散调色盘
func_getPalettes = function(continuousPalette_floor = "blue3",
                            continuousPalette_middle = "grey",
                            continuousPalette_celling = "red2",
                            continuous_floor = -10,
                            continuous_middle = 0,
                            continuous_celling = 10){
  pkg = c("RColorBrewer", "circlize")
  sapply(pkg, simplify = T, function(x){
    if(!require(x,character.only=T)){
      install.packages(x,ask = F,update = F)
      require(x,character.only=T) 
    }})
  library(RColorBrewer)
  library(circlize)
  
  ## 构建颜色向量表
  ### 获取调色盘
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  ###处理后有74种差异还比较明显的颜色
  discretePalette = unlist(mapply(brewer.pal,qual_col_pals$maxcolors, rownames(qual_col_pals))) 
  ### 直接从调色盘中选色
  #pie(rep(1,length(discretePalette)), col = util_continuousPalette, radius = 1.07)
  imgUrl_discretePalette = "https://haokun-img-storage-1302331098.cos.ap-chengdu.myqcloud.com/img/discrete_palette.png"
  getImg_discretePalette = "pie(rep(1,length(discretePalette)), col = util_continuousPalette, radius = 1.07)"
  ##获取连续调色盘
  continuousPalette = colorRamp2(
    c(continuous_floor, continuous_middle, continuous_celling),  #根据值的范围设置
    c(continuousPalette_floor, continuousPalette_middle, continuousPalette_celling)
  )
  return(list(
    discretePalette = discretePalette,
    continuousPalette = continuousPalette,
    imgUrl_discretePalette = imgUrl_discretePalette,
    getImg_discretePalette = getImg_discretePalette
  ))
}



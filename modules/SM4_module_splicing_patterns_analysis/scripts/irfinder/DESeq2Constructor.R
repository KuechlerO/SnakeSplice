DESeqDataSetFromIRFinder = function(filePaths,designMatrix,designFormula){
    res=c()
    libsz=c()
    spl=c()
    irtest=read.table(filePaths[1])
    if (irtest[1,1]=="Chr"){irtest=irtest[-1,]}
    irnames=unname(apply(as.matrix(irtest),1,FUN=function(x){return(paste0(x[4],"/",x[1],":",x[2],"-",x[3],":",x[6]))}))
    n=1
    for (i in filePaths){
        print(paste0("processing file ",n," at ",i))
        irtab=read.table(i)
        if (irtab[1,1]=="Chr"){irtab=irtab[-1,]}
        #rn=unname(apply(irtab,1,FUN=function(x){return(paste0(x[4],"/",x[1],":",x[2],"-",x[3],":",x[6]))}))
        #row.names(irtab)=rn
        #tmp1=round(as.numeric(as.vector(irtab[irnames,9])))
        #tmp2=as.numeric(as.vector(irtab[irnames,19]))
        tmp1=as.numeric(as.vector(irtab[,9]))
        tmp2=as.numeric(as.vector(irtab[,19]))
        tmp3=tmp1+tmp2
        tmp4=as.numeric(as.vector(irtab[,17]))
        tmp5=as.numeric(as.vector(irtab[,18]))
        tmp6=pmax(tmp4,tmp5, na.rm=T)
        res=cbind(res,tmp1)
        libsz=cbind(libsz,tmp2)
        spl=cbind(spl,tmp6)
        n=n+1
    }
    res.rd=round(res)
    libsz.rd=round(libsz)
    spl.rd=round(spl)
    colnames(res.rd)=paste("intronDepth",as.vector(designMatrix[,1]),sep=".")
    rownames(res.rd)=irnames
    colnames(libsz.rd)=paste("totalSplice",as.vector(designMatrix[,1]),sep=".")
    rownames(libsz.rd)=irnames
    colnames(spl.rd)=paste("maxSplice",as.vector(designMatrix[,1]),sep=".")
    rownames(spl.rd)=irnames

    ir=c(rep("IR",dim(designMatrix)[1]),rep("Splice",dim(designMatrix)[1]))
    group=rbind(designMatrix,designMatrix)
    group$IRFinder=ir
    group$IRFinder=factor(group$IRFinder,levels=c("Splice","IR"))

    #counts.IRFinder=cbind(res.rd,libsz.rd)
    counts.IRFinder=cbind(res.rd,spl.rd)

    dd = DESeqDataSetFromMatrix(countData = counts.IRFinder, colData = group, design = designFormula)
    sizeFactors(dd)=rep(1,dim(group)[1])
    rownames(dd)=irnames
    final=list(dd,res,libsz,spl)
    names(final)=c("DESeq2Object","IntronDepth","SpliceDepth","MaxSplice")
    return(final)
}

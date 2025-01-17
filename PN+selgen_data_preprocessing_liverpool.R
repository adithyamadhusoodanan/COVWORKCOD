  library(plyr)
  library(dplyr)
  library(reshape2)
  library(ape)
  library(seqinr)
  library(matrixStats)
  `%not in%` <- function (x, table) is.na(match(x, table, nomatch=NA_integer_))
  
  # Start h2o JVM
  
  # Read data from file
  f1<-read.csv(file="liverpool_features_genus.csv",header=T)
  f1 <- rename(f1, Genbank.accession=SeqName)
  f1 <- subset(f1, select = -c(X,Unnamed..0,genus)) # remove X 
  allP<-read.fasta(file ="WGS_liverpool_data.fas",  seqtype = "DNA", as.string = TRUE, seqonly = F, strip.desc = T)
  fis<-read.csv(file="featureImportance_reservoir.csv",header=T)
  
  fsars2<-read.csv(file="sars_liverpool_features.csv",header=T)
  fsars2 <- rename(fsars2, Genbank.accession=SeqName)
  fsars2 <- subset(fsars2, select = -c(X)) # remove X
  sars2fas<-read.fasta(file ="sars_wholegenome.fasta",  seqtype = "DNA", as.string = TRUE, seqonly = F, strip.desc = T)
  
  # Feature definition
  dinucs<-grep("[A|T|G|C|U]p[A|T|G|C|U]",names(f1),value=T)
  cps<-grep(".[A|C|D|E|F|G|H|I|K|L|M|N|P|Q|R|S|T|V|W|X|Y]..[A|T|G|C|U]",names(f1),value=T)
  aa.codon.bias<-grep(".Bias",names(f1),value=T)
  
  # Feature selection (simplify dataset to required columns)
  
  # Run the script for 25 50 100 and 200  features.
  nfeats<-25  # then tr y 100, 75, 25
  
  totalfeats<-length(fis$vimean)
  f<-seq(from = totalfeats-(nfeats-1),to = totalfeats, by=1)
  gen.feats<-as.character(fis$X[f])
  
  f1<-f1[,c("Genbank.accession","Reservoir","virus", gen.feats)]
  
  # Remove orphans
  f2<-subset(f1,f1$Reservoir!="Orphan")
  f<-droplevels(f2)
  
    # Group selection based on thresholds
  t<-15 # threshold for minimum sample size of groups
  s<-.85 # proportion in the training set
  
  host.counts<-table(f$Reservoir)
  
  min.t<-host.counts[host.counts>=t] # minimum number of viruses per host group
  f_st3<-f[f$Reservoir %in% c(names(min.t)),]
  f_st3<-droplevels(f_st3)
  f_st3$SeqName2<-do.call(rbind,strsplit(as.character(f_st3$Genbank.accession),"[.]"))[,1]
  
  # Rare hosts
  rare<-f[!f$Reservoir %in% c(names(min.t)),]
  rare<-droplevels(rare)
  rare$SeqName2<-do.call(rbind,strsplit(as.character(rare$Genbank.accession),"[.]"))[,1]
  
  # Number and names of host taxa
  ntax<-length(unique(f_st3$Reservoir))
  bp<-as.character(sort(unique(f_st3$Reservoir)))
  
  # sample split of training/test to get counts in each
  
  set.seed(78910)
  
  #trains<-f_st3 %>% group_by(Reservoir) %>%
  # filter(Genbank.accession %in% sample(unique(Genbank.accession), ceiling(s*length(unique(Genbank.accession)))))
  #testval<-subset(f_st3,!(f_st3$Genbank.accession %in% trains$Genbank.accession)) # ref numbers absent from training set
  
  #optims<-testval %>% group_by(Reservoir) %>%
  #  filter(Genbank.accession %in% sample(unique(Genbank.accession), floor(.5*length(unique(Genbank.accession)))))
  #tests<-subset(testval,!(testval$Genbank.accession %in% optims$Genbank.accession)) # ref numbers in testval set absent from test set
  #ntest<-dim(tests)[1]
  
  # write orphan sequences
  # orp<-allP[c(which(names(allP) %in% orphans$Genbank.accession))]
  # write.fasta(orp,names(orp),file.out="orphanDB.fasta", open = "w", nbchar = 100, as.string = T)
  
  # write rare host sequences
  #rar<-allP[c(which(names(allP) %in% rare$Genbank.accession))]
  #write.fasta(rar,names(rar),file.out="rareDB.fasta", open = "w", nbchar = 100, as.string = T)
  
  # Remove unneeded files 
  #     rm(f,f1,f2,fis,rar,orp)
  # Train many models
  set.seed(78910)
  
  uni_val<-unique(f_st3$virus)
  
  nloops<-length(uni_val)
  #lr<-c()
  #md<-c()
  #sr<-c()
  #csr<-c()
  #nt<-c()
  # mr<-c()
  # accuracy.st3<-c()
  
  # pc.accuracy<-matrix(nrow=nloops,ncol=ntax)
  # test.record<-matrix(nrow=ntest,ncol=nloops)
  nfeatures<-length(gen.feats)+ntax
  vimps<-matrix(nrow=nfeatures,ncol=nloops)
  for (i in 1:nloops) {
      # Stratified random sampling
      trains<-f_st3[f_st3$virus!=uni_val[i],]
      tests<-f_st3[f_st3$virus==uni_val[i],]

      trains<-droplevels(trains)
      fsars2<-droplevels(fsars2)
      tests<-droplevels(tests)
      #optims<-droplevels(optims)
      #   test.record[,i]<-as.character(tests$Genbank.accession)
  
      # Select and write sequences to local directory
      trainSeqs<-allP[c(which(names(allP) %in% trains$Genbank.accession))] # pick sequences in the training set
      testSeqs<-allP[c(which(names(allP) %in% tests$Genbank.accession))] # pick sequences in the validation set
      sars2Seqs<-sars2fas[c(which(names(sars2fas) %in% fsars2$Genbank.accession))]
      
      #optSeqs<-allP[c(which(names(allP) %in% optims$Genbank.accession))] # pick sequences in the optimization set
      write.fasta(testSeqs, names(testSeqs), file.out="testDB.fasta", open = "w", nbchar = 100, as.string = T)
      write.fasta(trainSeqs, names(trainSeqs), file.out="trainDB.fasta", open = "w", nbchar = 100, as.string = T)
      #write.fasta(optSeqs, names(optSeqs), file.out="optDB.fasta", open = "w", nbchar = 100, as.string = T)
      write.fasta(sars2Seqs, names(sars2Seqs), file.out="sars2DB.fasta", open = "w", nbchar = 100, as.string = T)
  
      # BLAST
      system("makeblastdb -in trainDB.fasta -dbtype nucl -parse_seqids -out allTrainingDB",intern=F)
  
      # Blast test against training
      system("blastn -db allTrainingDB -query testDB.fasta -out testOut.out -num_threads 4 -outfmt 10 -max_target_seqs=5 -max_hsps 1 -reward 2 -task blastn -evalue 10 -word_size 8 -gapopen 2 -gapextend 2",inter=F,wait=FALSE)
  
      # Blast validation against training
     # system("blastn -db allTrainingDB -query optDB.fasta -out optOut.out -num_threads 4 -outfmt 10 -max_target_seqs=5 -max_hsps 1 -reward 2 -task blastn -evalue 10 -word_size 8 -gapopen 2 -gapextend 2",inter=F,wait=FALSE)
  
      system("blastn -db allTrainingDB -query sars2DB.fasta -out sars2Out.out -num_threads 4 -outfmt 10 -max_target_seqs=5 -max_hsps 1 -reward 2 -task blastn -evalue 10 -word_size 8 -gapopen 2 -gapextend 2",inter=F,wait=FALSE)
      
      # Blast training against the training set (take top 5 hits)
      system("blastn -db allTrainingDB -query trainDB.fasta -out trainOut.out -num_threads 4 -outfmt 10 -max_target_seqs=5 -max_hsps 1 -reward 2 -task blastn -evalue 10 -word_size 8 -gapopen 2 -gapextend 2",inter=F,wait=TRUE)
  
      # Read in Blast output for training set
      allBlast<-read.csv(file="trainOut.out",col.names = c("query acc.", "subject acc.", "% identity", "alignment length", "mismatches", "gap opens", "q. start", "q. end"," s. start"," s. end"," evalue"," bit score"),header=F)
  
      # Summarize blast hits
      nvir<-length(unique(allBlast$query.acc.))
      virnames<-unique(allBlast$query.acc.)
      ecutoff<-1E-3
      j=1
      d<-subset(allBlast,allBlast$query.acc.==virnames[j])
      d2<-subset(d,d$X..identity<100)
      d2<-subset(d2,d2$X.evalue<ecutoff)
      # Assign equal probability across all hosts if there is no good blast hit
      for (z in 1:1){
          if (nrow(d2)==0){
              blast.uc<-rep(1/ntax,ntax)
              blast.uc<-data.frame(t(blast.uc))
              colnames(blast.uc)<-sort(unique(trains$Reservoir))
              id<-as.character(virnames[j])
              blast.uc<-cbind(id,blast.uc)}
          else {
              dhost<-merge(d2,trains,by.x="subject.acc.",by.y="Genbank.accession",all.x = T,all.y = F)
              dhost$rel.support<-dhost$X..identity/sum(dhost$X..identity)
              hosts<-tapply(dhost$rel.support,dhost$Reservoir,sum,na.rm=F)
              hosts[is.na(hosts)]<-0
              hosts<-t(data.frame(hosts))
              hosts<-data.frame(hosts)
              id<-as.character(virnames[j])
              blast.uc<-cbind(id,hosts)}}
  
      for (j in 2:nvir){
          d<-subset(allBlast,allBlast$query.acc.==virnames[j])
          d2<-subset(d,d$X..identity<100)
          d2<-subset(d2,d2$X.evalue<ecutoff)
          if (nrow(d2)==0){
              blast.uc.s<-rep(1/ntax,ntax)
              blast.uc.s<-data.frame(t(blast.uc.s))
              colnames(blast.uc.s)<-sort(unique(trains$Reservoir))
              id<-as.character(virnames[j])
              blast.uc.s<-cbind(id,blast.uc.s)}
          else {
              dhost<-merge(d2,trains,by.x="subject.acc.",by.y="Genbank.accession",all.x = T,all.y = F)
              dhost$rel.support<-dhost$X..identity/sum(dhost$X..identity)
              hosts<-tapply(dhost$rel.support,dhost$Reservoir,sum,na.rm=T)
              hosts[is.na(hosts)]<-0
              hosts<-t(data.frame(hosts))
              hosts<-data.frame(hosts)
              id<-as.character(virnames[j])
              blast.uc.s<-cbind(id,hosts)
              }
          blast.uc<-rbind.fill(blast.uc,blast.uc.s)
          blast.uc[is.na(blast.uc)] <- 0
          }
  
      f1_train<-merge(trains,blast.uc,by.x="Genbank.accession",by.y="id",all.x=F,all.y=T)
      set<-c("Reservoir","Genbank.accession",gen.feats,bp)
      myData=data.frame(matrix(nrow=1, ncol = length(set)))
      f1_train=bind_rows(myData,f1_train)
      f1_train<-f1_train[,c(set)] # this is the full training dataset with genomic features and blast probabilities
      f1_train=f1_train[-1,]
      f1_train[is.na(f1_train)]=0
  
      # Summarize blast hits from test set
      testBlast<-read.csv(file="testOut.out",col.names = c("query acc.", "subject acc.", "% identity", "alignment length", "mismatches", "gap opens", "q. start", "q. end"," s. start"," s. end"," evalue"," bit score"),header=F)
      nvir<-length(unique(testBlast$query.acc.))
      virnames<-unique(testBlast$query.acc.)
      ecutoff<-1E-3
      j=1
      d<-subset(testBlast,testBlast$query.acc.==virnames[j])
      d2<-subset(d,d$X.evalue<ecutoff)
  
      for (z in 1:1){
          if (nrow(d2)==0){
              blast.uc<-rep(1/ntax,ntax)
              blast.uc<-data.frame(t(blast.uc))
              colnames(blast.uc)<-sort(unique(trains$Reservoir))
              id<-as.character(virnames[j])
              blast.uc<-cbind(id,blast.uc)}
          else {
              dhost<-merge(d2,trains,by.x="subject.acc.",by.y="Genbank.accession",all.x = T,all.y = F)
              dhost$rel.support<-dhost$X..identity/sum(dhost$X..identity)
              hosts<-tapply(dhost$rel.support,dhost$Reservoir,sum,na.rm=F)
              hosts[is.na(hosts)]<-0
              hosts<-t(data.frame(hosts))
              hosts<-data.frame(hosts)
              id<-as.character(virnames[j])
              blast.uc<-cbind(id,hosts)}}
  
      for (j in 2:nvir){
          d<-subset(testBlast,testBlast$query.acc.==virnames[j])
          d2<-subset(d,d$X.evalue<ecutoff)
          if (nrow(d2)==0){
              blast.uc.s<-rep(1/ntax,ntax)
              blast.uc.s<-data.frame(t(blast.uc.s))
              colnames(blast.uc.s)<-sort(unique(trains$Reservoir))
              id<-as.character(virnames[j])
              blast.uc.s<-cbind(id,blast.uc.s) }
          else {
              dhost<-merge(d,trains,by.x="subject.acc.",by.y="Genbank.accession",all.x = T,all.y = F)
              dhost$rel.support<-dhost$X..identity/sum(dhost$X..identity)
              hosts<-tapply(dhost$rel.support,dhost$Reservoir,sum,na.rm=T)
              hosts[is.na(hosts)]<-0
              hosts<-t(data.frame(hosts))
              hosts<-data.frame(hosts)
              id<-as.character(d$query.acc.[1])
              blast.uc.s<-cbind(id,hosts)}
          blast.uc<-rbind.fill(blast.uc,blast.uc.s)
          blast.uc[is.na(blast.uc)] <- 0
          }
  
      f1_test<-merge(tests,blast.uc,by.x="Genbank.accession",by.y="id",all.x=F,all.y=T)
      testID<-f1_test$Virus.name
      # adding in the values for reservoir groups with no relative support
      set<-c("Reservoir","Genbank.accession",gen.feats,bp)
      myData=data.frame(matrix(nrow=1, ncol = length(set)))
      colnames(myData) = set
      f1_test=bind_rows(myData,f1_test)
      f1_test<-f1_test[,c(set)]
      f1_test=f1_test[-1,]
      f1_test[is.na(f1_test)]=0
      
      #cleaning the output for n=1 in testing set
      if((nrow(tests)==1)&(nrow(f1_test)==3)) {
      f1_test=f1_test[2,]
      }
      
      ' Summarize blast hits from optimization set
      optBlast<-read.csv(file="optOut.out",col.names = c("query acc.", "subject acc.", "% identity", "alignment length", "mismatches", "gap opens", "q. start", "q. end"," s. start"," s. end"," evalue"," bit score"),header=F)
      nvir<-length(unique(optBlast$query.acc.))
      virnames<-unique(optBlast$query.acc.)
      ecutoff<-1E-3
      j=1
      d<-subset(optBlast,optBlast$query.acc.==virnames[j])
      d2<-subset(d,d$X.evalue<ecutoff)
  
      for (z in 1:1){
        if (nrow(d2)==0){
          blast.uc<-rep(1/ntax,ntax)
          blast.uc<-data.frame(t(blast.uc))
          colnames(blast.uc)<-sort(unique(trains$Reservoir))
          id<-as.character(virnames[j])
          blast.uc<-cbind(id,blast.uc)}
        else {
          dhost<-merge(d2,trains,by.x="subject.acc.",by.y="Genbank.accession",all.x = T,all.y = F)
          dhost$rel.support<-dhost$X..identity/sum(dhost$X..identity)
          hosts<-tapply(dhost$rel.support,dhost$Reservoir,sum,na.rm=F)
          hosts[is.na(hosts)]<-0
          hosts<-t(data.frame(hosts))
          hosts<-data.frame(hosts)
          id<-as.character(virnames[j])
          blast.uc<-cbind(id,hosts)}}
  
      for (j in 2:nvir){
        d<-subset(optBlast,optBlast$query.acc.==virnames[j])
        d2<-subset(d,d$X.evalue<ecutoff)
        if (nrow(d2)==0){
          blast.uc.s<-rep(1/ntax,ntax)
          blast.uc.s<-data.frame(t(blast.uc.s))
          colnames(blast.uc.s)<-sort(unique(trains$Reservoir))
          id<-as.character(virnames[j])
          blast.uc.s<-cbind(id,blast.uc.s) }
        else {
          dhost<-merge(d,trains,by.x="subject.acc.",by.y="Genbank.accession",all.x = T,all.y = F)
          dhost$rel.support<-dhost$X..identity/sum(dhost$X..identity)
          hosts<-tapply(dhost$rel.support,dhost$Reservoir,sum,na.rm=T)
          hosts[is.na(hosts)]<-0
          hosts<-t(data.frame(hosts))
          hosts<-data.frame(hosts)
          id<-as.character(d$query.acc.[1])
          blast.uc.s<-cbind(id,hosts)}
        blast.uc<-rbind.fill(blast.uc,blast.uc.s)
        blast.uc[is.na(blast.uc)] <- 0
        }
  
      f1_opt<-merge(optims,blast.uc,by.x="Genbank.accession",by.y="id",all.x=F,all.y=T)
      optID<-f1_opt$Virus.name
      set<-c("Reservoir",gen.feats,bp)
      f1_opt<-f1_opt[,c(set)]'
  
      # save train, test, opt sets
      
      sars2Blast<-read.csv(file="sars2Out.out",col.names = c("query acc.", "subject acc.", "% identity", "alignment length", "mismatches", "gap opens", "q. start", "q. end"," s. start"," s. end"," evalue"," bit score"),header=F)
      nvir<-length(unique(sars2Blast$query.acc.))
      virnames<-unique(sars2Blast$query.acc.)
      ecutoff<-1E-6
      j=1
      d<-subset(sars2Blast,sars2Blast$query.acc.==virnames[j])
      d2<-subset(d,d$X.evalue<ecutoff)
      
      for (z in 1:1){
        if (nrow(d2)==0){
          blast.uc<-rep(1/ntax,ntax)
          blast.uc<-data.frame(t(blast.uc))
          colnames(blast.uc)<-sort(unique(trains$Reservoir))
          id<-as.character(virnames[j])
          blast.uc<-cbind(id,blast.uc)}
        else {
          dhost<-merge(d2,trains,by.x="subject.acc.",by.y="Genbank.accession",all.x = T,all.y = F)
          dhost$rel.support<-dhost$X..identity/sum(dhost$X..identity)
          hosts<-tapply(dhost$rel.support,dhost$Reservoir,sum,na.rm=F)
          hosts[is.na(hosts)]<-0
          hosts<-t(data.frame(hosts))
          hosts<-data.frame(hosts)
          id<-as.character(virnames[j])
          blast.uc<-cbind(id,hosts)}}
      
      for (j in 2:nvir){
        d<-subset(sars2Blast,sars2Blast$query.acc.==virnames[j])
        d2<-subset(d,d$X.evalue<ecutoff)
        if (nrow(d2)==0){
          blast.uc.s<-rep(1/ntax,ntax)
          blast.uc.s<-data.frame(t(blast.uc.s))
          colnames(blast.uc.s)<-sort(unique(trains$Reservoir))
          id<-as.character(virnames[j])
          blast.uc.s<-cbind(id,blast.uc.s) }
        else {
          dhost<-merge(d,trains,by.x="subject.acc.",by.y="Genbank.accession",all.x = T,all.y = F)
          dhost$rel.support<-dhost$X..identity/sum(dhost$X..identity)
          hosts<-tapply(dhost$rel.support,dhost$Reservoir,sum,na.rm=T)
          hosts[is.na(hosts)]<-0
          hosts<-t(data.frame(hosts))
          hosts<-data.frame(hosts)
          id<-as.character(d$query.acc.[1])
          blast.uc.s<-cbind(id,hosts)}
        blast.uc<-rbind.fill(blast.uc,blast.uc.s)
        blast.uc[is.na(blast.uc)] <- 0
      }
      
      f1_sars2<-merge(fsars2,blast.uc,by.x="Genbank.accession",by.y="id",all.x=F,all.y=T)
      sars2ID<-f1_sars2$Virus.name
      set<-c("Reservoir","Genbank.accession",gen.feats,bp)
      myData=data.frame(matrix(nrow=1, ncol = length(set)))
      colnames(myData) = set
      f1_sars2=bind_rows(myData,f1_sars2)
      f1_sars2<-f1_sars2[,c(set)]
      f1_sars2=f1_sars2[-1,]
      f1_sars2[is.na(f1_sars2)]=0
   
      write.csv(f1_train,file=paste("results/training_set",i,".csv",sep="_"))
      write.csv(f1_sars2,file=paste("results/sars_set",i,".csv",sep="_"))
      write.csv(f1_test,file=paste("results/test_set",i,".csv",sep="_"))
      print(i)
      system("git add .",intern=F)
      system("git commit",intern=F)
      system("git push https://adithyamadhusoodanan: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@github.com/adithyamadhusoodanan/COVWORKCOD.git",intern=F)
   
  }



###############################################################################
# TEST COMMUNITIES PARTITIONS
# Copyright (C) 2022
#
# This code is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. This code is distributed in the hope that it will be
# useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
# Public License for more details.
#
# Elaine Cecilia Gatto | Prof. Dr. Ricardo Cerri | Prof. Dr. Mauri Ferrandin
# Federal University of Sao Carlos (UFSCar: https://www2.ufscar.br/) Campus
# Sao Carlos Computer Department (DC: https://site.dc.ufscar.br/)
# Program of Post Graduation in Computer Science
# (PPG-CC: http://ppgcc.dc.ufscar.br/)
# Bioinformatics and Machine Learning Group
# (BIOMAL: http://www.biomal.ufscar.br/)
#
###############################################################################


FolderRoot = "~/TCP-KNN-H-Clus"
FolderScripts = paste(FolderRoot, "/R", sep="")


##########################################################################
# FUNCTION COMPUTE SILHOUETE
##########################################################################
silhouete <- function(parameters){

  f = 1
  silhoueteParalel <- foreach(f = 1:parameters$Number.Folds) %dopar% {
  #while(f<=number_folds){

    cat("\n#=========================================================")
    cat("\n# FOLD = \t", f)
    cat("\n#=========================================================")

    #########################################################################
    FolderRoot = "~/TCP-KNN-H-Clus"
    FolderScripts = paste(FolderRoot, "/R", sep="")

    setwd(FolderScripts)
    source("libraries.R")

    setwd(FolderScripts)
    source("utils.R")

    ##########################################################################
    constroiParticoes <- function(TotalParticoes){

      data <- data.frame(matrix(NA,    # Create empty data frame
                                nrow = TotalParticoes,
                                ncol = 2))

      names(data) = c("numberPartition", "numberGroup")

      i = 1
      a = 1
      while(i<=nrow(data)){
        data[i,1] = a + 1
        data[i,2] = a + 1
        i = i + 1
        a = a + 1
        gc()
      }

      return(data)

    }

    fold = c(0)
    knn = c(0)
    part = c(0)
    maximo = c(0)
    minimo = c(0)
    mediana = c(0)
    media = c(0)
    primeiroQuadrante = c(0)
    terceiroQuadrante = c(0)
    valueSilhouete = c(0)
    bestPartition = data.frame(fold, knn, part, maximo, minimo,
                               mediana, media, primeiroQuadrante,
                               terceiroQuadrante, valueSilhouete)

    ########################################################################
    #  "/dev/shm/j-GpositiveGO/Validation/Split-1"
    FolderSplitVal = paste(parameters$Folders$folderValSilho,
                           "/Split-", f, sep="")
    if(dir.exists(FolderSplitVal)==FALSE){dir.create(FolderSplitVal)}

    # "/dev/shm/j-GpositiveGO/Communities/Split-1"
    FolderSplitComm = paste(parameters$Folders$folderCommunities,
                            "/Split-", f, sep="")

    # "/dev/shm/j-GpositiveGO/Partitions/Split-1"
    FolderPSplit = paste(parameters$Folders$folderPartitions,
                         "/Split-",f, sep="")

    ########################################################################
    # get the space label
    espacoDeRotulos = data.frame(parameters$LabelSpace$Classes[f])
    espacoDeRotulos2 = data.frame(t(espacoDeRotulos))
    labels = rownames(espacoDeRotulos2)
    espacoDeRotulos2 = cbind(labels, espacoDeRotulos2)
    espacoDeRotulos2 = data.frame(espacoDeRotulos2[order(espacoDeRotulos2$labels,
                                                         decreasing = FALSE),])

    ########################################################################
    #cat("\n\nGrupos por particão")
    knn_H = data.frame(read.csv(paste(FolderPSplit,
                                      "/fold-", f,
                                      "-knn-h-choosed.csv", sep="")))
    total_knn_H = nrow(knn_H)

    at = data.frame(fold = f, sparc = total_knn_H)

    write.csv(at, paste(parameters$Folders$folderValSilho,
                               "/fold-", f, "-sparcification.csv",
                        sep=""),
              row.names = FALSE)

    # do primeiro knn ao ultimo
    k = 1
    while(k<=total_knn_H){

      cat("\n#=========================================================")
      cat("\n# Knn = ", k)
      cat("\n#=========================================================")

      # "/dev/shm/j-GpositiveGO/Partitions/Split-1/knn-1"
      FolderPartKnn = paste(FolderPSplit, "/knn-", k, sep="")

      # "/dev/shm/j-GpositiveGO/Communities/Split-1/knn-1"
      FolderPartComm = paste(FolderSplitComm, "/knn-", k, sep="")

      # "/dev/shm/j-GpositiveGO/Partitions/Split-1"
      FolderPSplit = paste(parameters$Folders$folderPartitions, "/Split-",f, sep="")

      # "/dev/shm/j-GpositiveGO/Validation/Split-1/knn-1"
      FolderKnnVal = paste(FolderSplitVal, "/knn-", k, sep="")
      if(dir.exists(FolderKnnVal)==FALSE){dir.create(FolderKnnVal)}

      knn_H = data.frame(read.csv(paste(FolderPSplit,
                                        "/fold-", f,
                                        "-knn-h-choosed.csv", sep="")))
      total_knn_H = nrow(knn_H)

      knn_H = knn_H[k,]

      particoes = data.frame(read.csv(paste(FolderPartComm,
                                            "/knn-", k,
                                            "-", knn_H$method,
                                            "-partitions-hierarchical.csv",
                                            sep="")))
      TotalParticoes = ncol(particoes)-1
      numPart = (parameters$Dataset.Info$Labels-1)

      ########################################################################
      fold = c(0)
      part = c(0)
      knn = c (0)
      maximo = c(0)
      minimo = c(0)
      mediana = c(0)
      media = c(0)
      primeiroQuadrante = c(0)
      terceiroQuadrante = c(0)
      valueSilhouete = c(0)
      Silhouete = data.frame(fold, part, knn, maximo, minimo, mediana, media,
                             primeiroQuadrante, terceiroQuadrante,
                             valueSilhouete)

      p = 2
      while(p<=numPart){

        cat("\n#=========================================================")
        cat("\n# Partition = ", p)
        cat("\n#=========================================================")

        # "/dev/shm/j-GpositiveGO/Validation/Split-1/knn-1/Partition-2"
        FolderPartVal = paste(FolderKnnVal, "/Partition-", p, sep="")
        if(dir.exists(FolderPartVal)==FALSE){dir.create(FolderPartVal)}

        ########################################################################
        cat("\n get the number of groups for this partition")
        particao = particoes[,c(1,p)]
        names(particao) = c("labels", "groups")

        res = constroiParticoes(TotalParticoes)
        res = filter(res, numberPartition == p)
        numGroups = as.numeric(res$numberGroup)

        #######################################################################
        if(numGroups==1){
          cat("\nOnly one group of labels (global partition)")
          fold = f
          part = p
          knn = k
          maximo = NA
          minimo = NA
          mediana = NA
          media = NA
          primeiroQuadrante = NA
          terceiroQuadrante = NA
          valueSilhouete = NA
          Silhouete = rbind(Silhouete, data.frame(fold, part, knn, maximo,
                                                  minimo, mediana, media,
                                                  primeiroQuadrante,
                                                  terceiroQuadrante,
                                                  valueSilhouete))
          write.csv(Silhouete[-1,], paste(FolderKnnVal, "/fold-", f, "-knn-",
                                          k, "-silho.csv", sep=""),
                    row.names = FALSE)

        } else {
          #cat("\ntwo or more labels in the group")
          groups_label_space = cbind(particao, espacoDeRotulos2)
          groups_label_space = groups_label_space[,c(-1,-3)]
          a = dist(groups_label_space)
          b = as.dist(a)
          sil = silhouette(groups_label_space[,1], b)
          sil = sortSilhouette(sil)
          write.csv(sil, paste(FolderPartVal, "/silho-fold", f, "-knn-", k,
                               "-part-", p, ".csv", sep=""),
                    row.names = FALSE)

          if(all(is.na(sil))==TRUE){
            #cat("\nOne label per group (local partition)\n")
            fold = f
            part = p
            knn = k
            maximo = NA
            minimo = NA
            mediana = NA
            media = NA
            primeiroQuadrante = NA
            terceiroQuadrante = NA
            valueSilhouete = NA
            Silhouete = rbind(Silhouete, data.frame(fold, part, knn, maximo,
                                                    minimo, mediana, media,
                                                    primeiroQuadrante,
                                                    terceiroQuadrante,
                                                    valueSilhouete))
            write.csv(Silhouete[-1,], paste(FolderKnnVal, "/fold-", f,
                                            "-knn-", k,
                                            "-silho.csv", sep=""),
                      row.names = FALSE)

          } else {
            #cat("\nMore than one label per group\n")
            pdf(paste(FolderPartVal, "/silho-fold-", f, "-knn-", k,
                      "-part-", p, ".pdf", sep=""), width = 10, height = 8)
            print(plot(sil))
            dev.off()
            cat("\n")

            pdf(paste(FolderPartVal, "/fviz-silh-fold-", f, "-knn-", k,
                      "-part-", p, ".pdf", sep=""), width = 10, height = 8)
            print(fviz_silhouette(sil))
            dev.off()
            cat("\n")

            # Summary of silhouette analysis
            si.sum = summary(sil)
            res.si.sum = unlist(si.sum)

            fold = f
            part = p
            knn = k
            maximo = res.si.sum$si.summary.Max.
            minimo = res.si.sum$si.summary.Min.
            mediana = res.si.sum$si.summary.Median
            media = res.si.sum$si.summary.Mean
            primeiroQuadrante = res.si.sum$`si.summary.1st Qu.`
            terceiroQuadrante = res.si.sum$`si.summary.3rd Qu.`
            valueSilhouete = res.si.sum$avg.width
            Silhouete = rbind(Silhouete, data.frame(fold, part, knn, maximo,
                                                    minimo, mediana, media,
                                                    primeiroQuadrante,
                                                    terceiroQuadrante,
                                                    valueSilhouete))
            write.csv(Silhouete[-1,], paste(FolderKnnVal, "/fold-", f,
                                            "-knn-", k, "-silho.csv", sep=""),
                      row.names = FALSE)
          } # fim do if

        } # fim do if

        p = p + 1
        gc()
      } # fim da partição

      Silhouete = Silhouete[-1,]
      indice = as.numeric(which.max(Silhouete$valueSilhouete))
      silhouete2 = Silhouete[indice,]
      bestPartition = rbind(bestPartition, silhouete2)

      k = k + 1
      gc()

    } # fim do knn

    write.csv(bestPartition[-1,], paste(FolderSplitVal, "/fold-", f,
                                        "-best-silho.csv", sep=""),
              row.names = FALSE)

    #f = f + 1
    gc()

  } # fim do fold

  gc()
  cat("\n#################################################################")
  cat("\n# END COMPUTE SILHOUETE                                         #")
  cat("\n#################################################################")
  cat("\n\n\n\n")
}


######################################################################
#
######################################################################
silho.choosed <- function(parameters){

  retorno = list()

  # generating partitions and groups
  partitions = seq(1,parameters$Dataset.Info$Labels,by=1)
  groups = seq(1,parameters$Dataset.Info$Labels,by=1)
  all.partitions = data.frame(partitions, groups)

  partitions = seq(2,(parameters$Dataset.Info$Labels-1),by=1)
  groups = seq(2,(parameters$Dataset.Info$Labels-1),by=1)
  all.hybrid.partitions = data.frame(partitions, groups)

  # saving information
  setwd(parameters$Folders$folderReports)
  write.csv(all.partitions, "all-partitions.csv", row.names = FALSE)

  # data frames
  todos.knn.h = data.frame()
  todos.eb = data.frame()
  todos.fg = data.frame()
  todos.wt = data.frame()
  todos = data.frame()

  # vectors
  total = c(0)
  nomes = c()

  # for all folds
  f = 1
  while(f<=number_folds){
    cat("\n#=======================================================#")
    cat("\n# FOLD ", f, "                                          #")
    cat("\n#=======================================================#\n")

    # folder name
    #  "/dev/shm/j-GpositiveGO/Partitions/Split-1"
    FolderSplit = paste(parameters$Folders$folderPartitions,
                        "/Split-", f,
                        sep="")

    # file name
    # "/dev/shm/j-GpositiveGO/Partitions/Split-1/fold-1-knn-h-choosed.csv"
    knn.h = paste(FolderSplit, "/fold-", f, "-knn-h-choosed.csv", sep="")

    # open file with all methods chosen for FOLD=1
    a.knn.h = data.frame(read.csv(knn.h))

    # how many methods there are?
    total[f] = nrow(a.knn.h)

    # what is it folds?
    nomes[f] = paste("fold-",f,sep="")

    # save the data frame with that information
    todos.knn.h = rbind(todos.knn.h, a.knn.h)

    # through all knn sparcification
    k = 1
    while(k<=nrow(a.knn.h)){
      cat("\n\tKNN ", k)

      # "/dev/shm/j-GpositiveGO/Partitions/Split-1/knn-1"
      FolderKnn = paste(FolderSplit, "/knn-", k, sep="")

      # file names
      eb = paste(FolderKnn, "/knn-", k, "-eb-partitions-hierarchical.csv",
                 sep="")

      fg = paste(FolderKnn, "/knn-", k, "-fg-partitions-hierarchical.csv",
                 sep="")

      wt = paste(FolderKnn, "/knn-", k, "-wt-partitions-hierarchical.csv",
                 sep="")

      # open files
      a.eb = data.frame(read.csv(eb))
      a.fg = data.frame(read.csv(fg))
      a.wt = data.frame(read.csv(wt))

      # gather info with the specific fold and knn
      d.eb = data.frame(fold = f, knn = k, a.eb)
      d.fg = data.frame(fold = f, knn = k, a.fg)
      d.wt = data.frame(fold = f, knn = k, a.wt)

      todos.eb = rbind(todos.eb, d.eb)
      todos.fg = rbind(todos.fg, d.fg)
      todos.wt = rbind(todos.wt, d.wt)

      m.eb = data.frame(fold = f, knn = k, method = "eb", a.eb)
      m.fg = data.frame(fold = f, knn = k, method = "fg", a.fg)
      m.wt = data.frame(fold = f, knn = k, method = "wt", a.wt)

      res = rbind(m.eb, m.fg, m.wt)
      todos = rbind(todos, res)

      k = k + 1
      gc()
    } # fim do knn

    f = f + 1
    gc()
  } # fim do fold

  setwd(parameters$Folders$folderReports)

  # saving sparcification
  sparcification = data.frame(nomes, total)
  write.csv(sparcification, "sparcification.csv", row.names = FALSE)

  # saving all information
  write.csv(todos, "all-partitions.choosed.csv", row.names = FALSE)
  write.csv(todos.knn.h, "all-knn-h-choosed.csv", row.names = FALSE)
  write.csv(todos.eb, "all-eb-partitions.csv", row.names = FALSE)
  write.csv(todos.fg, "all-fg-partitions.csv", row.names = FALSE)
  write.csv(todos.wt, "all-wt-partitions.csv", row.names = FALSE)

  # return
  retorno$all.partitions.choosed = todos
  retorno$all.methods.choosed = todos.knn.h
  retorno$all.eb.partitions = todos.eb
  retorno$all.fg.partitions = todos.fg
  retorno$all.wt.partitions = todos.wt
  retorno$sparcification = sparcification
  retorno$all.partitions = all.partitions
  retorno$all.hybrid.partitions = all.hybrid.partitions
  return(retorno)

  cat("\n\n##########################################################")
  cat("\n# FINISH CHOOSED                                           #")
  cat("\n############################################################\n\n")

}


######################################################################
#
######################################################################
silho.best.partitions <- function(parameters){

  retorno = list()

  all.silho = data.frame()
  sparcification = data.frame()

  cat("\n#=========================================================")
  f = 1
  while(f<=parameters$Number.Folds){

    cat("\n# FOLD = \t", f)

    FolderSplit = paste(parameters$Folders$folderValSilho,
                        "/Split-", f, sep="")

    # fold-1-best-silho.csv
    nome = paste(FolderSplit, "/fold-", f, "-best-silho.csv", sep="")

    resultado = data.frame(read.csv(nome))
    all.silho = rbind(all.silho, resultado)

    # fold-1-sparcification.csv
    nome2 = paste(parameters$Folders$folderValSilho
                  ,"/fold-", f, "-sparcification.csv", sep="")

    res.2 = data.frame(read.csv(nome2))
    sparcification = rbind(sparcification, res.2)

    unlink(nome2)

    f = f + 1
    gc()
  }
  cat("\n#=========================================================")

  retorno$all.silhouette = all.silho
  retorno$sparcification = sparcification

  setwd(parameters$Folders$folderReports)
  write.csv(all.silho, "best-silhouete.csv", row.names = FALSE)

  setwd(parameters$Folders$folderReports)
  write.csv(sparcification, "silho-sparcification.csv", row.names = FALSE)

  n = mean(sparcification$sparc)

  res.4 = data.frame()
  k = 1
  while(k<=n){
    res.1 = data.frame(filter(all.silho, knn == k))
    res.2 = data.frame(count(res.1, part))
    names(res.2) = c("partitions", "frequency")
    res.3 = cbind(knn = k, res.2)
    res.4 = rbind(res.4, res.3)
    k = k + 1
    gc()
  }

  setwd(parameters$Folders$folderReports)
  write.csv(res.4, "silho-frequency-partitions-knn.csv",
            row.names = FALSE)

  return(retorno)

  gc()
  cat("\n#################################################################")
  cat("\n# END BEST PARTITIONS SILHOUETE                                 #")
  cat("\n#################################################################")
  cat("\n\n\n\n")

}


#####################################################################
# Please, any errors, contact us: elainececiliagatto@gmail.com      #
# Thank you very much!                                              #
#####################################################################

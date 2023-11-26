get.cell.type.sig <- function(r = NULL, subtype.flag = F) {
  if (is.null(r)) {
    r <- readRDS("../Data/scData/Mel.all.data.QC.rds")
    cell.type.de <- get.cell.type.sig(r)
    return(cell.type.de)
  }
  cell.types.u <- sort(unique(r$cell.types))
  # Make all pairwise comparisons
  cell.type.pairs <- t(combn(unique(r$cell.types), 2))
  de.zscores <- apply(cell.type.pairs, 1, function(x) {
    print(paste("Comparing", x[1], "to", x[2]))
    b1 <- is.element(r$cell.types, x[1])
    b2 <- is.element(r$cell.types, x[2])
    de <- apply.ttest(r$tpm[, b1 | b2], b1[b1 | b2], ranksum.flag = T)[, "zscores"]
    return(de)
  })
  rownames(de.zscores) <- r$genes
  colnames(de.zscores) <- paste(cell.type.pairs[, 1], cell.type.pairs[, 2], sep = "_")
  results <- list(cell.type.pairs = cell.type.pairs, cell.types = cell.types.u, de.zscores = de.zscores)
  results$de.sum <- lapply(cell.types.u, function(x) {
    b1 <- is.element(cell.type.pairs[, 1], x)
    b2 <- is.element(cell.type.pairs[, 2], x)
    de <- cbind.data.frame(de.zscores[, b1], -de.zscores[, b2])
    colnames(de) <- c(cell.type.pairs[b1, 2], cell.type.pairs[b2, 1])
    de$Min <- rowMins(as.matrix(de))
    return(de)
  })
  names(results$de.sum) <- cell.types.u
  results$sig <- lapply(results$de.sum, function(x) sort(rownames(x)[x$Min > 5]))
  if (!subtype.flag) {
    results <- get.t.cell.sig(results)
  }
  results$sig <- remove.ribo(results$sig)

  cell.type.de <- results
  cell.sig <- results$sig
  if (!subtype.flag) {
    save(cell.type.de, file = "../Results/CellTypes/cell.type.de.RData")
    save(cell.sig, file = "../Results/CellTypes/cell.type.sig.RData")
    save(cell.sig, file = "../Results/Signatures/cell.type.sig.RData")
    cell.type.de <- get.cell.supertype.sig(de = cell.type.de)
    cell.sig <- get.all.t.subset.sig(rF = r, r4 = NULL, r8 = NULL)
    load("../Data/PublicData/mouse_human_mapping.RData")
    cell.sig.mouse <- lapply(cell.sig, function(x) convert.to.mice(mouse.human = mouse.human, x))
    save(cell.sig.mouse, file = "../Results/Signatures/cell.type.sig.full.mouse.RData")
    print("Writing TableS4: Denovo cell sub/type signatures.")
    write.csv(list.2.mat(cell.sig),
      row.names = F,
      file = "../Output/Tables/TableS4_denovo.cell.subtype.sig.csv"
    )
    type.OE <- compute.cell.type.OE(r = r, cell.sig = cell.sig)
  }
  return(cell.type.de)
}

get.t.cell.sig <- function(results) {
  print("Generating CD4 and CD8 T cell signatures.")
  f <- function(m, cell.types) {
    cell.types <- intersect(colnames(m), cell.types)
    m$Min_noT <- rowMins(as.matrix(m[, setdiff(cell.types, c("T.CD4", "T.CD8", "NK"))]))
    m$Min_noNK <- rowMins(as.matrix(m[, setdiff(cell.types, "NK")]))
    return(m)
  }
  results$de.sum$T.CD4 <- f(results$de.sum$T.CD4, results$cell.types)
  results$de.sum$T.CD8 <- f(results$de.sum$T.CD8, results$cell.types)
  genes <- rownames(results$de.zscores)
  results$sig$T.cell <- sort(genes[results$de.sum$T.CD4$Min_noT > 5 & results$de.sum$T.CD8$Min_noT > 5])
  results$sig$T.CD8.NK <- sort(genes[results$de.sum$T.CD8$Min_noNK > 5])
  results$sig$T.CD8 <- get.top.elements(-results$de.sum$T.CD8, no.elm = 50)[["Min"]]
  results$sig$T.CD4 <- get.top.elements(-results$de.sum$T.CD4, no.elm = 50)[["Min"]]
  return(results)
}

get.cell.supertype.sig <- function(de) {
  print("Generating supertype signatures (stroma cells, lymphocytes, and immune cells).")
  f <- function(de, cell.types, supertype) {
    other.types <- setdiff(de$cell.types, cell.types)
    m <- get.mat(rownames(de$de.zscores), m.cols = cell.types)
    for (i in cell.types) {
      m[, i] <- rowMins(as.matrix(de$de.sum[[i]][, other.types]))
    }
    m <- cbind.data.frame(m, min = rowMins(as.matrix(m)))
    de$de.sum[[supertype]] <- m
    de$sig[[supertype]] <- get.top.elements(-m, no.elm = 100, min.cf = (-3))$min
    return(de)
  }
  supertypes <- list(
    stroma = c("CAF", "Endo."),
    lymphocyte = c("B.cell", "T.CD4", "T.CD8"),
    immune = c("B.cell", "Macrophage", "NK", "T.CD4", "T.CD8")
  )
  for (x in names(supertypes)) {
    de <- f(de = de, cell.types = supertypes[[x]], supertype = x)
  }
  cell.type.de <- de
  cell.sig <- de$sig
  save(cell.type.de, file = "../Results/CellTypes/cell.type.de.full.RData")
  save(cell.sig, file = "../Results/CellTypes/cell.type.sig.full.RData")
  save(cell.sig, file = "../Results/Signatures/cell.type.sig.full.RData")
  return(de)
}

get.all.t.subset.sig <- function(rF = NULL, r4 = NULL, r8 = NULL) {
  print("Generating T cell state signatures.")
  t.markers <- readRDS("../Data/PublicData/T.cell.state.markers.rds")
  if (is.null(rF)) {
    rF <- readRDS("../Data/scData/Mel.all.data.QC.rds")
  }
  if (is.null(r4)) {
    r4 <- readRDS("../Data/scData/Mel.T.CD4.QC.rds")
  }
  if (is.null(r8)) {
    r8 <- readRDS("../Data/scData/Mel.T.CD8.QC.rds")
  }

  de4 <- get.t.cell.subset.sig(rF, r4, t.markers, T.type = "T.CD4")
  de8 <- get.t.cell.subset.sig(rF, r8, t.markers, T.type = "T.CD8")
  t.sig4 <- readRDS("../Results/CellTypes/T.CD4.subtype.sig.rds")
  t.sig8 <- readRDS("../Results/CellTypes/T.CD8.subtype.sig.rds")
  load("../Results/Signatures/cell.type.sig.full.RData")
  t.sig <- c(cell.sig[c("T.CD4", "T.CD8", "T.cell")], t.sig4, t.sig8)
  cell.sig <- c(cell.sig[setdiff(names(cell.sig), c(names(t.sig4), names(t.sig8)))], t.sig4, t.sig8)
  saveRDS(t.sig, file = "../Results/CellTypes/T.subtype.sig.rds")
  save(cell.sig, file = "../Results/CellTypes/cell.type.sig.full.RData")
  save(cell.sig, file = "../Results/Signatures/cell.type.sig.full.RData")
  return(cell.sig)
}

get.t.cell.subset.sig <- function(rF, r, t.markers, T.type = "T.CD8") {
  print(paste("Generating", T.type, "state signatures."))
  if (is.null(r)) {
    t.markers <- readRDS("../Data/PublicData/T.cell.state.markers.rds")
    r <- readRDS(paste0("../Data/scData/Mel.", T.type, ".QC.rds"))
    rF <- readRDS("../Data/scData/Mel.all.data.QC.rds")
    de <- get.t.cell.subset.sig(rF, r, t.markers, T.type = ifelse(grepl("CD4", r$data.name), "T.CD4", "T.CD8"))
  }

  print(paste("1. Assign", T.type, "to cell subtypes according to well-established markers."))
  r$scores <- get.OE.sc(r, t.markers)
  B <- discretize.mat.q(r$scores, q1 = 0.9)
  L <- get.mat(m.rows = r$cells, m.cols = names(t.markers), data = "")
  for (i in names(t.markers)) {
    L[, i] <- ifelse(B[, i], i, "")
  }
  r$cell.types <- remove.redundant.dashs(paste(L[, 1], L[, 2], L[, 3], sep = "_"))
  if (T.type == "T.CD4") {
    r <- find.Tregs(r)
    B <- cbind.data.frame(B, Treg = r$B.treg$Treg)
    r$cell.types[r$B.treg$Treg] <- paste0(r$cell.types[r$B.treg$Treg], "_Treg")
    r$cell.types <- remove.redundant.dashs(r$cell.types)
  }
  r$cell.types[r$cell.types == ""] <- "other"
  r1 <- set.list(r, rowSums(B) < 2)
  print("Cell subtypes:")
  print(table(r1$cell.types))

  print("2. Get subtype signatures.")
  de <- get.cell.type.sig(r1, subtype.flag = T)
  de$sig <- remove.ribo(lapply(de$de.sum, function(x) sort(rownames(x)[x$Min > 2])))

  print("3. Ensure the genes in the subtype signatures are underexpressed in other cell types.")
  u.cell.types <- sort(setdiff(unique(r1$cell.types), "other"))
  de$other.cell.types <- list()
  de$sig.specific <- list()
  cell.types.no.nk <- setdiff(sort(unique(rF$cell.types)), c(T.type, "NK"))
  for (x in u.cell.types) {
    genes <- de$sig[[x]]
    b <- is.element(rF$cells, r$cells[r$cell.types == x]) & rF$cell.types == T.type
    m <- cmp.multi.refs(rF$tpm[genes, ], b, ref.groups = rF$cell.types)$up
    m <- cbind.data.frame(max.p.no.nk = rowMaxs(as.matrix(m[, cell.types.no.nk])), m)
    de$sig.specific[[x]] <- rownames(m)[m$max.p < 0.05]
    de$sig.specific.no.nk[[x]] <- rownames(m)[m$max.p.no.nk < 0.05]
    de$other.cell.types[[x]] <- m
  }
  t.sig <- de$sig.specific
  names(t.sig) <- paste0(T.type, ".", names(t.sig))
  t.sig <- t.sig[laply(t.sig, length) > 5]

  print("Saving the results...")
  saveRDS(r, file = paste0("../Data/scData/Mel.", T.type, ".QC.rds"))
  saveRDS(de, file = paste0("../Results/CellTypes/", T.type, ".subtype.de.rds"))
  saveRDS(t.sig, file = paste0("../Results/CellTypes/", T.type, ".subtype.sig.rds"))
  saveRDS(t.sig, file = paste0("../Results/Signatures/", T.type, ".subtype.sig.rds"))
  return(de)
}

find.Tregs <- function(r) {
  if (is.null(r)) {
    r <- readRDS("../Data/scData/Mel.T.CD4.QC.rds")
  }
  Treg.genes <- c("FOXP3", "IL2RA") # IL2RA = CD25
  if (is.null(r$imputed)) {
    D1 <- get.seurat.dataset(r$genes, round(r$tpm, 2), D1.name = r$data.name)
    D1 <- AddImputedScore(D1, genes.fit = Treg.genes)
    if (!identical(r$cells, D1@cell.names)) {
      return()
    }
    r$imputed <- as.matrix(D1@imputed)
    saveRDS(r, file = "../Data/scData/Mel.T.CD4.QC.rds")
  }

  B.treg <- get.mat(r$cells, Treg.genes)
  for (x in Treg.genes) {
    par(mfrow = c(2, 2))
    b1 <- r$tpm[x, ] > quantile(r$tpm[x, ], 0.75, na.rm = T)
    b2 <- r$imputed[x, ] > quantile(r$imputed[x, ], 0.75, na.rm = T)
    B.treg[, x] <- b1 & b2
  }
  r$B.treg <- cbind.data.frame(Treg = rowSums(B.treg) == 2, B.treg)
  return(r)
}

compute.cell.type.OE <- function(r, cell.sig) {
  type.OE <- get.OE.sc(r, cell.sig)
  saveRDS(type.OE, file = "../Results/CellTypes/cell.type.OE.rds")
  generate.fig1.S1FG(r, type.OE)
  return(type.OE)
}

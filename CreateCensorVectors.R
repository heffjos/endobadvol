
tmp <- matrix(data=0, nrow=120, ncol=1)
tmp[42, 1] <- 1
write.table(tmp, file='./abs13ins20024_04204/func/verbal/run_04/CensorVectors.csv',
  quote=F, row.names=F, col.names=F)

tmp <- matrix(data=0, nrow=120, ncol=3)
tmp[3, 1] <- 1
tmp[4, 2] <- 1
tmp[41, 3] <- 1
write.table(tmp, file='./abs13ins20024_04204/func/verbal/run_03/CensorVectors.csv',
  quote=F, row.names=F, col.names=F)

tmp <- matrix(data=0, nrow=120, ncol=4)
tmp[6, 1] <- 1
tmp[62, 2] <- 1
tmp[86, 3] <- 1
tmp[87, 4] <- 1
write.table(tmp, file='./abs13ins20024_04204/func/verbal/run_02/CensorVectors.csv',
  quote=F, row.names=F, col.names=F)

tmp <- matrix(data=0, nrow=180, ncol=1)
tmp[148, 1] <- 1
write.table(tmp, file='./abs13ins20024_04204/func/visual/run_03/CensorVectors.csv',
  quote=F, row.names=F, col.names=F)

tmp <- matrix(data=0, nrow=180, ncol=3)
tmp[13, 1] <- 1
tmp[51, 2] <- 1
tmp[52, 3] <- 1
write.table(tmp, file='./abs13ins20024_04204/func/visual/run_02/CensorVectors.csv',
  quote=F, row.names=F, col.names=F)

tmp <- matrix(data=0, nrow=180, ncol=3)
tmp[42, 1] <- 1
tmp[43, 2] <- 1
tmp[70, 3] <- 1
write.table(tmp, file='./abs13ins20024_04204/func/visual/run_01/CensorVectors.csv',
  quote=F, row.names=F, col.names=F)



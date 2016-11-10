

beta <- rep(1/ns,ns)  # Dir prior parameters for sources
data <- list("IASP","IGLN","IGLT","IGLY","IPGM","ITKT","IUNC",
"sourcesASP","sourcesGLN","sourcesGLT","sourcesGLY","sourcesPGM","sourcesTKT","sourcesUNC",
"ns","nat","beta","Nisolates","FULL")

parameters <- c("qASP","qGLN","qGLT","qGLY","qPGM","qTKT","qUNC","phi","P","Z")

inits <- function(){list(g0=rgamma(ns,1,1),Z=round(runif(Nisolates,0.5/ns,1+0.5/ns)*ns),
qASP=structure(.Data=rep(1/nat[1],ns*nat[1]),.Dim=c(ns,nat[1])),
qGLN=structure(.Data=rep(1/nat[2],ns*nat[2]),.Dim=c(ns,nat[2])),
qGLT=structure(.Data=rep(1/nat[3],ns*nat[3]),.Dim=c(ns,nat[3])),
qGLY=structure(.Data=rep(1/nat[4],ns*nat[4]),.Dim=c(ns,nat[4])),
qPGM=structure(.Data=rep(1/nat[5],ns*nat[5]),.Dim=c(ns,nat[5])),
qTKT=structure(.Data=rep(1/nat[6],ns*nat[6]),.Dim=c(ns,nat[6])),
qUNC=structure(.Data=rep(1/nat[7],ns*nat[7]),.Dim=c(ns,nat[7]))
)}

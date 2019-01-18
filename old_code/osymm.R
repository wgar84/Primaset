## This function now works with both 3D and 2D; no need to specify which it is
## X is a matrix of k landmarks by m dimensions of one specimen
# "midline", "right", and "left" are vectors indicating the names (appear 
# as row names of X; could also work with indices but not recommended) 
# of the midline, right and left landmarks, respectively.
# The function uses object symmetry to find the symmetry plane following
# Klingenberg et al. 2002 (Evolution 56:1909-1920)
# missing data should be designated with NA
# The output includes the original configuration ($rec.orig), the symmetric consensus
# configuration ($symmconf), and the reflected configuration with missing landmarks
# reconstructed ($rec.ref). All three option are with missing landmarks reconstructed
# based on the other side if available. 
# If both sides are missing for a certain landmark or for a midline landmark,
# that landmark is simply ignored and returned as NA.

# NOTE: the resulting configuration will be re-organized to have all the midline 
# landmarks first, then the right lateral ones, and then the left lateral ones. 
# It is therefore HIGHLY RECOMMENDED that landmarks will be designated by names in X
# (i.e., rownames(X) is a charcter vector rather than NULL)

# an example for an input file and a protocol that implements this function can be downloaded from 
# http://home.uchicago.edu/~annat/
# Please email me with comments and questions annat22@gmail.com
# Last updated Feb 28th 2011

OSymm <- function(X, midline, right, left) {
	ncl <- ncol(X)
	Xr <- cbind(X[,-ncl], -X[,ncl])
	Xow <- Xo <- rbind(X[c(midline, right, left),])
	Xrw <- Xr <- rbind(Xr[c(midline, left, right),])
	rownames(Xrw) <- rownames(Xr) <- rownames(X)
	Xo[which(is.na(Xr))] <- NA
	Xr[which(is.na(Xo))] <- NA
	mo <- matrix(apply(Xo, 2, mean, na.rm=TRUE), byrow=TRUE, nr=nrow(Xow), nc=ncol(Xow))
	mr <- matrix(apply(Xr, 2, mean, na.rm=TRUE), byrow=TRUE, nr=nrow(Xrw), nc=ncol(Xrw))
	Xrwc <- Xrw-mr
	SVD <- svd(t(na.omit(Xr-mr)) %*% na.omit(Xo-mo))
	L <- diag(SVD$d)
	S <- ifelse(L<0, -1, L)
	S <- ifelse(L>0, 1, L)
	RM <- SVD$v %*% S %*% t(SVD$u)
	Xrot <- (Xow-mo) %*% RM                                                                                      
	SC <- apply(array(c(Xrwc,Xrot),
                          dim=c(nrow(Xrot),ncol(Xrot),2),
                          dimnames=list(rownames(Xrot),colnames(Xrot))), 1:2, mean, na.rm=TRUE)
	Xrot[which(is.na(Xrot))] <- Xrwc[which(is.na(Xrot))]
	Xrwc[which(is.na(Xrwc))] <- Xrot[which(is.na(Xrwc))]
	list(rec.orig=Xrot, symmconf=SC, rec.ref=Xrwc)}

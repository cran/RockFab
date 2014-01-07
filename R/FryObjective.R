FryObjective <-
function(object.data, n.pass = 15, pie.step = 12, expansion = 1.5, pie.pts = 1, section.name, norm = TRUE){ 
  #Function to fit 2D ellipse to points
  EllipFit <- function(dat, node.start = 0, node.stop = 2 * pi, npts = 180){
    D1 <- cbind(dat$x * dat$x, dat$x * dat$y, dat$y * dat$y)
    D2 <- cbind(dat$x, dat$y, 1)
    S1 <- t(D1) %*% D1
    S2 <- t(D1) %*% D2
    S3 <- t(D2) %*% D2
    T <- -solve(S3) %*% t(S2)
    M <- S1 + S2 %*% T
    M <- rbind(M[3,] / 2, -M[2,], M[1,] / 2)
    evec <- eigen(M)$vec
    cond <- 4 * evec[1,] * evec[3,] - evec[2,]^2
    a1 <- evec[, which(cond > 0)]
    f <- c(a1, T %*% a1)
    names(f) <- letters[1:6]
    
    #Calculate the center and lengths of the semi-axes
    b2 <- f[2]^2 / 4
    center <- c((f[3] * f[4] / 2 - b2 * f[5]), (f[1] * f[5] / 2 - f[2] *
    f[4] / 4)) / (b2 - f[1] * f[3])
    names(center) <- c("x", "y")
    
    num <- 2 * (f[1] * f[5]^2 / 4 + f[3] * f[4]^2 / 4 + f[6] * b2 -
    f[2]*f[4]*f[5]/4 - f[1]*f[3]*f[6])
    den1 <- (b2 - f[1]*f[3])
    den2 <- sqrt((f[1] - f[3])^2 + 4*b2)
    den3 <- f[1] + f[3]
    
    semi.axes <- sqrt(c( num / (den1 * (den2 - den3)),  num / (den1 *
      (-den2 - den3)) ))
  
    #calculate the angle of rotation
    term <- (f[1] - f[3]) / f[2]
    angle <- atan(1 / term) / 2
    
    #Determine XY plotting coords
    tt <- seq(from = node.start, to = node.stop, length = npts)
    sa <- sin(angle)
    ca <- cos(angle)
    ct <- cos(tt)
    st <- sin(tt)
    x <- max(semi.axes) * ct * ca - min(semi.axes) * st * sa
    y <- max(semi.axes) * ct * sa + min(semi.axes) * st * ca
    
    angle <- angle * -1
    angle <- ifelse(angle < 0, angle + pi, angle)
    
    my.ellipse <- list()
    my.ellipse[[1]] <- c(max(semi.axes), min(semi.axes))
    my.ellipse[[2]] <- unname(angle)
    my.ellipse[[3]] <- data.frame(x, y)
    return(my.ellipse)
  }

  #Define empty objects for coordinates
  x.coords <- NULL
  y.coords <- NULL
  radii <- NULL
  
  #Determine object radius based on circle with equal area
  object.data$radius <- sqrt(object.data$m.area / pi)
  
  #Loop through each point to determine Fry cooordinates
  for(j in 1:length(object.data$m.cx)){
    x.coords <- c(x.coords, object.data$m.cx[j] - object.data$m.cx)
y.coords <- c(y.coords, object.data$m.cy[j] - object.data$m.cy)
radii <- c(radii, object.data$radius[j] + object.data$radius)
  }
  
  #Construct data frame of Fry coordinates and remove origin points
  my.data <- data.frame(x.coords, y.coords, radii)
  my.data <- my.data[which(my.data$x.coords != 0 & my.data$y.coords != 0),]
  
  #Determine center to point distances and sort by distance
  my.data$dist <- sqrt(my.data$x.coords^2 + my.data$y.coords^2)
  if(norm){
    norm.dist <- my.data$dist / my.data$radii
    my.data$x.coords <- my.data$x.coords * (norm.dist / my.data$dist) * sqrt(mean(object.data$m.area)/pi)
    my.data$y.coords <- my.data$y.coords * (norm.dist / my.data$dist) * sqrt(mean(object.data$m.area)/pi)
    my.data$dist <- norm.dist * sqrt(mean(object.data$m.area)/pi)
  }
  my.data <- my.data[order(my.data$dist),]
  
  #Determine center to point angle
  my.data$angle <- abs(atan(my.data$y.coords / my.data$x.coords) * (180 / pi))
  my.data$angle[which(my.data$x.coords >= 0 & my.data$y.coords > 0)] <- 270 + my.data$angle[which(my.data$x.coords >= 0 & my.data$y.coords > 0)]
  my.data$angle[which(my.data$x.coords < 0 & my.data$y.coords >= 0)] <- 180 + my.data$angle[which(my.data$x.coords < 0 & my.data$y.coords >= 0)]
  my.data$angle[which(my.data$x.coords <= 0 & my.data$y.coords < 0)] <- 90 + my.data$angle[which(my.data$x.coords <= 0 & my.data$y.coords < 0)]
  
  for(k in 1:n.pass){
x.coord <- vector()
y.coord <- vector()
min.point <- vector()
if(k == 1){
  wedges <- seq(from = 0, to = 360, by = pie.step)
  wedges <- wedges[-length(wedges)]
}
for(j in 1:length(wedges)){
  if(j == 1){
    temp.data <- my.data[which(my.data$angle < wedges[j] | my.data$angle > wedges[length(wedges)]),]
  } else{
    temp.data <- my.data[which(my.data$angle < wedges[j] & my.data$angle > wedges[j - 1]),]
  }
  x.coord <- c(x.coord, temp.data$x.coords[1:(pie.pts)])
  y.coord <- c(y.coord, temp.data$y.coords[1:(pie.pts)])
  min.temp <- temp.data$dist[1]
  min.point <- c(min.point, min.temp)
}
dat <- data.frame(x = x.coord, y = y.coord)
dat <- dat[complete.cases(dat),]

temp.fit <- EllipFit(dat = dat)
dat <- EllipFit(dat = dat, npts = length(wedges) + 1)
node.data <- dat[[3]]
node.data$angles <- abs(atan(node.data$y / node.data$x) * (180 / pi))
node.data$angles[which(node.data$x >= 0 & node.data$y > 0)] <- 270 + node.data$angles[which(node.data$x >= 0 & node.data$y > 0)]
    node.data$angles[which(node.data$x < 0 & node.data$y >= 0)] <- 180 + node.data$angles[which(node.data$x < 0 & node.data$y >= 0)]
    node.data$angles[which(node.data$x <= 0 & node.data$y < 0)] <- 90 + node.data$angles[which(node.data$x <= 0 & node.data$y < 0)]

    wedges <- sort(node.data$angles)
wedges <- wedges[-length(wedges)]
  }
  
  
  fry.data <- new("FRY")
  fry.data@rsAxes <- c(2 * dat[[1]][1], 2 * dat[[1]][2])
  fry.data@strainRatio <- dat[[1]][1] / dat[[1]][2]
  fry.data@meanObjectArea <- dat[[1]][1] * dat[[1]][2] * pi
  fry.data@vectorMean <-  dat[[2]] * (180 / pi)
  fry.data@sectionName <- section.name
  fry.data@sampleSize <- length(object.data$m.cx)
  fry.data@fryParams <- my.data
  fry.data@voidScale <- expansion * max(min.point, na.rm = TRUE)
  
  return(fry.data)
}

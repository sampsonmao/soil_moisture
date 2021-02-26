filename.data <- "SMAP_week.txt"
mydata <- read.table(filename.data)

utmcoor <- SpatialPoints(cbind(mydata$V1,mydata$V2), 
                         proj4string = CRS("+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs"))

lonlatcoor <- spTransform(utmcoor,
                          CRS("+proj=longlat"))

# replace the first two columns with the long lat coordinates
mydata[, 1:2] <- data.frame(lonlatcoor)